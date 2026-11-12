"""Coordinate-derived sequence correspondence; never renumbers or superposes atoms."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB import MMCIFParser, PDBParser
from MDAnalysis.lib.util import inverse_aa_codes


# Coordinate files from MD retain protonation/disulfide names. Normalize only
# their sequence identity, never their stored residue names or coordinates.
_SEQUENCE_LETTERS = {
    **protein_letters_3to1_extended,
    **inverse_aa_codes,
    "CYX": "C", "CYM": "C", "HIP": "H", "HSP": "H",
}


ALIGNMENT_THRESHOLDS = {
    "minimum_coverage": .50, "minimum_identity": .30,
    "excellent_identity": .98, "excellent_coverage": .98,
    "good_identity": .95, "good_coverage": .95, "good_missing": 2,
}


@dataclass(frozen=True)
class ResidueRef:
    chain_id: str
    resid: int
    insertion_code: str = ""
    resname: Optional[str] = None

    def __str__(self):
        return f"{self.chain_id}:{self.resid}{self.insertion_code}"


@dataclass
class ReferenceStructure:
    chain_id: str
    sequence: str
    residues: list[ResidueRef]

    @property
    def records(self):
        return [{**asdict(r), "index": i, "residue_key": f"res_{r.resid}{r.insertion_code}"}
                for i, r in enumerate(self.residues)]


@dataclass
class AlignmentResult:
    reference: ReferenceStructure
    target: ReferenceStructure
    metadata: dict


@dataclass
class SelectionMappingResult:
    mapped: dict[str, Optional[str]]
    missing: list[str]


class AlignmentError(ValueError):
    """Actionable failure; candidates expose choices without registering a state."""
    def __init__(self, message, candidates=None):
        super().__init__(message)
        self.candidates = candidates or []


def protein_chains(path: Path) -> list[ReferenceStructure]:
    parser = MMCIFParser(QUIET=True) if Path(path).suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("structure", str(path))
    except Exception as exc:
        raise AlignmentError(f"Cannot parse protein coordinates in {Path(path).name}: {exc}") from exc
    model = next(iter(structure), None)
    chains = []
    for chain in model or []:
        residues, sequence = [], []
        for residue in chain:
            name = residue.resname.upper().strip()
            if name not in _SEQUENCE_LETTERS or "CA" not in residue:
                continue
            residues.append(ResidueRef(str(chain.id).strip(), int(residue.id[1]), str(residue.id[2]).strip(), name))
            letter = _SEQUENCE_LETTERS[name]
            sequence.append(letter if letter in "ARNDCQEGHILKMFPSTWYVBZX" else "X")
        if residues:
            if len(set(map(str, residues))) != len(residues):
                raise AlignmentError(f"Chain {chain.id}: duplicate residue identities in coordinates.")
            chains.append(ReferenceStructure(str(chain.id).strip(), "".join(sequence), residues))
    if not chains:
        raise AlignmentError("No protein chain with coordinate C-alpha atoms was found.")
    return chains


def build_reference_residue_map(reference_pdb: Path, *, chain_id=None) -> ReferenceStructure:
    chain_id = "" if chain_id == "(blank)" else chain_id
    chains = protein_chains(reference_pdb)
    if chain_id is not None:
        matches = [chain for chain in chains if chain.chain_id == chain_id]
        if not matches:
            raise AlignmentError(f"Protein chain {chain_id!r} not found; available: {[c.chain_id for c in chains]}.")
        return matches[0]
    return max(chains, key=lambda chain: len(chain.residues))


def alignment_quality(metadata):
    t = ALIGNMENT_THRESHOLDS
    if metadata.get("alternate_optimal_mapping"):
        return "warning"
    identity, coverage = metadata["identity"], metadata["coverage_reference"]
    missing = len(metadata.get("missing_selected_residues", []))
    if identity >= t["excellent_identity"] and coverage >= t["excellent_coverage"] and missing == 0:
        return "excellent"
    if identity >= t["good_identity"] and coverage >= t["good_coverage"] and missing <= t["good_missing"]:
        return "good"
    return "warning"


def align_pdb_to_reference(reference_pdb: Path, target_pdb: Path, *, reference_chain=None,
                           target_chain=None, gap_open=-10., gap_extend=-.5) -> AlignmentResult:
    reference = build_reference_residue_map(reference_pdb, chain_id=reference_chain)
    target_chain = "" if target_chain == "(blank)" else target_chain
    targets = protein_chains(target_pdb)
    if target_chain is not None:
        targets = [chain for chain in targets if chain.chain_id == target_chain]
        if not targets:
            raise AlignmentError(f"Target protein chain {target_chain!r} not found.")
    matrix = substitution_matrices.load("BLOSUM62")
    aligner = PairwiseAligner(mode="global", substitution_matrix=matrix,
                              open_gap_score=gap_open, extend_gap_score=gap_extend)
    candidates = sorted([(float(aligner.score(reference.sequence, chain.sequence)), chain) for chain in targets],
                        key=lambda value: value[0], reverse=True)
    if len(candidates) > 1 and np.isclose(candidates[0][0], candidates[1][0]):
        raise AlignmentError("Ambiguous target chains: specify a target chain explicitly. Candidates: "
                             + ", ".join(f"{chain.chain_id or '(blank)'} (score {score:g})" for score, chain in candidates),
                             [{"chain_id": chain.chain_id, "score": score} for score, chain in candidates])
    score, target = candidates[0]
    alignments = aligner.align(reference.sequence, target.sequence)
    aligned = alignments[0]
    mapping = {record["residue_key"]: None for record in reference.records}
    inverse, rows, inserted, deleted = {}, [], [], []
    matches = positives = aligned_count = 0
    for r, t in aligned.indices.T:
        if r < 0:
            inserted.append(int(t)); continue
        key = reference.records[r]["residue_key"]
        if t < 0:
            deleted.append(int(r)); continue
        target_id = str(target.residues[t])
        mapping[key] = target_id
        inverse[target_id] = key
        same = reference.sequence[r] == target.sequence[t]
        matches += int(same)
        positives += int(matrix[reference.sequence[r], target.sequence[t]] > 0)
        aligned_count += 1
        rows.append({"reference_residue": key, "target_residue": target_id,
                     "status": "match" if same else "mismatch"})
    coverage = aligned_count / len(reference.residues)
    identity = matches / max(1, aligned_count)
    if coverage < ALIGNMENT_THRESHOLDS["minimum_coverage"] or identity < ALIGNMENT_THRESHOLDS["minimum_identity"]:
        raise AlignmentError(f"Insufficient sequence correspondence: reference chain {reference.chain_id!r} "
                             f"length {len(reference.residues)}, target {target.chain_id!r} length {len(target.residues)}; "
                             f"coverage {coverage:.1%}, identity {identity:.1%}. Minimum coverage "
                             f"{ALIGNMENT_THRESHOLDS['minimum_coverage']:.0%}, identity {ALIGNMENT_THRESHOLDS['minimum_identity']:.0%}.")
    # Keep alternate optimal mappings explicit rather than pretending the path is unique.
    alternate = False
    try:
        alternate = not np.array_equal(aligned.indices, alignments[1].indices)
    except IndexError:
        pass
    meta = {"schema_version": 1, "method": "global_sequence_alignment", "matrix": "BLOSUM62",
            "quality_thresholds": dict(ALIGNMENT_THRESHOLDS),
            "gap_open": gap_open, "gap_extend": gap_extend, "reference_chain": reference.chain_id,
            "target_chain": target.chain_id, "score": score, "alignment_score": score, "identity": identity,
            "similarity": positives / max(1, aligned_count), "coverage_reference": coverage,
            "coverage_target": aligned_count / len(target.residues), "reference_length": len(reference.residues),
            "target_length": len(target.residues), "matches": matches, "mismatches": aligned_count - matches,
            "insertions": inserted, "deletions": deleted, "insertion_count": len(inserted), "deletion_count": len(deleted),
            "residue_mapping": mapping, "target_to_reference": inverse,
            "target_residue_records": target.records, "alternate_optimal_mapping": alternate,
            "extra_target_residues": [str(target.residues[t]) for t in inserted],
            "correspondence": rows + [{"reference_residue": reference.records[r]["residue_key"],
                                       "target_residue": None, "status": "deletion"} for r in deleted]}
    meta["quality"] = alignment_quality(meta)
    return AlignmentResult(reference, target, meta)


def map_selection_to_target(reference_selection: Sequence[ResidueRef], alignment: AlignmentResult) -> SelectionMappingResult:
    keys = {str(residue): record["residue_key"] for residue, record in zip(alignment.reference.residues, alignment.reference.records)}
    selected = [keys[str(residue)] for residue in reference_selection]
    mapped = {key: alignment.metadata["residue_mapping"][key] for key in selected}
    return SelectionMappingResult(mapped, [key for key, value in mapped.items() if value is None])
