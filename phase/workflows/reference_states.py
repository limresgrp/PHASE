"""Reference-first ingestion shared by HTTP and offline clients."""
from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import tempfile
import uuid
from urllib.error import HTTPError
from urllib.request import urlopen

import MDAnalysis as mda
import numpy as np

from phase.io.alignment import align_pdb_to_reference, alignment_quality, build_reference_residue_map, ResidueRef
from phase.services.project_store import DescriptorState, _write_json


def _hash(path):
    with Path(path).open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest() if hasattr(hashlib, "file_digest") else hashlib.sha256(source.read()).hexdigest()


def _residue_identity(residue):
    chains = set(str(value).strip() for value in getattr(residue.atoms, "chainIDs", [""]))
    if len(chains) != 1:
        raise ValueError(f"Residue {residue.resid}: inconsistent chain IDs.")
    return str(ResidueRef(next(iter(chains)), int(residue.resid), str(getattr(residue, "icode", "")).strip()))


def selection_mapping(topology, identity_mapping):
    """Resolve complete chain/resid/insertion identities to exact topology indices."""
    universe = mda.Universe(str(topology))
    try:
        identities = {}
        for residue in universe.residues:
            identities.setdefault(_residue_identity(residue), []).append(int(residue.ix))
        output = {}
        for key, identity in identity_mapping.items():
            if identity is None:
                continue
            indices = identities.get(identity, [])
            if len(indices) != 1:
                raise ValueError(f"{key} -> {identity}: expected one topology residue, found {len(indices)}.")
            output[key] = f"resindex {indices[0]}"
        return output
    finally:
        universe.trajectory.close()


def set_reference(store, project_id, system_id, pdb_path, *, chain_id=None, selection="protein", provenance=None):
    system = store.get_system(project_id, system_id)
    if system.states or system.metastable_clusters:
        raise ValueError("Cannot change the reference/selection after states or clusters exist. Create a new system.")
    reference = build_reference_residue_map(Path(pdb_path), chain_id=chain_id)
    records = reference.records
    mapping = selection_mapping(pdb_path, {r["residue_key"]: str(ref) for r, ref in zip(records, reference.residues)})
    universe = mda.Universe(str(pdb_path))
    try:
        try:
            chosen_residues = universe.select_atoms(f"protein and ({selection or 'protein'})").residues
        except Exception as exc:
            raise ValueError(f"Invalid reference residue selection: {exc}") from exc
        all_mapped_indices = {int(value.split()[1]) for value in mapping.values()}
        unavailable = [f"{_residue_identity(r)} ({r.resname})" for r in chosen_residues
                       if reference.chain_id in {str(c).strip() for c in getattr(r.atoms, "chainIDs", [""])}
                       and int(r.ix) not in all_mapped_indices]
        if unavailable:
            raise ValueError("Selected reference residues have no usable sequence/C-alpha identity "
                             "(unrecognized residue name or missing CA): " + ", ".join(unavailable))
        chosen = set(int(i) for i in chosen_residues.ix)
    finally:
        universe.trajectory.close()
    selected = [r["residue_key"] for r in records if int(mapping[r["residue_key"]].split()[1]) in chosen]
    if not selected:
        raise ValueError("Reference selection matched no residues in the selected protein chain.")
    base = store.resolve_path(project_id, system_id, ".")
    directory = base / "reference"
    directory.mkdir(exist_ok=True)
    destination = directory / "reference.pdb"
    if Path(pdb_path).suffix.lower() not in {".pdb", ".ent"}:
        raise ValueError("Reference ingestion currently requires PDB format. Convert mmCIF to PDB first.")
    metadata = {"schema_version": 1, "chain_id": reference.chain_id, "sequence": reference.sequence,
                "residue_records": records, "selected_residue_keys": selected, "selection": selection,
                "sha256": _hash(pdb_path), "sequence_source": "coordinate residues with CA, not SEQRES",
                "source": provenance or {"source": "upload", "original_filename": Path(pdb_path).name}}
    if Path(pdb_path).resolve() != destination.resolve():
        shutil.copy2(pdb_path, destination)
    _write_json(directory / "reference_metadata.json", metadata)
    system.reference_pdb_file = "reference/reference.pdb"
    system.reference_alignment_metadata_file = "reference/reference_metadata.json"
    system.reference_chain_id = reference.chain_id
    system.reference_sequence = reference.sequence
    system.reference_residue_records = records
    system.reference_residue_keys = selected
    system.reference_selection = selection or "protein"
    system.descriptor_keys = list(selected)
    store.save_system(system)
    return system


def preview_state_alignment(store, project_id, system_id, pdb_path, *, chain_id=None, provenance=None):
    system = store.get_system(project_id, system_id)
    if not system.reference_pdb_file:
        raise ValueError("Set a system reference structure first.")
    reference = store.resolve_path(project_id, system_id, system.reference_pdb_file)
    result = align_pdb_to_reference(reference, Path(pdb_path), reference_chain=system.reference_chain_id, target_chain=chain_id)
    metadata = result.metadata
    selected = system.reference_residue_keys
    metadata.update(selected_reference_residues=list(selected),
                    mapped_selected_residues=[key for key in selected if metadata["residue_mapping"].get(key) is not None],
                    missing_selected_residues=[key for key in selected if metadata["residue_mapping"].get(key) is None],
                    reference_sha256=_hash(reference), target_sha256=_hash(pdb_path),
                    source=provenance or {"source": "upload", "original_filename": Path(pdb_path).name})
    metadata["missing_selected_residue_count"] = len(metadata["missing_selected_residues"])
    metadata["mismatched_selected_residues"] = [row["reference_residue"] for row in metadata["correspondence"]
                                               if row["status"] == "mismatch" and row["reference_residue"] in selected]
    metadata["quality"] = alignment_quality(metadata)
    metadata["alignment_quality_score"] = metadata["identity"] * metadata["coverage_reference"]
    metadata["alignment_status"] = "partial" if metadata["missing_selected_residues"] else "aligned"
    if metadata["alternate_optimal_mapping"]:
        metadata["quality"] = "warning"
    metadata["correspondence"] = sorted(metadata["correspondence"], key=lambda row: next(i for i, record in enumerate(result.reference.records) if record["residue_key"] == row["reference_residue"]))
    # Check that the descriptor reader can resolve the proposed mapping before publication.
    selection_mapping(pdb_path, {key: metadata["residue_mapping"][key] for key in selected})
    return metadata


def load_alignment(store, system, state):
    alignment = state.reference_alignment
    if not alignment:
        return None
    sidecar = alignment.get("metadata_file")
    if sidecar:
        import json
        with store.resolve_path(system.project_id, system.system_id, sidecar).open() as source:
            return json.load(source)
    return alignment


def canonical_mapping(store, system, state):
    alignment = load_alignment(store, system, state)
    if alignment:
        reference = store.resolve_path(system.project_id, system.system_id, system.reference_pdb_file)
        topology = store.resolve_path(system.project_id, system.system_id, state.pdb_file)
        if _hash(reference) != alignment["reference_sha256"] or _hash(topology) != alignment["target_sha256"]:
            raise ValueError("Reference or target structure changed since alignment; re-import the state.")
        identities = {key: alignment["residue_mapping"][key] for key in system.reference_residue_keys}
    elif state.alignment_status == "reference":
        metadata_path = store.resolve_path(system.project_id, system.system_id, system.reference_alignment_metadata_file)
        if _hash(store.resolve_path(system.project_id, system.system_id, state.pdb_file)) != json.loads(metadata_path.read_text())["sha256"]:
            raise ValueError("Reference structure changed on disk since initialization; create a new system.")
        identities = {record["residue_key"]: str(ResidueRef(record["chain_id"], record["resid"], record["insertion_code"]))
                      for record in system.reference_residue_records if record["residue_key"] in system.reference_residue_keys}
        topology = store.resolve_path(system.project_id, system.system_id, state.pdb_file)
    else:
        raise ValueError("State has no reference correspondence; import it through reference-aware ingestion.")
    mapping = selection_mapping(topology, identities)
    return {key: mapping[key] for key in system.reference_residue_keys if key in mapping}


def validate_trajectory(topology, trajectory):
    if Path(trajectory).suffix.lower() in {".pdb", ".ent", ".gro"}:
        reference, target = mda.Universe(str(topology)), mda.Universe(str(trajectory))
        try:
            if len(reference.atoms) != len(target.atoms) or any(not np.array_equal(getattr(reference.atoms, key), getattr(target.atoms, key))
                    for key in ("names", "resids", "resnames")):
                raise ValueError("Trajectory topology/atom ordering does not match the state's PDB. "
                                 "Import an independent PDB as a new aligned state instead of attaching it "
                                 "to an existing topology.")
        finally:
            reference.trajectory.close(); target.trajectory.close()
    universe = mda.Universe(str(topology), str(trajectory))
    try:
        if not len(universe.trajectory):
            raise ValueError("Trajectory contains no frames.")
        universe.trajectory[0]
        return len(universe.trajectory)
    finally:
        universe.trajectory.close()


def create_reference_state(store, project_id, system_id, *, name, state_id=None, pdb_path=None, source_state_id=None,
                           chain_id=None, trajectory=None, copy_traj=True, slice_spec=None, build_descriptors=False,
                           provenance=None, trajectory_source=None):
    from phase.workflows.macro_states import allocate_state_storage_key, refresh_system_metadata, build_state_descriptors
    system = store.get_system(project_id, system_id)
    if not system.reference_pdb_file:
        raise ValueError("Set a system reference first.")
    state_id = state_id or str(uuid.uuid4())
    if state_id in system.states or not name.strip():
        raise ValueError("A unique state ID and nonempty name are required.")
    source = system.states.get(source_state_id) if source_state_id else None
    if pdb_path and source_state_id:
        raise ValueError("Choose a new PDB or an existing topology, not both.")
    if source_state_id and source is None:
        raise ValueError("Source state does not exist.")
    # A self-contained PDB supplied without an explicit topology is a new
    # structure, not coordinates in the reference's atom order.
    if not pdb_path and source is None and trajectory and Path(trajectory).suffix.lower() in {".pdb", ".ent"}:
        pdb_path = trajectory
        provenance = provenance or {"source": "upload", "original_filename": trajectory_source or Path(trajectory).name}
    metadata = None
    if pdb_path:
        metadata = preview_state_alignment(store, project_id, system_id, pdb_path, chain_id=chain_id, provenance=provenance)
    elif source:
        metadata = deepcopy(load_alignment(store, system, source))
    topology = Path(pdb_path) if pdb_path else store.resolve_path(project_id, system_id, source.pdb_file if source else system.reference_pdb_file)
    if trajectory:
        validate_trajectory(topology, trajectory)
    from phase.common.slice_utils import parse_slice_spec
    slice_value, stride = parse_slice_spec(slice_spec) if slice_spec else (None, 1)
    storage_key = allocate_state_storage_key(system, name, state_id)
    dirs = store.ensure_state_directories(project_id, system_id, state_id, storage_key=storage_key)
    directory, base = dirs["state_dir"], dirs["system_dir"]
    try:
        if pdb_path:
            destination = base / "reference" / "topologies" / f"{_hash(topology)}.pdb"
            destination.parent.mkdir(exist_ok=True)
            if not destination.exists():
                temporary = destination.with_name(f".{destination.name}.{uuid.uuid4().hex}.tmp")
                try:
                    shutil.copy2(topology, temporary)
                    temporary.replace(destination)
                finally:
                    temporary.unlink(missing_ok=True)
            topology = destination
        state = DescriptorState(state_id=state_id, name=name.strip(), storage_key=storage_key,
                                pdb_file=str(topology.relative_to(base)), slice_spec=slice_value, stride=stride,
                                residue_keys=list(system.reference_residue_keys),
                                alignment_status=(metadata["alignment_status"] if metadata else "reference"))
        if metadata:
            path = directory / "alignment_metadata.json"
            _write_json(path, metadata)
            state.reference_alignment = {**metadata, "metadata_file": str(path.relative_to(base))}
        state.residue_mapping = canonical_mapping(store, system, state) if state.alignment_status != "partial" else selection_mapping(topology, {
            key: metadata["residue_mapping"][key] for key in system.reference_residue_keys})
        if trajectory:
            destination = directory / f"trajectory{Path(trajectory).suffix}"
            if copy_traj:
                shutil.copy2(trajectory, destination)
            state.trajectory_file = str(destination.relative_to(base)) if copy_traj else str(Path(trajectory).resolve())
            state.source_traj = trajectory_source or str(trajectory)
        system.states[state_id] = state
        if build_descriptors:
            build_state_descriptors(store, project_id, system, state)
        else:
            refresh_system_metadata(system)
            store.save_system(system)
        return state
    except Exception:
        shutil.rmtree(directory)
        system.states.pop(state_id, None)
        refresh_system_metadata(system)
        store.save_system(system)
        raise


def attach_trajectory(store, project_id, system_id, state_id, trajectory, *, copy_traj=True, build_descriptors=True, slice_spec=None, trajectory_source=None):
    from phase.workflows.macro_states import build_state_descriptors, refresh_system_metadata
    system = store.get_system(project_id, system_id)
    state = system.states[state_id]
    validate_trajectory(store.resolve_path(project_id, system_id, state.pdb_file), trajectory)
    directory = store.ensure_state_directories(project_id, system_id, state_id, storage_key=state.storage_key)["state_dir"]
    if copy_traj:
        destination = directory / f"trajectory_{uuid.uuid4().hex}{Path(trajectory).suffix}"
        shutil.copy2(trajectory, destination)
        state.trajectory_file = str(destination.relative_to(store.resolve_path(project_id, system_id, ".")))
    else:
        state.trajectory_file = str(Path(trajectory).resolve())
    state.source_traj = trajectory_source or str(trajectory)
    from phase.common.slice_utils import parse_slice_spec
    state.slice_spec, state.stride = parse_slice_spec(slice_spec) if slice_spec else (None, 1)
    state.descriptor_file = state.descriptor_metadata_file = None
    state.metastable_labels_file = state.metastable_metadata_file = None
    system.metastable_states = [item for item in system.metastable_states if str(item.get("macro_state_id")) != state_id]
    metastable_dir = directory / "metastable"
    if metastable_dir.exists():
        shutil.rmtree(metastable_dir)
    state.n_frames = 0
    # Persist the invalidation even when a subsequent build fails: old labels
    # must never appear as descriptors of a different trajectory.
    refresh_system_metadata(system)
    store.save_system(system)
    if build_descriptors:
        build_state_descriptors(store, project_id, system, state)
    return system


def download_pdb_structure(pdb_id, destination):
    if not re.fullmatch(r"[0-9][A-Za-z0-9]{3}", pdb_id):
        raise ValueError(f"Malformed PDB ID: {pdb_id!r}.")
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    original = None
    try:
        response = urlopen(url, timeout=30)
    except HTTPError as exc:
        if exc.code != 404:
            raise
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
        response = urlopen(url, timeout=30)
        original = destination.with_suffix(".cif")
    with response:
        content = response.read(50 * 1024 * 1024 + 1)
    if len(content) > 50 * 1024 * 1024:
        raise ValueError("Structure download exceeds 50 MiB.")
    if original:
        from Bio.PDB import MMCIFParser, PDBIO
        original.write_bytes(content)
        writer = PDBIO()
        writer.set_structure(MMCIFParser(QUIET=True).get_structure("download", str(original)))
        try:
            writer.save(str(destination))
        except Exception as exc:
            raise ValueError("Downloaded mmCIF cannot be represented losslessly by the current PDB/descriptor reader (chain/atom/residue limits). Import a PDB-compatible selected-chain structure instead.") from exc
    else:
        destination.write_bytes(content)
    build_reference_residue_map(destination)
    _write_json(destination.with_suffix(".download.json"), {"source": "rcsb", "pdb_id": pdb_id.upper(),
                "download_url": url, "retrieved_at": datetime.now(timezone.utc).isoformat(),
                "original_sha256": hashlib.sha256(content).hexdigest(), "converted_from_mmcif": bool(original)})
    return destination


def import_pdb_ids(store, project_id, system_id, pdb_ids, *, preview=False, chain_ids=None,
                   name_template="PDB_{pdb_id}", build_descriptors=False, expected_hashes=None):
    results = []
    for pdb_id in pdb_ids:
        try:
            with tempfile.TemporaryDirectory(prefix="phase-pdb-") as temp:
                path = download_pdb_structure(str(pdb_id).strip(), Path(temp) / "structure.pdb")
                provenance = {"source": "rcsb", "pdb_id": pdb_id.upper(),
                              "download_url": f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb",
                              "retrieved_at": datetime.now(timezone.utc).isoformat()}
                if path.with_suffix(".download.json").exists():
                    provenance = json.loads(path.with_suffix(".download.json").read_text())
                if not preview and path.with_suffix(".cif").exists():
                    original = path.with_suffix(".cif")
                    saved = store.resolve_path(project_id, system_id, "reference/originals") / f"{_hash(original)}.cif"
                    saved.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(original, saved)
                    provenance["original_structure_file"] = str(saved.relative_to(store.resolve_path(project_id, system_id, ".")))
                metadata = preview_state_alignment(store, project_id, system_id, path,
                    chain_id=(chain_ids or {}).get(pdb_id), provenance=provenance)
                if expected_hashes and expected_hashes.get(pdb_id) != metadata["target_sha256"]:
                    raise ValueError("Downloaded structure changed since preview. Preview again before creating this state.")
                state = None if preview else create_reference_state(store, project_id, system_id,
                    name=name_template.format(pdb_id=pdb_id.upper()), pdb_path=path,
                    chain_id=metadata["target_chain"], provenance=provenance, build_descriptors=build_descriptors)
                results.append({"pdb_id": pdb_id, "status": "preview" if preview else "created",
                                "state_id": state.state_id if state else None, "alignment": metadata})
        except Exception as exc:
            results.append({"pdb_id": pdb_id, "status": "failed", "error": str(exc)})
    return results
