from pathlib import Path
import shutil

import numpy as np
import pytest

from phase.io.alignment import align_pdb_to_reference, build_reference_residue_map, map_selection_to_target
from phase.services.project_store import ProjectStore
from phase.workflows.reference_states import set_reference, create_reference_state, canonical_mapping, import_pdb_ids


def pdb(path, sequence, *, start=1, chain="A", insertion_codes=None, append=False):
    from Bio.Data.PDBData import protein_letters_1to3
    lines = []
    for i, letter in enumerate(sequence):
        resid, icode = (insertion_codes[i] if insertion_codes else (start + i, ""))
        for j, name in enumerate(["N", "CA", "C", "O", "CB"]):
            serial = i * 5 + j + 1
            x, y, z = float(i * 3 + j * .4), float(j % 2), float(i % 3)
            lines.append(f"ATOM  {serial:5d} {name:^4s} {protein_letters_1to3[letter].upper():3s} {chain}{resid:4d}{icode:1s}   {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           {name[0]:>2s}\n")
    with Path(path).open("a" if append else "w") as out:
        out.writelines(lines + ["TER\n"])
    return Path(path)


def test_identical_shifted_numbering_and_nonlexical_order(tmp_path):
    ref = pdb(tmp_path / "ref.pdb", "ACDEFGHIK", start=8)
    target = pdb(tmp_path / "target.pdb", "ACDEFGHIK", start=108)
    result = align_pdb_to_reference(ref, target)
    assert result.metadata["identity"] == 1
    assert list(result.metadata["residue_mapping"])[:4] == ["res_8", "res_9", "res_10", "res_11"]
    assert result.metadata["residue_mapping"]["res_8"] == "A:108"


@pytest.mark.parametrize("alias,parent", [
    ("CYX", "C"), ("CYM", "C"), ("HIE", "H"), ("HID", "H"),
    ("HIP", "H"), ("HSD", "H"), ("HSE", "H"), ("HSP", "H"),
    ("ASH", "D"), ("GLH", "E"), ("LYN", "K"),
])
def test_force_field_names_keep_coordinates_and_parent_sequence(tmp_path, alias, parent):
    reference = pdb(tmp_path / "canonical.pdb", "A" + parent + "G", start=2, chain="X")
    target = pdb(tmp_path / "md.pdb", "A" + parent + "G", start=2, chain="X")
    target.write_text("".join(
        line[:17] + alias + line[20:]
        if line.startswith("ATOM") and int(line[22:26]) == 3 else line
        for line in target.read_text().splitlines(keepends=True)))
    original = target.read_bytes()
    result = align_pdb_to_reference(reference, target)
    assert result.metadata["identity"] == 1
    assert result.target.sequence == "A" + parent + "G"
    assert result.target.residues[1].resname == alias
    store = ProjectStore(tmp_path / "projects")
    store.create_project("Project", project_id="p")
    store.create_system("p", name="System", system_id="s")
    system = set_reference(store, "p", "s", target, selection="resid 2 to 4")
    assert system.reference_chain_id == "X"
    assert system.reference_residue_keys == ["res_2", "res_3", "res_4"]
    assert target.read_bytes() == original


def test_insertions_deletions_and_selection(tmp_path):
    ref = pdb(tmp_path / "r.pdb", "ACDEFGHIKLMNPQRSTVWY")
    target = pdb(tmp_path / "t.pdb", "ACDEGHIKLMNPQRSTVWY", start=20)
    result = align_pdb_to_reference(ref, target)
    assert result.metadata["residue_mapping"]["res_5"] is None
    assert result.metadata["deletions"] == [4]
    selected = map_selection_to_target([result.reference.residues[i] for i in [1, 4, 5]], result)
    assert list(selected.mapped) == ["res_2", "res_5", "res_6"]
    assert selected.missing == ["res_5"]
    inserted = pdb(tmp_path / "i.pdb", "ACDEWFGHIKLMNPQRSTVWY")
    assert align_pdb_to_reference(ref, inserted).metadata["insertion_count"] == 1


def test_insertion_codes_are_distinct_and_mismatches_quantified(tmp_path):
    ref = pdb(tmp_path / "r.pdb", "ACD", insertion_codes=[(45, ""), (45, "A"), (46, "")])
    assert [r["residue_key"] for r in build_reference_residue_map(ref).records] == ["res_45", "res_45A", "res_46"]
    t = pdb(tmp_path / "t.pdb", "ACE")
    assert align_pdb_to_reference(ref, t).metadata["identity"] == pytest.approx(2 / 3)


def test_best_chain_not_first_and_ambiguity_requires_choice(tmp_path):
    ref = pdb(tmp_path / "r.pdb", "ACDEFGHIK")
    t = pdb(tmp_path / "t.pdb", "WWWWWWWWW")
    pdb(t, "ACDEFGHIK", chain="B", append=True)
    assert align_pdb_to_reference(ref, t).metadata["target_chain"] == "B"
    pdb(t, "ACDEFGHIK", chain="C", append=True)
    with pytest.raises(ValueError, match="Ambiguous target chains"):
        align_pdb_to_reference(ref, t)
    assert align_pdb_to_reference(ref, t, target_chain="C").metadata["identity"] == 1


@pytest.fixture
def reference_system(tmp_path):
    store = ProjectStore(tmp_path / "projects")
    store.create_project("Project", project_id="p")
    store.create_system("p", name="System", system_id="s")
    ref = pdb(tmp_path / "r.pdb", "ACDEFGHIKLMNPQRSTVWY", start=8)
    set_reference(store, "p", "s", ref, selection="resid 9 to 20")
    return store, ref


def test_pdb_only_build_reference_and_shifted_state_same_order(reference_system, tmp_path):
    store, ref = reference_system
    a = create_reference_state(store, "p", "s", name="Reference", build_descriptors=True)
    t = pdb(tmp_path / "t.pdb", "ACDEFGHIKLMNPQRSTVWY", start=108)
    b = create_reference_state(store, "p", "s", name="Shifted", pdb_path=t, build_descriptors=True)
    system = store.get_system("p", "s")
    assert a.residue_keys == b.residue_keys == system.reference_residue_keys
    assert a.n_frames == b.n_frames == 1
    with np.load(store.resolve_path("p", "s", b.descriptor_file)) as data:
        assert data.files == system.reference_residue_keys
    assert b.residue_mapping["res_9"] == "resindex 1"
    assert system.states[b.state_id].reference_alignment["quality"] == "excellent"
    with pytest.raises(ValueError, match="Cannot change"):
        set_reference(store, "p", "s", ref)


def test_partial_state_preserves_keys_and_refuses_build(reference_system, tmp_path):
    store, _ = reference_system
    target = pdb(tmp_path / "t.pdb", "ACDEGHIKLMNPQRSTVWY", start=100)
    b = create_reference_state(store, "p", "s", name="Partial", pdb_path=target)
    system = store.get_system("p", "s")
    assert b.alignment_status == "partial"
    assert b.residue_keys == system.reference_residue_keys
    assert b.reference_alignment["missing_selected_residues"] == ["res_12"]
    with pytest.raises(ValueError, match="res_12"):
        canonical_mapping(store, system, b)
    with pytest.raises(ValueError, match="missing selected"):
        create_reference_state(store, "p", "s", name="Not created", pdb_path=target, build_descriptors=True)
    assert len(store.get_system("p", "s").states) == 1


def test_bad_alignment_never_registers_state(reference_system, tmp_path):
    store, _ = reference_system
    t = pdb(tmp_path / "t.pdb", "WWWW")
    with pytest.raises(ValueError, match="Insufficient"):
        create_reference_state(store, "p", "s", name="Bad", pdb_path=t)
    assert not store.get_system("p", "s").states
    assert not list(store.resolve_path("p", "s", "states").glob("*/state_metadata.json"))


def test_replica_survives_source_deletion(reference_system, tmp_path):
    from phase.scripts.offline_system import delete_state
    store, _ = reference_system
    t = pdb(tmp_path / "t.pdb", "ACDEFGHIKLMNPQRSTVWY", start=108)
    a = create_reference_state(store, "p", "s", name="Source", pdb_path=t)
    b = create_reference_state(store, "p", "s", name="Replica", source_state_id=a.state_id)
    assert a.pdb_file == b.pdb_file
    delete_state(store, "p", "s", a.state_id)
    assert store.resolve_path("p", "s", b.pdb_file).exists()


def test_pdb_batch_failure_is_per_id(reference_system, monkeypatch):
    store, ref = reference_system
    def download(pid, dest):
        if pid == "XXXX":
            raise ValueError("Malformed PDB ID")
        shutil.copy2(ref, dest)
        return dest
    monkeypatch.setattr("phase.workflows.reference_states.download_pdb_structure", download)
    rows = import_pdb_ids(store, "p", "s", ["1ABC", "XXXX"])
    assert [r["status"] for r in rows] == ["created", "failed"]


def test_download_validation_and_404_fallback_are_offline(tmp_path, monkeypatch):
    from urllib.error import HTTPError
    from phase.workflows.reference_states import download_pdb_structure
    with pytest.raises(ValueError, match="Malformed"):
        download_pdb_structure("XXXX", tmp_path / "download.pdb")
    urls = []
    def missing(url, **kwargs):
        urls.append(url)
        raise HTTPError(url, 404, "not found", None, None)
    monkeypatch.setattr("phase.workflows.reference_states.urlopen", missing)
    with pytest.raises(HTTPError):
        download_pdb_structure("1ABC", tmp_path / "download.pdb")
    assert urls == ["https://files.rcsb.org/download/1ABC.pdb", "https://files.rcsb.org/download/1ABC.cif"]


def test_trajectory_replacement_invalidates_descriptors_and_uses_only_traj(reference_system, tmp_path):
    import MDAnalysis as mda
    from phase.workflows.reference_states import attach_trajectory
    store, ref = reference_system
    state = create_reference_state(store, "p", "s", name="Reference", build_descriptors=True)
    universe = mda.Universe(str(ref))
    traj = tmp_path / "test.xtc"
    with mda.Writer(str(traj), universe.atoms.n_atoms) as out:
        for _ in range(3):
            out.write(universe.atoms)
    system = attach_trajectory(store, "p", "s", state.state_id, traj, build_descriptors=False)
    assert system.states[state.state_id].descriptor_file is None
    assert store.get_system("p", "s").states[state.state_id].descriptor_file is None
    from phase.workflows.macro_states import build_state_descriptors
    system = build_state_descriptors(store, "p", system, system.states[state.state_id])
    assert system.states[state.state_id].n_frames == 3  # Not PDB + trajectory.
    assert system.states[state.state_id].residue_keys == system.reference_residue_keys
    replica = create_reference_state(store, "p", "s", name="Replica", source_state_id=state.state_id,
                                     trajectory=traj, build_descriptors=True)
    assert replica.n_frames == 3


def test_descriptor_failure_rolls_back_and_clone_preserves_reference(reference_system, monkeypatch):
    store, _ = reference_system
    from phase.workflows import macro_states
    def fail(*args, **kwargs):
        raise ValueError("Synthetic descriptor failure")
    with monkeypatch.context() as patch:
        patch.setattr(macro_states, "build_state_descriptors", fail)
        with pytest.raises(ValueError, match="Synthetic"):
            create_reference_state(store, "p", "s", name="No leftovers", build_descriptors=True)
    assert not store.get_system("p", "s").states
    assert not list(store.resolve_path("p", "s", "states").glob("*/state_metadata.json"))
    state = create_reference_state(store, "p", "s", name="Reference")
    clone = store.clone_system("p", "s", name="Clone", system_id="clone")
    assert clone.reference_residue_keys == store.get_system("p", "s").reference_residue_keys
    assert canonical_mapping(store, clone, clone.states[state.state_id]) == state.residue_mapping


def test_legacy_shift_still_works_without_reference(tmp_path):
    from phase.workflows.macro_states import add_state
    store = ProjectStore(tmp_path / "projects")
    store.create_project("Project", project_id="p")
    store.create_system("p", name="Legacy", system_id="s")
    topology = pdb(tmp_path / "legacy.pdb", "ACDEFGHIK", start=10)
    state = add_state(store, "p", "s", "legacy", "Legacy", topology, None, "resid 11 to 16", False, True, None, 100)
    assert state.resid_shift == 100
    assert set(state.residue_keys) == {f"res_{i}" for i in range(111, 117)}
    assert store.get_system("p", "s").reference_pdb_file is None


def test_missing_provenance_does_not_hide_system_or_reuse_embedded_mapping(reference_system, tmp_path):
    store, _ = reference_system
    target = pdb(tmp_path / "t.pdb", "ACDEFGHIKLMNPQRSTVWY", start=108)
    state = create_reference_state(store, "p", "s", name="Target", pdb_path=target)
    path = store.resolve_path("p", "s", state.reference_alignment["metadata_file"])
    path.unlink()
    system = store.get_system("p", "s")
    assert system.states[state.state_id].alignment_status == "failed"
    assert "error" in system.states[state.state_id].reference_alignment
    with pytest.raises(FileNotFoundError):
        canonical_mapping(store, system, system.states[state.state_id])


def test_reference_torsions_use_bonded_neighbors_not_residue_number_arithmetic(tmp_path):
    import MDAnalysis as mda
    from phase.features.extraction import FeatureExtractor
    path = pdb(tmp_path / "insertion.pdb", "ACD", insertion_codes=[(45, ""), (45, "A"), (46, "")])
    universe = mda.Universe(str(path))
    for i, residue in enumerate(universe.residues):
        for name, position in [("N", [i * 3., 0, 0]), ("CA", [i * 3. + 1, .5, 0]), ("C", [i * 3. + 2, 0, 0])]:
            residue.atoms.select_atoms(f"name {name}").positions = position
    selections = FeatureExtractor(use_coordinate_neighbors=True)._build_dihedral_selections(universe.residues)
    assert selections["phi"][1].resindices.tolist() == [0, 1, 1, 1]
    assert selections["psi"][0].resindices.tolist() == [0, 0, 0, 1]  # 45 -> 45A, not 46.
    assert selections["omega"][1].resindices.tolist() == [1, 1, 2, 2]
    universe.residues[1].atoms.positions += [20, 0, 0]
    broken = FeatureExtractor(use_coordinate_neighbors=True)._build_dihedral_selections(universe.residues)
    assert broken["phi"][1] is None
    assert broken["psi"][0] is None


def test_clustering_keeps_reference_order_including_insertion_codes():
    from phase.workflows.clustering import _coerce_residue_keys
    from phase.services.project_store import DescriptorState
    keys = ["res_8", "res_9", "res_10", "res_10A", "res_11"]
    state = DescriptorState(state_id="aligned", name="Aligned", alignment_status="aligned", residue_keys=keys)
    assert _coerce_residue_keys([], {key: np.zeros((2, 1, 10)) for key in keys}, state) == keys
    with pytest.raises(ValueError, match="Canonical descriptors missing"):
        _coerce_residue_keys([], {"res_8": np.zeros((2, 1, 10))}, state)


def test_descriptor_frame_count_obeys_actual_python_slice_bounds():
    from types import SimpleNamespace
    from phase.features.extraction import FeatureExtractor
    extractor = FeatureExtractor()
    universe = SimpleNamespace(trajectory=range(5))
    assert extractor._get_sliced_length(universe, slice(0, 100, 2)) == 3
    assert extractor._get_sliced_length(universe, slice(0, 0)) == 0
    assert extractor._get_sliced_length(universe, slice(-2, None)) == 2


def test_offline_cli_reference_create_and_inspection(reference_system, capsys):
    import json
    from phase.scripts.offline_system import main
    store, ref = reference_system
    root = ["--root", str(store.base_dir.parent)]
    assert main(root + ["create-system", "--project-id", "p", "--name", "CLI", "--use-slug-ids",
                        "--reference", str(ref), "--selection", "resid 9 to 20"]) == 0
    capsys.readouterr()
    assert main(root + ["add-state", "--project-id", "p", "--system-id", "cli", "--state-id", "ref",
                        "--use-reference", "--no-build-descriptors"]) == 0
    capsys.readouterr()
    assert main(root + ["alignments", "--project-id", "p", "--system-id", "cli", "--json"]) == 0
    rows = json.loads(capsys.readouterr().out)
    assert rows[0]["status"] == "reference"
    assert rows[0]["state_id"] == "ref"


def test_reference_selection_cannot_silently_drop_a_residue_without_ca(tmp_path):
    store = ProjectStore(tmp_path / "projects")
    store.create_project("Project", project_id="p")
    store.create_system("p", name="System", system_id="s")
    reference = pdb(tmp_path / "ref.pdb", "ACDEFGHIK")
    reference.write_text("".join(line for line in reference.read_text().splitlines(keepends=True)
                                 if not (line.startswith("ATOM") and line[12:16].strip() == "CA" and int(line[22:26]) == 5)))
    with pytest.raises(ValueError, match="no usable sequence/C-alpha identity"):
        set_reference(store, "p", "s", reference)
