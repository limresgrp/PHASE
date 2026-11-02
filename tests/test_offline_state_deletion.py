import pytest

from phase.scripts.offline_system import delete_state
from phase.services.project_store import DescriptorState, ProjectStore


def test_delete_state_removes_hydration_metadata_and_preserves_same_named_state(tmp_path):
    store = ProjectStore(tmp_path / "projects")
    store.create_project("Project", project_id="p")
    system = store.create_system("p", name="System", system_id="s")
    root = tmp_path / "projects/p/systems/s"
    external = tmp_path / "source.xtc"
    external.write_bytes(b"external original")
    system.states = {
        "uuid-to-delete": DescriptorState(state_id="uuid-to-delete", name="FANG", storage_key="fang",
            pdb_file="states/fang/structure.pdb", trajectory_file=str(external)),
        "FANG": DescriptorState(state_id="FANG", name="FANG", storage_key="fang_2", pdb_file="states/fang_2/structure.pdb"),
    }
    store.save_system(system)
    for folder in ("fang", "fang_2"):
        (root / "states" / folder / "structure.pdb").write_text("END\n")
    delete_state(store, "p", "s", "uuid-to-delete")
    assert not (root / "states/fang").exists()
    assert (root / "states/fang_2/state_metadata.json").is_file()
    assert (root / "states/fang_2/structure.pdb").read_text() == "END\n"
    assert external.read_bytes() == b"external original"
    # Repeated reloads must not rediscover the deleted state.
    for _ in range(2):
        assert set(ProjectStore(tmp_path / "projects").get_system("p", "s").states) == {"FANG"}


def test_delete_metadata_only_ghost_and_its_metastable_entries(tmp_path):
    store = ProjectStore(tmp_path / "projects")
    store.create_project("Project", project_id="p")
    system = store.create_system("p", name="System", system_id="s")
    system.states["ghost"] = DescriptorState(state_id="ghost", name="FANG", storage_key="fang", pdb_file="states/fang/structure.pdb")
    system.metastable_states = [{"metastable_id": "ghost__m0", "macro_state_id": "ghost"}]
    store.save_system(system)
    delete_state(store, "p", "s", "ghost")
    restored = store.get_system("p", "s")
    assert not restored.states
    assert not restored.metastable_states
    assert not restored.analysis_states


def test_delete_legacy_flat_files_does_not_rediscover_state(tmp_path):
    store = ProjectStore(tmp_path / "projects")
    store.create_project("Project", project_id="p")
    system = store.create_system("p", name="System", system_id="s")
    root = tmp_path / "projects/p/systems/s"
    (root / "structures/old.pdb").write_text("END\n")
    system.states["old"] = DescriptorState(state_id="old", name="Old", pdb_file="structures/old.pdb")
    store.save_system(system)
    delete_state(store, "p", "s", "old")
    assert not store.get_system("p", "s").states
    assert not (root / "structures/old.pdb").exists()


def test_delete_reports_storage_permission_errors(tmp_path, monkeypatch):
    store = ProjectStore(tmp_path / "projects")
    store.create_project("Project", project_id="p")
    system = store.create_system("p", name="System", system_id="s")
    system.states["s1"] = DescriptorState(state_id="s1", name="State", storage_key="state")
    store.save_system(system)
    def denied(*args, **kwargs):
        raise PermissionError("State folder is not writable")
    monkeypatch.setattr("phase.scripts.offline_system.shutil.rmtree", denied)
    with pytest.raises(PermissionError, match="not writable"):
        delete_state(store, "p", "s", "s1")
    assert "s1" in store.get_system("p", "s").states
