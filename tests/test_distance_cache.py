from dataclasses import replace
import json
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pytest

from phase.potts.metrics import js_divergence, pairwise_joints_on_edges
from phase.scripts.distance_cache import main as distance_cache_main
from phase.services.project_store import DescriptorState, ProjectStore
from phase.workflows.distances import (
    build_distance_cache, distance_bin_indices, load_distance_cache,
    pair_distance_histogram, read_cluster_frames, validate_distance_alignment,
    validate_distance_preset,
)


def _trajectory(path: Path, distances, *, shift=0, periodic=False):
    path.mkdir(parents=True)
    u = mda.Universe.empty(3, n_residues=3, atom_resindex=np.arange(3), trajectory=True)
    u.add_TopologyAttr("names", ["CA"] * 3)
    u.add_TopologyAttr("resnames", ["ALA", "GLY", "SER"])
    u.add_TopologyAttr("resids", np.asarray([10, 20, 30]) + shift)
    u.atoms.positions = [[0, 0, 0], [distances[0], 0, 0], [40, 0, 0]]
    u.atoms.write(str(path / "structure.pdb"))
    with mda.Writer(str(path / "trajectory.xtc"), n_atoms=3) as writer:
        for distance in distances:
            u.atoms.positions = [[0, 0, 0], [distance, 0, 0], [40, 0, 0]]
            if periodic:
                u.dimensions = [50, 50, 50, 90, 90, 90]
            writer.write(u.atoms)


@pytest.fixture
def dataset(tmp_path):
    store = ProjectStore(tmp_path / "projects")
    store.create_project("Toy", project_id="p")
    system = store.create_system("p", name="Toy", system_id="s")
    system_dir = tmp_path / "projects/p/systems/s"
    _trajectory(system_dir / "inputs/a", [2, 4, 7, 10, 12])
    _trajectory(system_dir / "inputs/b", [5, 6, 8, 9, 11], shift=100)
    system.states = {
        "a": DescriptorState(
            state_id="a", name="A", pdb_file="inputs/a/structure.pdb", trajectory_file="inputs/a/trajectory.xtc",
            slice_spec="1:5:2", stride=2, n_frames=2, residue_selection="resid 10 20 30",
            residue_mapping={"res_10": "resid 10", "res_20": "resid 20", "res_30": "resid 30"},
        ),
        "b": DescriptorState(
            state_id="b", name="B", pdb_file="inputs/b/structure.pdb", trajectory_file="inputs/b/trajectory.xtc",
            stride=2, resid_shift=-100, n_frames=3, residue_selection="resid 110 120 130",
        ),
    }
    store.save_system(system)
    cluster_dir = system_dir / "clusters/c"
    cluster_dir.mkdir(parents=True)
    (cluster_dir / "cluster_metadata.json").write_text(json.dumps({
        "cluster_id": "c", "path": "clusters/c/cluster.npz", "max_cluster_frames": 1,
    }))
    cluster_path = cluster_dir / "cluster.npz"
    np.savez_compressed(
        cluster_path, residue_keys=np.asarray(["res_10", "res_20", "res_30"]),
        merged__labels_assigned=np.asarray([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [1, 1, 0], [-1, 0, 0]]),
        merged__frame_state_ids=np.asarray(["b", "a", "a", "b", "a"]),
        merged__frame_metastable_ids=np.asarray(["mb", "ma", "ma", "mb", "overlap"]),
        merged__frame_indices=np.asarray([2, 1, 0, 0, 0]),
        contact_edge_index=np.asarray([[0], [1]]),
    )
    return store, cluster_path


def test_sparse_distances_exact_frame_alignment_slice_shift_and_roundtrip(dataset):
    store, path = dataset
    original = path.read_bytes()
    result = build_distance_cache("p", "s", "c", store=store, max_distance_cutoff=6, periodic=False)
    cache = load_distance_cache(path)
    assert path.read_bytes() == original
    np.testing.assert_array_equal(cache.edge_index, [[0], [1]])
    # Source frames: b[4], a[3], a[1], b[0], a[1]. Fitting subsampling
    # and invalid cluster labels do not remove any geometry rows.
    np.testing.assert_allclose(cache.distances[:, 0], [11, 10, 4, 5, 4], atol=1e-5)
    np.testing.assert_array_equal(cache.frame_indices, [2, 1, 0, 0, 0])
    assert cache.distances.dtype == np.float32
    assert result["size_bytes"] > 0
    assert cache.metadata["schema_version"] == 1
    assert cache.metadata["subsampling"] is None
    assert cache.metadata["stored_contact_edges_count"] == 1
    histogram = pair_distance_histogram(path, 1, 0, state_ids=["a"], metastable_ids=["ma"], thresholds=[5, 6])
    assert histogram["frame_count"] == 2
    assert histogram["bin_counts"] == [1, 0, 1]
    assert sum(histogram["histogram"]["counts"]) == 2
    assert histogram["summary"]["min"] == pytest.approx(4)
    with pytest.raises(ValueError, match="Unknown frame filter"):
        pair_distance_histogram(path, 0, 1, state_ids=["unknown"])
    with pytest.raises(ValueError, match="candidate-edge set"):
        pair_distance_histogram(path, 0, 2)


def test_alignment_cannot_be_checked_by_length_only(dataset):
    store, path = dataset
    build_distance_cache("p", "s", "c", store=store, max_distance_cutoff=6)
    cluster = read_cluster_frames(path)
    cache = load_distance_cache(path)
    for field in ("frame_state_ids", "frame_metastable_ids", "frame_indices", "residue_keys"):
        wrong = replace(cache, **{field: getattr(cache, field)[::-1]})
        with pytest.raises(ValueError, match="alignment mismatch"):
            validate_distance_alignment(cluster, wrong)


def test_boundary_ownership_and_missing_distance_values():
    np.testing.assert_array_equal(distance_bin_indices(np.asarray([5.999, 6, 7.999, 8, np.nan]), [6, 8]), [0, 1, 1, 2, -1])
    for thresholds in ([6, 6], [8, 6], [np.nan], [-1]):
        with pytest.raises(ValueError, match="strictly increasing"):
            distance_bin_indices(np.asarray([5]), thresholds)


def test_legacy_absence_is_one_bin_and_preserves_unconditional_pair_statistics(dataset):
    _, path = dataset
    cache = load_distance_cache(path)
    assert not cache.available
    assert not path.with_name("distance_cache.npz").exists()
    indices = distance_bin_indices(None, [6, 8, 10], frame_count=5)
    assert np.all(indices == 0)
    result = pair_distance_histogram(path, 0, 1, thresholds=[6, 8])
    assert result["bin_labels"] == ["All distances"]
    assert result["bin_counts"] == [5]
    assert result["summary"] is None
    assert result["histogram"]["counts"] == []
    cluster = read_cluster_frames(path)
    valid = np.all(cluster.labels >= 0, axis=1)
    labels = cluster.labels[valid]
    original = pairwise_joints_on_edges(labels, [2, 2, 1], [(0, 1)])[(0, 1)]
    single_bin = pairwise_joints_on_edges(labels[indices[valid] == 0], [2, 2, 1], [(0, 1)])[(0, 1)]
    np.testing.assert_array_equal(single_bin, original)
    assert js_divergence(original.ravel(), single_bin.ravel()) == pytest.approx(0)


def test_preset_compatibility_and_corrupt_cache_are_explicit(dataset):
    store, path = dataset
    build_distance_cache("p", "s", "c", store=store, max_distance_cutoff=6)
    cache = load_distance_cache(path)
    validate_distance_preset(cache, [5, 6])
    with pytest.raises(ValueError, match="candidate-edge cutoff"):
        validate_distance_preset(cache, [5, 8])
    with pytest.raises(ValueError, match="ca_ca"):
        validate_distance_preset(cache, [5], distance_type="com_com")
    path.with_name("distance_cache.npz").write_bytes(b"broken")
    with pytest.raises(ValueError):
        load_distance_cache(path)


def test_failed_refresh_preserves_previous_cache(dataset):
    store, path = dataset
    result = build_distance_cache("p", "s", "c", store=store, max_distance_cutoff=6)
    cache_path = Path(result["path"])
    original = cache_path.read_bytes()
    with pytest.raises(FileNotFoundError):
        build_distance_cache("p", "s", "c", store=store, trajectory_overrides={"b": path.parent / "missing.xtc"})
    assert cache_path.read_bytes() == original


def test_empty_candidate_set(dataset):
    store, path = dataset
    result = build_distance_cache("p", "s", "c", store=store, max_distance_cutoff=0.1)
    assert result["metadata"]["edge_count"] == 0
    cache = load_distance_cache(path)
    assert cache.edge_index.shape == (2, 0)
    assert cache.distances.shape == (5, 0)


def test_periodic_distances_use_minimum_image(dataset):
    store, path = dataset
    state_path = path.parents[2] / "inputs/periodic"
    _trajectory(state_path, [49], periodic=True)
    # Exercise the same builder with a one-frame cluster and stored topology only.
    system = store.get_system("p", "s")
    system.states["a"].pdb_file = "inputs/periodic/structure.pdb"
    system.states["a"].trajectory_file = "inputs/periodic/trajectory.xtc"
    system.states["a"].slice_spec = None
    system.states["a"].stride = 1
    store.save_system(system)
    np.savez_compressed(path, residue_keys=np.asarray(["res_10", "res_20", "res_30"]),
        merged__labels=np.zeros((1, 3), dtype=np.int32),
        merged__frame_state_ids=np.asarray(["a"]), merged__frame_indices=np.asarray([0]))
    build_distance_cache("p", "s", "c", store=store, max_distance_cutoff=2)
    cache = load_distance_cache(path)
    np.testing.assert_array_equal(cache.edge_index, [[0], [1]])
    np.testing.assert_allclose(cache.distances, [[1]], atol=1e-5)
    system.states["a"].trajectory_file = None
    store.save_system(system)
    build_distance_cache("p", "s", "c", store=store, max_distance_cutoff=50, periodic=False)
    assert load_distance_cache(path).distances[0, 0] == pytest.approx(49)


def test_cli_build_status_and_filtered_histogram(dataset, capsys):
    store, _ = dataset
    base = ["--root", str(store.base_dir.parent), "--project-id", "p", "--system-id", "s", "--cluster-id", "c"]
    assert distance_cache_main([*base, "status"]) == 0
    assert json.loads(capsys.readouterr().out)["fallback"] == "single_unconditional_bin"
    assert distance_cache_main([*base, "build", "--max-distance-cutoff", "6", "--no-pbc"]) == 0
    assert json.loads(capsys.readouterr().out)["metadata"]["frame_count"] == 5
    assert distance_cache_main([*base, "pair", "--pair", "res_10", "res_20", "--state-ids", "b", "--thresholds", "5,6"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["bin_counts"] == [0, 1, 1]
    assert payload["frame_count"] == 2
