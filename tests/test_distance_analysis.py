import json
from pathlib import Path

import numpy as np
import pytest

from phase.services.project_store import ProjectStore
from phase.workflows.distance_analysis import (
    distance_pair_table, pair_detail, run_distance_js_analysis,
)


@pytest.fixture
def pair_dataset(tmp_path):
    store = ProjectStore(tmp_path / "projects")
    store.create_project("Project", project_id="p")
    system = store.create_system("p", name="System", system_id="s")
    cluster = tmp_path / "projects/p/systems/s/clusters/c"
    cluster.mkdir(parents=True)
    (cluster / "cluster_metadata.json").write_text(json.dumps({"cluster_id": "c", "path": "clusters/c/cluster.npz"}))
    # Global marginals are identical. In A short distance means c0, whereas
    # in B it means c1. Distance conditioning reveals a maximally different
    # joint state in both bins despite zero unconditional pair JS.
    labels = np.asarray([[0, 0], [0, 0], [1, 0], [1, 0], [1, 0], [1, 0], [0, 0], [0, 0]])
    states = np.asarray(["a"] * 4 + ["b"] * 4)
    metas = np.asarray(["near", "near", "far", "far"] * 2)
    frames = np.asarray([0, 1, 2, 3] * 2)
    path = cluster / "cluster.npz"
    np.savez_compressed(path, residue_keys=np.asarray(["res_1", "res_2"]), merged__labels_assigned=labels,
        merged__frame_state_ids=states, merged__frame_metastable_ids=metas,
        merged__frame_indices=frames, cluster_counts=np.asarray([2, 1]), contact_edge_index=np.asarray([[0], [1]]))
    np.savez_compressed(cluster / "distance_cache.npz", residue_keys=np.asarray(["res_1", "res_2"]),
        edge_index=np.asarray([[0], [1]]), frame_state_ids=states, frame_metastable_ids=metas, frame_indices=frames,
        distances=np.asarray([[5], [5], [9], [9]] * 2, dtype=np.float32),
        metadata_json=np.asarray(json.dumps({"schema_version": 1, "distance_type": "ca_ca", "distance_unit": "angstrom", "max_distance_cutoff": 12})))
    return store, path


def _run(store, **kwargs):
    return run_distance_js_analysis("p", "s", "c", store=store,
        condition_a={"state_ids": ["a"]}, condition_b={"state_ids": ["b"]}, **kwargs)


def test_conditioning_reveals_separation_hidden_by_global_js(pair_dataset):
    store, path = pair_dataset
    original = path.read_bytes()
    output = _run(store, thresholds=[6, 8])
    assert path.read_bytes() == original
    with np.load(output["analysis_npz"]) as result:
        assert result["edge_bin_js"].shape == (3, 1)
        np.testing.assert_allclose(result["edge_global_js"], [0], atol=1e-10)
        assert result["edge_bin_js"][0, 0] == pytest.approx(np.log(2), abs=1e-9)
        assert result["edge_bin_js"][2, 0] == pytest.approx(np.log(2), abs=1e-9)
        assert np.isnan(result["edge_bin_js"][1, 0])
        np.testing.assert_array_equal(result["edge_bin_counts_a"][:, 0], [2, 0, 2])
        np.testing.assert_allclose(result["edge_bin_occupancy_b"][:, 0], [.5, 0, .5])
    assert "/analyses/distance_js/" in output["analysis_npz"]
    meta = json.loads(Path(output["analysis_npz"]).with_name("analysis_metadata.json").read_text())
    assert meta["analysis_type"] == "distance_js"
    assert meta["js_units"] == "nats"
    assert meta["summary"]["overlap_frames"] == 0


def test_pair_heatmaps_filtering_and_boundary_counts(pair_dataset):
    _, path = pair_dataset
    data = pair_detail(path, 0, 1, condition={"state_ids": ["a"]}, thresholds=[5, 8])
    assert data["bin_counts"] == [0, 2, 2]  # 5 belongs to [5,8).
    np.testing.assert_allclose(data["panels"][0]["probabilities"], [[.5], [.5]])
    np.testing.assert_allclose(data["panels"][2]["probabilities"], [[1], [0]])
    assert data["panels"][1]["count"] == 0
    far = pair_detail(path, 0, 1, condition={"state_ids": ["a"], "metastable_ids": ["far"]})
    assert far["frame_count"] == 2
    assert sum(far["histogram"]["counts"]) == 2
    table = distance_pair_table(path, state_ids=["a"], contact_cutoff=6)
    assert table["pairs"][0]["contact_occupancy"] == .5
    assert table["pairs"][0]["median"] == 7
    assert table["pairs"][0]["range"] == 4
    assert distance_pair_table(path, search="missing")["total"] == 0
    assert distance_pair_table(path, contact_cutoff=6, min_contact_fraction=.5)["total"] == 1
    assert distance_pair_table(path, contact_cutoff=6, min_contact_fraction=.6)["total"] == 0
    assert distance_pair_table(path, contact_cutoff=5, min_contact_fraction=0)["total"] == 0
    with pytest.raises(ValueError, match="between 0 and 1"):
        distance_pair_table(path, min_contact_fraction=2)


def test_null_control_distance_independent_joint_states(pair_dataset):
    store, path = pair_dataset
    with np.load(path) as data:
        payload = {key: data[key] for key in data.files}
    payload["merged__labels_assigned"][:, 0] = [0, 1, 0, 1] * 2
    np.savez_compressed(path, **payload)
    output = _run(store, thresholds=[6, 8])
    with np.load(output["analysis_npz"]) as data:
        np.testing.assert_allclose(data["edge_bin_js"][[0, 2]], [[0], [0]], atol=1e-10)


def test_single_bin_fallback_matches_global_and_invalid_labels_are_pair_local(pair_dataset):
    store, path = pair_dataset
    path.with_name("distance_cache.npz").unlink()
    with pytest.raises(ValueError, match="requires a distance cache"):
        distance_pair_table(path, min_contact_fraction=.1)
    with np.load(path) as source:
        payload = {key: source[key] for key in source.files}
    payload["residue_keys"] = np.asarray(["res_1", "res_2", "res_3"])
    payload["merged__labels_assigned"] = np.column_stack([payload["merged__labels_assigned"], np.full(8, -1)])
    payload["cluster_counts"] = np.asarray([2, 1, 1])
    np.savez_compressed(path, **payload)
    output = _run(store, thresholds=[6, 8])
    assert output["metadata"]["distance_available"] is False
    with np.load(output["analysis_npz"]) as data:
        np.testing.assert_array_equal(data["distance_bin_labels"], ["All distances"])
        np.testing.assert_allclose(data["edge_bin_js"][0], data["edge_global_js"])
        assert data["edge_valid_counts_a"][0] == 4
    data = pair_detail(path, 0, 1, thresholds=[6, 8])
    assert data["bin_counts"] == [8]
    assert not data["histogram"]["counts"]


def test_empty_or_small_bins_are_unavailable_and_bad_conditions_fail(pair_dataset):
    store, _ = pair_dataset
    output = _run(store, thresholds=[6], min_bin_frames=3)
    with np.load(output["analysis_npz"]) as data:
        assert np.isnan(data["edge_bin_js"]).all()
        assert np.isfinite(data["edge_global_js"]).all()
    with pytest.raises(ValueError, match="Unknown frame filter"):
        run_distance_js_analysis("p", "s", "c", store=store,
            condition_a={"state_ids": ["newly_assigned_not_in_merged"]}, condition_b={"state_ids": ["b"]})
    with pytest.raises(ValueError, match="at least one"):
        run_distance_js_analysis("p", "s", "c", store=store,
            condition_a={"state_ids": []}, condition_b={"state_ids": ["b"]})


def test_no_geometry_or_contact_graph_requires_cache_build(pair_dataset):
    store, path = pair_dataset
    path.with_name("distance_cache.npz").unlink()
    with np.load(path) as data:
        payload = {key: data[key] for key in data.files if key != "contact_edge_index"}
    np.savez_compressed(path, **payload)
    with pytest.raises(ValueError, match="Build the distance cache first"):
        _run(store)


def test_a_bin_present_in_only_one_condition_is_not_zero_js(pair_dataset):
    store, path = pair_dataset
    cache = path.with_name("distance_cache.npz")
    with np.load(cache) as data:
        payload = {key: data[key] for key in data.files}
    payload["distances"][4:] = 5
    np.savez_compressed(cache, **payload)
    result = _run(store, thresholds=[6])
    with np.load(result["analysis_npz"]) as data:
        assert data["edge_bin_counts_a"][1, 0] == 2
        assert data["edge_bin_counts_b"][1, 0] == 0
        assert np.isnan(data["edge_bin_js"][1, 0])
