import json
import os
from types import SimpleNamespace

import pytest

os.environ.setdefault("PHASE_DATA_ROOT", "/tmp/phase-distance-api-tests")
fastapi = pytest.importorskip("fastapi")
pytest.importorskip("httpx")
from fastapi.testclient import TestClient

from backend.api.v1.routes import distances
from backend import distance_tasks
from tests.test_distance_analysis import pair_dataset
from tests.test_distance_cache import dataset


@pytest.fixture
def client(pair_dataset, monkeypatch):
    store, path = pair_dataset
    monkeypatch.setattr(distances, "project_store", store)
    app = fastapi.FastAPI()
    app.include_router(distances.router)
    enqueued = []
    class Queue:
        def enqueue(self, function, **kwargs):
            enqueued.append((function, kwargs))
            return SimpleNamespace(id=kwargs["job_id"])
    app.dependency_overrides[distances.get_queue] = lambda: Queue()
    with TestClient(app) as instance:
        yield instance, enqueued, store, path


BASE = "/projects/p/systems/s/metastable/clusters/c/distance"


def test_metadata_pair_query_and_persistent_presets(client):
    http, _, _, path = client
    assert http.get(f"{BASE}/cache").json()["distance_available"] is True
    pairs = http.get(f"{BASE}/pairs", params={"state_ids": "a", "limit": 1}).json()
    assert pairs["total"] == 1
    assert pairs["pairs"][0]["contact_occupancy"] == .5
    assert http.get(f"{BASE}/pairs", params={"contact_cutoff": 6, "min_contact_fraction": .6}).json()["total"] == 0
    detail = http.post(f"{BASE}/pairs/0/1", json={"condition": {"state_ids": ["a"]}, "thresholds": [6, 8]}).json()
    assert detail["bin_counts"] == [2, 0, 2]
    saved = http.post(f"{BASE}/presets", json={"name": "Near/far", "thresholds": [6, 8]}).json()
    assert (path.parent / "ui_setups" / f"{saved['setup_id']}.json").is_file()
    assert http.get(f"{BASE}/presets").json()["presets"][0]["compatible"] is True
    assert http.post(f"{BASE}/presets", json={"name": "Too far", "thresholds": [20]}).status_code == 400
    assert http.post(f"{BASE}/pairs/0/1", json={"condition": {"state_ids": ["missing"]}}).status_code == 400
    assert http.delete(f"{BASE}/presets/{saved['setup_id']}").status_code == 200
    assert not http.get(f"{BASE}/presets").json()["presets"]


def test_structure_mapping_uses_state_offsets_and_actual_topology(dataset, monkeypatch):
    store, _ = dataset
    monkeypatch.setattr(distances, "project_store", store)
    app = fastapi.FastAPI()
    app.include_router(distances.router)
    with TestClient(app) as http:
        response = http.get(f"{BASE}/structure-map/b")
        assert response.status_code == 200
        residues = response.json()["residues"]
        assert [r["key"] for r in residues] == ["res_10", "res_20", "res_30"]
        assert [r["resid"] for r in residues] == [110, 120, 130]
        assert [r["resname"] for r in residues] == ["ALA", "GLY", "SER"]
        assert residues[1]["position"] == [5, 0, 0]
        assert http.get(f"{BASE}/structure-map/missing").status_code == 400


def test_jobs_are_enqueued_without_processing_trajectories_in_http(client):
    http, queued, _, _ = client
    response = http.post(f"{BASE}/cache", json={"max_distance_cutoff": 12})
    assert response.status_code == 200
    response = http.post(f"{BASE}/analyses", json={"condition_a": {"state_ids": ["a"]}, "condition_b": {"state_ids": ["b"]}, "thresholds": [6, 8]})
    assert response.status_code == 200
    assert len(queued) == 2
    assert all(item[0] == "backend.distance_tasks.run_distance_job" for item in queued)
    assert queued[0][1]["args"][2] == "distance_cache"
    assert queued[1][1]["args"][2] == "distance_js"


def test_missing_cache_and_contact_edges_rejected_before_enqueue(client):
    http, queued, _, path = client
    path.with_name("distance_cache.npz").unlink()
    import numpy as np
    with np.load(path) as data:
        payload = {key: data[key] for key in data.files if key != "contact_edge_index"}
    np.savez_compressed(path, **payload)
    assert http.get(f"{BASE}/cache").json()["pair_count"] == 0
    response = http.post(f"{BASE}/analyses", json={
        "condition_a": {"state_ids": ["a"]}, "condition_b": {"state_ids": ["b"]}})
    assert response.status_code == 400
    assert "Build / refresh cache" in response.json()["detail"]
    assert not queued


def test_rq_task_runs_on_one_worker_and_saves_native_job_result(pair_dataset, monkeypatch):
    store, _ = pair_dataset
    monkeypatch.setattr(distance_tasks, "ProjectStore", lambda: store)
    class Job:
        id = "distance-js-test"
        meta = {}
        def save_meta(self):
            pass
    job = Job()
    monkeypatch.setattr(distance_tasks, "get_current_job", lambda: job)
    result = distance_tasks.run_distance_job("job-test", {"project_id": "p", "system_id": "s", "cluster_id": "c"}, "distance_js",
        {"condition_a": {"state_ids": ["a"]}, "condition_b": {"state_ids": ["b"]}, "thresholds": [6, 8]})
    assert result["status"] == "finished"
    assert job.meta["progress"] == 100
    assert result["results"]["metadata"]["analysis_type"] == "distance_js"
    saved = store.ensure_results_directories("p", "s")["jobs_dir"] / "job-test.json"
    assert json.loads(saved.read_text())["status"] == "finished"
