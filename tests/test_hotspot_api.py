import json
import os
from types import SimpleNamespace

import pytest

os.environ.setdefault('PHASE_DATA_ROOT', '/tmp/phase-hotspot-api-tests')
fastapi = pytest.importorskip('fastapi')
pytest.importorskip('httpx')
from fastapi.testclient import TestClient
from backend.api.v1.routes import hotspots
from backend import hotspot_tasks
from backend.api.v1.analysis_cleanup import _analysis_is_orphan
from tests.test_model_pair_hotspots import store

BASE = '/projects/p/systems/s/metastable/clusters/c/hotspots'


@pytest.fixture
def client(store, monkeypatch):
    monkeypatch.setattr(hotspots, 'project_store', store)
    app = fastapi.FastAPI()
    app.include_router(hotspots.router)
    queued = []
    class Queue:
        def enqueue(self, function, **kwargs):
            queued.append((function, kwargs))
            return SimpleNamespace(id=kwargs['job_id'])
    app.dependency_overrides[hotspots.get_queue] = lambda: Queue()
    with TestClient(app) as http:
        yield http, queued, store


def test_context_and_submit_only_enqueue(client, monkeypatch):
    http, queued, store = client
    store.model_entries[0]['params'] = {'state_ids': ['a']}
    response = http.get(BASE + '/context')
    assert response.status_code == 200
    assert response.json()['models'][0]['suggested_references'] == ['sa']
    response = http.post(BASE + '/analyses', json={'model_a_id': 'a', 'model_b_id': 'b', 'reference_a': ['sa'], 'reference_b': ['sb']})
    assert response.status_code == 200
    assert queued[0][0] == 'backend.hotspot_tasks.run_hotspot_job'
    assert queued[0][1]['args'][2]['block_length'] == 50
    assert not list(store.root.glob('analyses/*/*/analysis_metadata.json'))
    assert http.post(BASE + '/analyses', json={'model_a_id': 'a', 'model_b_id': 'a'}).status_code == 400
    assert http.post(BASE + '/analyses', json={'model_a_id': 'a', 'model_b_id': 'b', 'bootstrap': 5}).status_code == 400
    assert http.post(BASE + '/analyses', json={'model_a_id': 'a', 'model_b_id': 'b', 'block_length': 0}).status_code == 422


def test_one_worker_job_persists_success_failure_and_progress(store, monkeypatch):
    jobs = store.root / 'jobs'
    jobs.mkdir()
    store.ensure_results_directories = lambda *args: {'jobs_dir': jobs}
    monkeypatch.setattr(hotspot_tasks, 'ProjectStore', lambda: store)
    job = SimpleNamespace(id='hotspots-test', meta={}, save_meta=lambda: None)
    monkeypatch.setattr(hotspot_tasks, 'get_current_job', lambda: job)
    dataset = {'project_id': 'p', 'system_id': 's', 'cluster_id': 'c'}
    params = {'model_a_id': 'a', 'model_b_id': 'b', 'reference_a': ['sa'], 'reference_b': ['sb'], 'block_length': 2}
    result = hotspot_tasks.run_hotspot_job('success', dataset, params)
    assert result['status'] == 'finished'
    assert job.meta['progress'] == 100
    assert json.loads((jobs / 'success.json').read_text())['results']['analysis_id']
    with pytest.raises(ValueError):
        hotspot_tasks.run_hotspot_job('failure', dataset, {**params, 'model_b_id': 'a'})
    saved = json.loads((jobs / 'failure.json').read_text())
    assert saved['status'] == 'failed' and saved['error']


def test_cleanup_tracks_models_references_and_preference_samples():
    meta = {'analysis_type': 'model_pair_hotspots', 'model_a_id': 'a', 'model_b_id': 'b',
            'reference_sample_ids_a': ['sa'], 'reference_sample_ids_b': ['sb'], 'sample_ids': ['sa', 'sb', 'other']}
    assert not _analysis_is_orphan(meta, {'sa', 'sb', 'other'}, {'a', 'b'})
    assert _analysis_is_orphan(meta, {'sa', 'sb'}, {'a', 'b'})
    assert _analysis_is_orphan(meta, {'sa', 'sb', 'other'}, {'a'})


def test_structure_mapping_respects_exact_identity_and_missing_coordinates(client):
    http, _, store = client
    pdb = store.root / 'state.pdb'
    pdb.write_text('ATOM      1  CA  SER A 102       1.000   2.000   3.000  1.00  0.00           C\nEND\n')
    state = SimpleNamespace(pdb_file='state.pdb', residue_mapping={'res_2_SER': 'resindex 0'})
    store.get_system = lambda *args: SimpleNamespace(states={'offset': state}, reference_pdb_file=None)
    response = http.get(BASE + '/structure-map/offset')
    assert response.status_code == 200
    residues = response.json()['residues']
    assert residues[0]['resid'] == 102 and residues[0]['chain'] == 'A'
    assert residues[1] is None
    assert http.get(BASE + '/structure-map/not-a-state').status_code == 400
