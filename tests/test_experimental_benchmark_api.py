import io
import json
import os
from types import SimpleNamespace
import zipfile

import pytest

os.environ.setdefault('PHASE_DATA_ROOT', '/tmp/phase-benchmark-api-tests')
fastapi = pytest.importorskip('fastapi')
pytest.importorskip('httpx')
from fastapi.testclient import TestClient
from backend.api.v1.routes import experimental_benchmarks as routes
from backend import experimental_benchmark_tasks as tasks
from backend.api.v1.analysis_cleanup import _analysis_is_orphan
from tests.test_experimental_benchmark import store, run, scientific_fixture

BASE = '/projects/p/systems/s/experimental-benchmarks'


@pytest.fixture
def client(store, monkeypatch):
    monkeypatch.setattr(routes, 'project_store', store)
    app = fastapi.FastAPI(); app.include_router(routes.router)
    queued = []
    class Queue:
        def enqueue(self, function, **kwargs):
            queued.append((function, kwargs))
            return SimpleNamespace(id=kwargs['job_id'])
    app.dependency_overrides[routes.get_queue] = lambda: Queue()
    with TestClient(app) as http:
        yield http, queued, store


def request():
    _, _, states, scheme = scientific_fixture()
    return {'cluster_id': 'c', 'state_ids': states, 'label_scheme': scheme,
            'contrasts': [{'model_a_id': 'a', 'model_b_id': 'b'}]}


def test_system_scheme_crud_does_not_edit_states_and_explicit_order_roundtrips(client):
    http, _, store = client
    scheme = scientific_fixture()[3]
    created = http.post(BASE + '/label-schemes', json=scheme)
    assert created.status_code == 200, created.text
    saved = created.json()
    assert saved['categories'] == ['third', 'first', 'second']
    assert http.get(BASE + '/label-schemes').json()['schemes'][0]['scheme_id'] == saved['scheme_id']
    updated = http.put(BASE + '/label-schemes/' + saved['scheme_id'], json={**scheme, 'name': 'revision'}).json()
    assert updated['revision'] == 2
    assert store.states['s0'].name == 'S0'
    assert http.delete(BASE + '/label-schemes/' + saved['scheme_id']).status_code == 200
    assert not http.get(BASE + '/label-schemes').json()['schemes']
    assert http.get(BASE + '/context?cluster_id=c').json()['states'][0]['assigned_samples']


def test_submit_only_enqueues_and_freezes_label_snapshot(client):
    http, queued, store = client
    scheme = http.post(BASE + '/label-schemes', json=scientific_fixture()[3]).json()
    body = request(); del body['label_scheme']; body['scheme_id'] = scheme['scheme_id']
    response = http.post(BASE + '/analyses', json=body)
    assert response.status_code == 200, response.text
    assert queued[0][0] == 'backend.experimental_benchmark_tasks.run_benchmark_job'
    params = queued[0][1]['args'][2]
    assert params['scheme_id'] is None and params['label_scheme']['revision'] == 1
    assert not list(store.root.glob('analyses/*/*'))
    http.put(BASE + '/label-schemes/' + scheme['scheme_id'], json={**scientific_fixture()[3], 'categories': ['second', 'first', 'third']})
    assert params['label_scheme']['categories'] == ['third', 'first', 'second']
    assert http.post(BASE + '/analyses', json={**request(), 'positive_label': 'first'}).status_code == 400
    assert http.post(BASE + '/analyses', json={**request(), 'bootstrap': 1}).status_code == 400


def test_single_worker_job_persists_progress_result_and_download_artifacts(client, monkeypatch):
    http, _, store = client
    monkeypatch.setattr(tasks, 'ProjectStore', lambda: store)
    job = SimpleNamespace(id='benchmark-job', meta={}, save_meta=lambda: None)
    monkeypatch.setattr(tasks, 'get_current_job', lambda: job)
    body = request(); body.pop('cluster_id')
    result = tasks.run_benchmark_job('job', {'project_id': 'p', 'system_id': 's', 'cluster_id': 'c'}, body)
    assert result['status'] == 'finished' and job.meta['progress'] == 100
    identifier = result['results']['analysis_id']
    response = http.get(f'{BASE}/analyses/{identifier}?cluster_id=c')
    assert response.status_code == 200, response.text
    assert 'fingerprints_raw' not in response.json()
    downloaded = http.get(f'{BASE}/analyses/{identifier}/download?cluster_id=c')
    assert downloaded.status_code == 200
    with zipfile.ZipFile(io.BytesIO(downloaded.content)) as archive:
        assert {'analysis.npz', 'scores.csv', 'labels.json', 'supports.json', 'provenance.json', 'hotspots.csv'} <= set(archive.namelist())
    assert json.loads((store.root / 'jobs/job.json').read_text())['status'] == 'finished'
    with pytest.raises(ValueError):
        tasks.run_benchmark_job('failure', {'project_id': 'p', 'system_id': 's', 'cluster_id': 'c'}, {**body, 'state_ids': ['s0']})
    assert json.loads((store.root / 'jobs/failure.json').read_text())['status'] == 'failed'


def test_cleanup_references_models_states_and_samples_but_not_deleted_label_scheme():
    metadata = {'analysis_type': 'experimental_structure_benchmark', 'model_ids': ['a', 'b'], 'sample_ids': ['sample'],
                'state_ids': ['state'], 'label_scheme_snapshot': {'scheme_id': 'removed-scheme'}}
    assert not _analysis_is_orphan(metadata, {'sample'}, {'a', 'b'}, {'state'})
    assert _analysis_is_orphan(metadata, set(), {'a', 'b'}, {'state'})
    assert _analysis_is_orphan(metadata, {'sample'}, {'a'}, {'state'})
    assert _analysis_is_orphan(metadata, {'sample'}, {'a', 'b'}, set())


def test_arbitrary_annotation_identifiers_do_not_become_scientific_state_dependencies():
    metadata = {'analysis_type': 'experimental_structure_benchmark', 'state_ids': ['kept'],
                'label_scheme_snapshot': {'labels': {'state_ids': 'arbitrary annotation', 'removed-unselected-state': 'label'}}}
    assert not _analysis_is_orphan(metadata, set(), set(), {'kept'})
