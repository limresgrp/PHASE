import json

import numpy as np
import pytest

from phase.potts.masking import MISSING_LABEL
from phase.potts.potts_model import PottsModel, save_potts_model
from phase.potts.analysis_run import compute_sample_energies, _compute_delta_energy_components, _delta_energy_sample_worker, compute_md_vs_sample_metrics
from phase.potts.sample_io import load_sample_npz, save_sample_npz
from phase.workflows.clustering import evaluate_state_with_models, _predict_labels_with_model
from phase.workflows.reference_states import create_reference_state
from tests.test_reference_alignment import pdb, reference_system


def model():
    return PottsModel(h=[np.array([1., 3.]), np.array([5., 7.]), np.array([11., 13.])],
                      J={(0, 1): np.full((2, 2), 17.), (0, 2): np.full((2, 2), 19.)}, edges=[(0, 1), (0, 2)])


def test_missing_energies_zero_only_incident_terms():
    potts = model()
    labels = np.array([[1, MISSING_LABEL, 0], [MISSING_LABEL, 0, 1], [-2, -2, -2]])
    result = compute_sample_energies(potts, labels)
    np.testing.assert_allclose(result['energies'], [3 + 11 + 19, 5 + 13, 0])
    np.testing.assert_allclose(result['node_energies'], [[3, 0, 11], [0, 5, 13], [0, 0, 0]])
    np.testing.assert_allclose(result['edge_energies'], [[0, 19], [0, 0], [0, 0]])
    np.testing.assert_allclose(potts.energy_batch(labels), result['energies'])
    assert potts.energy(labels[0]) == result['energies'][0]
    assert result['observed_node_counts'].tolist() == [2, 2, 0]
    assert result['observed_edge_counts'].tolist() == [1, 0, 0]
    with pytest.raises(ValueError, match='Invalid cluster label'):
        compute_sample_energies(potts, np.array([[0, -1, 0]]))


def test_missing_delta_components_and_worker(tmp_path):
    a, b = model(), model()
    b.h = [field * 2 for field in b.h]
    labels = np.array([[1, -2, 0], [0, 1, -2]])
    node, edge, edges = _compute_delta_energy_components(labels, a, b)
    np.testing.assert_allclose(node, [[-3, 0, -11], [-1, -7, 0]])
    np.testing.assert_allclose(edge, 0)
    save_potts_model(a, tmp_path/'a.npz')
    save_potts_model(b, tmp_path/'b.npz')
    save_sample_npz(tmp_path/'sample.npz', labels=labels, frame_indices=np.array([8, 19]))
    result = _delta_energy_sample_worker(dict(model_a_path=tmp_path/'a.npz', model_b_path=tmp_path/'b.npz',
                                             sample_path=tmp_path/'sample.npz', sample_id='pose'))
    assert result['used_frame_count'] == 2
    assert result['source_frame_indices'].tolist() == [8, 19]
    assert result['missing_mask'].tolist() == [[False, True, False], [False, False, True]]


def test_sample_mask_roundtrip_and_validation(tmp_path):
    path = tmp_path/'sample.npz'
    save_sample_npz(path, labels=np.array([[0, -2]]), labels_halo=np.array([[-1, -2]]))
    loaded = load_sample_npz(path)
    assert loaded.missing_mask.tolist() == [[False, True]]
    assert loaded.invalid_mask is None
    with pytest.raises(ValueError, match='missing_mask'):
        save_sample_npz(path, labels=np.array([[0, -2]]), missing_mask=np.array([[False, False]]))


def test_js_does_not_invent_distribution_for_missing_residues():
    result = compute_md_vs_sample_metrics(np.array([[0, -2, 1], [1, -2, 0]]),
                                          np.array([[0, 1, 1], [1, 0, 0]]),
                                          K=[2, 2, 2], edges=[(0, 1), (0, 2)])
    assert np.isnan(result['node_js'][1])
    assert np.isnan(result['edge_js'][0])
    assert result['node_js_mean'] == pytest.approx(0)
    assert result['edge_js_mean'] == pytest.approx(0)
    assert result['combined_distance'] == pytest.approx(0)


def test_assignment_preserves_incomplete_pdb_frames(reference_system, tmp_path, monkeypatch):
    store, _ = reference_system
    state = create_reference_state(store, 'p', 's', name='Partial',
        pdb_path=pdb(tmp_path/'partial.pdb', 'ACDEGHIKLMNPQRSTVWY', start=100), build_descriptors=True)
    keys = store.get_system('p', 's').reference_residue_keys
    directory = store.ensure_cluster_directories('p', 's', 'c')['cluster_dir']
    np.savez(directory/'cluster.npz', metadata_json=json.dumps({'residue_keys': keys}))
    fake = dict(kind='frozen_gmm', means=[[0]], covariances=[[[1]]], weights=[1])
    monkeypatch.setattr('phase.workflows.clustering._load_residue_models_from_metadata', lambda **kwargs: [fake]*len(keys))
    def predict(model, samples):
        assert np.all(np.isfinite(samples))
        return np.zeros(len(samples), dtype=np.int32), np.zeros(len(samples), dtype=np.int32)
    monkeypatch.setattr('phase.workflows.clustering._predict_cluster_frozen_gmm', predict)
    entry = evaluate_state_with_models('p', 's', 'c', state.state_id, store=store, workers=1)
    sample = load_sample_npz(store.resolve_path('p', 's', entry['path']))
    assert sample.labels.shape == (1, len(keys))
    assert sample.labels[0, keys.index('res_12')] == MISSING_LABEL
    assert sample.missing_mask.sum() == 1
    assert entry['missing_residue_keys'] == ['res_12']


def test_frustration_and_delta_store_missing_masks(tmp_path, monkeypatch):
    from phase.potts.analysis_run import upsert_delta_energy_analysis, upsert_endpoint_frustration_analysis
    from phase.services.project_store import ProjectStore
    from tests.test_endpoint_frustration_analysis import _write_model, _write_sample
    monkeypatch.setenv('PHASE_DATA_ROOT', str(tmp_path))
    store = ProjectStore(tmp_path/'projects')
    store.create_project('Test', project_id='p')
    system = store.create_system('p', name='Test', system_id='s')
    dirs = store.ensure_cluster_directories('p', 's', 'c')
    np.savez(dirs['cluster_dir']/'cluster.npz', residue_keys=['res_1', 'res_2', 'res_3'], cluster_counts=[2, 2, 2])
    a = model()
    b = model()
    b.h = [field * -1 for field in b.h]
    for name, potts in [('a', a), ('b', b)]:
        _write_model(dirs['potts_models_dir'], dirs['system_dir'], name, name, h=potts.h, J=potts.J, edges=potts.edges)
    _write_sample(dirs['system_dir'], 'c', 'partial', {'type': 'md_eval', 'name': 'Partial'}, labels=np.array([[0, -2, 0], [1, -2, 1]]))
    common = dict(project_id='p', system_id='s', cluster_id='c', model_a_ref='a', model_b_ref='b', sample_ids=['partial'])
    delta = upsert_delta_energy_analysis(**common)
    assert delta['metadata']['sample_missing_residue_indices'] == {'partial': [1]}
    component_path = dirs['system_dir']/delta['metadata']['component_sample_paths']['partial']
    with np.load(component_path) as data:
        assert data['missing_mask'][:, 1].all()
        assert np.all(data['delta_node_energy'][:, 1] == 0)
        assert np.all(data['delta_edge_energy'][:, 0] == 0)
    frustration = upsert_endpoint_frustration_analysis(**common)
    with np.load(frustration['analysis_npz']) as data:
        assert np.isnan(data['q_residue_all'][0, 1])
        assert np.isnan(data['frustration_node_sym_mean'][0, 1])
        assert np.isfinite(data['frustration_node_sym_mean'][0, 0])


def test_delta_js_partial_references_pool_only_observed_frames(tmp_path, monkeypatch):
    from phase.potts.analysis_run import upsert_delta_js_analysis
    from phase.services.project_store import ProjectStore
    from tests.test_endpoint_frustration_analysis import _write_sample
    monkeypatch.setenv('PHASE_DATA_ROOT', str(tmp_path))
    store = ProjectStore(tmp_path/'projects')
    store.create_project('Test', project_id='p')
    store.create_system('p', name='Test', system_id='s')
    dirs = store.ensure_cluster_directories('p', 's', 'c')
    np.savez(dirs['cluster_dir']/'cluster.npz', residue_keys=['res_1', 'res_2'], cluster_counts=[2, 2])
    for name, labels in [('partial', [[0, -2], [0, -2]]), ('a', [[0, 0]]), ('b', [[1, 1]]), ('other', [[0, -2]])]:
        _write_sample(dirs['system_dir'], 'c', name, {'type': 'md_eval', 'name': name}, labels=np.asarray(labels))
    result = upsert_delta_js_analysis(project_id='p', system_id='s', cluster_id='c', sample_ids=['other'],
        reference_sample_ids_a=['partial', 'a'], reference_sample_ids_b=['b'], edge_mode='all_vs_all')
    with np.load(result['analysis_npz']) as data:
        np.testing.assert_allclose(data['p_node_ref_a'][1, :2], [1, 0])
        assert np.isnan(data['js_node_a'][0, 1])
        assert np.isnan(data['js_edge_a'][0, 0])
