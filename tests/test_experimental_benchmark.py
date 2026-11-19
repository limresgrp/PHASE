import json
from types import SimpleNamespace

import numpy as np
import pytest

from phase.analysis.experimental_benchmark import (percentiles, signed_auc, bh_adjust, pcoa, compute_benchmark,
                                                   run_experimental_benchmark, clean_json, correlation)
from phase.potts.potts_model import PottsModel, save_potts_model
from phase.services.analysis_annotations import save_scheme, list_schemes, get_scheme, delete_scheme, validate_scheme


def scientific_fixture():
    models = {
        'a': PottsModel([np.array([1., -1.]), np.array([2., -2.]), np.array([.4, -.4])], {(0, 2): np.array([[2., -2.], [-2., 2.]])}, [(0, 2)]),
        'b': PottsModel([np.array([0., 0.]), np.array([1., -1.]), np.array([.1, -.1])], {}, []),
        'c': PottsModel([np.array([-2., 2.]), np.array([0., 0.]), np.array([.3, -.3])], {(0, 1): np.array([[1., -1.], [-1., 1.]])}, [(0, 1)]),
    }
    labels = np.array([[0, 0, -2], [1, 0, 1], [1, 1, 0], [0, 1, 1]])
    states = ['s0', 's1', 's2', 's3']
    scheme = {'name': 'Supplied labels', 'kind': 'ordinal', 'categories': ['third', 'first', 'second'],
              'labels': {'s0': 'third', 's1': 'first', 's2': 'second', 's3': 'third'}}
    return models, labels, states, scheme


def compute(**options):
    models, labels, states, scheme = scientific_fixture()
    return compute_benchmark(models, [('a', 'b'), ('b', 'c'), ('c', 'a')], labels, states, scheme, **options)


def test_midrank_ties_missing_and_singleton():
    np.testing.assert_allclose(percentiles([1, 1, 3, np.nan]), [.25, .25, 1, np.nan], equal_nan=True)
    np.testing.assert_allclose(percentiles([3, 3]), [.5, .5])
    assert percentiles([3])[0] == .5
    assert np.isnan(percentiles([np.nan])[0])


def test_signed_auc_excludes_other_groups_and_preserves_orientation():
    labels = ['positive', 'positive', 'negative', 'negative', 'other']
    result = signed_auc([0, 1, 2, 3, -999], labels, 'positive', 'negative', True)
    assert result['auc'] == 0 and result['best_oriented_auc'] == 1 and result['n'] == 4
    assert signed_auc([1, 1], ['p', 'n'], 'p', 'n')['auc'] == .5
    assert np.isnan(signed_auc([1], ['p'], 'p', 'n')['auc'])
    with pytest.raises(ValueError):
        signed_auc([1], ['p'], 'p', 'p')


def test_bh_and_spearman_constants():
    np.testing.assert_allclose(bh_adjust([.01, .04, .03, np.nan]), [.03, .04, .04, np.nan], equal_nan=True)
    assert np.isnan(correlation([0, 0, 0], [1, 2, 3])['coefficient'])
    assert correlation([0, 1, 2], [2, 1, 0])['coefficient'] == -1


def test_common_support_missing_masks_and_exact_axis_algebra():
    result, arrays, _ = compute()
    assert arrays['node_observed'][0].tolist() == [True, True, False]
    assert np.isnan(arrays['term_contributions'][:, 0, 2]).all()
    assert arrays['common_support'].tolist() == [True, True, False, True, False]
    assert result['coverage_qc'][0]['common_nodes'] == 2
    assert result['coverage_qc'][0]['common_edges'] == 0
    assert result['axis_algebra'][0]['passed'] is True
    assert result['axis_algebra'][0]['max_residual'] < 1e-12
    assert arrays['common_scores'].shape == (4, 3)
    # Every stored local feature uses one fixed common term, never a missing term.
    assert arrays['common_support'][arrays['feature_term']].all()
    assert np.isfinite(arrays['fingerprints_raw']).all()


def test_ordinal_order_is_explicit_categorical_never_gets_spearman():
    models, labels, states, scheme = scientific_fixture()
    first, _, _ = compute_benchmark(models, [('a', 'b')], labels, states, scheme)
    reverse = {**scheme, 'categories': list(reversed(scheme['categories']))}
    second, _, _ = compute_benchmark(models, [('a', 'b')], labels, states, reverse)
    assert first['statistics'][0]['spearman']['coefficient'] == -second['statistics'][0]['spearman']['coefficient']
    categorical = {**scheme, 'kind': 'categorical'}
    third, _, _ = compute_benchmark(models, [('a', 'b')], labels, states, categorical)
    assert all('spearman' not in row for row in third['statistics'])


def test_distances_standardization_pcoa_and_raw_fingerprint_option():
    result, arrays, _ = compute()
    fp = arrays['fingerprints_raw']
    expected = np.divide(fp - fp.mean(axis=0), fp.std(axis=0), out=np.zeros_like(fp), where=fp.std(axis=0) > 1e-12)
    np.testing.assert_allclose(arrays['fingerprints_standardized'], expected)
    for key in ['euclidean_distance', 'correlation_distance']:
        np.testing.assert_allclose(arrays[key], arrays[key].T, equal_nan=True)
        np.testing.assert_array_equal(np.diag(arrays[key]), 0)
    assert result['pcoa']['available']
    coords = result['pcoa']['coordinates']
    np.testing.assert_allclose(np.linalg.norm(coords[:, None] - coords[None, :], axis=2), arrays['euclidean_distance'], atol=1e-8)
    _, raw, _ = compute(fingerprint_mode='raw')
    np.testing.assert_allclose(raw['euclidean_distance'], arrays['euclidean_distance'])


def test_undefined_correlation_is_not_imputed_and_pcoa_reports_failure():
    embedding = pcoa([[0, np.nan], [np.nan, 0]])
    assert not embedding['available']
    assert pcoa(np.zeros((3, 3)))['degenerate']
    with pytest.raises(ValueError):
        pcoa([[0, -1], [-1, 0]])


def test_no_common_support_preserves_raw_results_without_fake_zero_distance():
    models = {'a': PottsModel([np.array([1., -1.]), np.array([2., -2.])], {}, []),
              'b': PottsModel([np.zeros(2), np.zeros(2)], {}, [])}
    result, arrays, _ = compute_benchmark(models, [('a', 'b')], np.array([[0, -2], [-2, 1]]), ['x', 'y'],
                                        {'kind': 'categorical', 'categories': ['c'], 'labels': {}})
    assert np.isfinite(arrays['raw_scores']).all()
    assert np.isnan(arrays['common_scores']).all()
    assert np.isnan(arrays['euclidean_distance'][0, 1])
    assert not result['pcoa']['available']


def test_discordance_thresholds_and_redundancy_exports():
    result, arrays, _ = compute(global_threshold=1, local_threshold=0)
    assert all(pair['threshold_match'] for pair in result['discordance'])
    assert len(result['discordance']) == 6
    assert arrays['rank_spearman_matrix'].shape == (3, 3)
    assert arrays['common_score_pearson_matrix'].shape == (3, 3)
    with pytest.raises(ValueError):
        compute(global_threshold=2)


class Store:
    def __init__(self, root):
        self.root = root
        models, labels, state_ids, _ = scientific_fixture()
        self.states = {sid: SimpleNamespace(name=sid.upper(), n_frames=1, descriptor_file='descriptors.npz', pdb_file=f'{sid}.pdb') for sid in state_ids}
        self.model_entries, self.sample_entries = [], []
        for identifier, model in models.items():
            save_potts_model(model, str(root / f'{identifier}.npz'))
            self.model_entries.append({'model_id': identifier, 'path': f'{identifier}.npz', 'name': identifier.upper()})
        for sid, frame in zip(state_ids, labels):
            (root / f'{sid}.pdb').write_text('static source')
            np.savez(root / f'{sid}.npz', labels=frame[None, :], frame_indices=[0], frame_state_ids=[sid])
            self.sample_entries.append({'sample_id': sid, 'state_id': sid, 'name': sid.upper(), 'type': 'md_eval', 'path': f'{sid}.npz'})
        np.savez(root / 'cluster.npz', residue_keys=['res_2_SER', 'res_3_ALA', 'res_4_TYR'])
    def get_system(self, *args):
        return SimpleNamespace(states=self.states)
    def _system_dir(self, *args):
        return self.root
    def resolve_path(self, p, s, path):
        return self.root / path
    def get_cluster_entry(self, *args):
        return {'path': 'cluster.npz'}
    def list_potts_models(self, *args):
        return self.model_entries
    def list_samples(self, *args):
        return self.sample_entries
    def ensure_results_directories(self, *args):
        path = self.root / 'jobs'; path.mkdir(exist_ok=True)
        return {'jobs_dir': path}


@pytest.fixture
def store(tmp_path):
    return Store(tmp_path)


def run(store, **kwargs):
    _, _, state_ids, scheme = scientific_fixture()
    return run_experimental_benchmark('p', 's', 'c', store=store, state_ids=state_ids,
                                     contrasts=[{'model_a_id': 'a', 'model_b_id': 'b'}, {'model_a_id': 'b', 'model_b_id': 'c'}, {'model_a_id': 'a', 'model_b_id': 'c'}],
                                     label_scheme=scheme, **kwargs)


def test_full_workflow_immutable_atomic_artifacts_and_no_state_edits(store):
    before = [state.__dict__.copy() for state in store.states.values()]
    first = run(store)
    second = run(store)
    assert first['analysis_id'] != second['analysis_id']
    directory = store.root / 'analyses' / 'experimental_structure_benchmark' / first['analysis_id']
    assert (directory / 'scores.csv').is_file()
    assert (directory / 'hotspots.csv').is_file()
    assert (directory / 'euclidean_distance.csv').is_file()
    summary = json.loads((directory / 'summary.json').read_text())
    assert summary['axis_algebra'][0]['passed']
    assert summary['metadata']['immutable']
    assert all(len(m['sha256']) == 64 for m in summary['metadata']['models'])
    with np.load(directory / 'analysis.npz', allow_pickle=False) as arrays:
        assert arrays['term_contributions'].shape[0] == 3
        assert arrays['fingerprints_standardized'].shape[0] == 4
        assert 'hotspot_0_node' in arrays
    assert before == [state.__dict__ for state in store.states.values()]
    assert not list(store.root.glob('.experimental-benchmark-*'))


def test_label_schemes_are_persistent_and_snapshot_independent(store):
    _, _, _, scheme = scientific_fixture()
    saved = save_scheme(store, 'p', 's', scheme)
    update = save_scheme(store, 'p', 's', {**scheme, 'name': 'Changed'}, saved['scheme_id'])
    assert update['revision'] == 2 and list_schemes(store, 'p', 's')[0]['name'] == 'Changed'
    assert get_scheme(store, 'p', 's', saved['scheme_id'])['categories'] == scheme['categories']
    delete_scheme(store, 'p', 's', saved['scheme_id'])
    assert not list_schemes(store, 'p', 's')
    with pytest.raises(ValueError):
        validate_scheme({**scheme, 'categories': ['duplicate', 'duplicate']}, store.states)
    with pytest.raises(ValueError):
        validate_scheme({**scheme, 'labels': {'not-a-state': 'first'}}, store.states)


def test_trajectory_inputs_unassigned_and_raw_delta_models_are_rejected(store):
    np.savez(store.root / 's0.npz', labels=np.array([[0, 0, 0], [1, 1, 1]]))
    with pytest.raises(ValueError, match='exactly one'):
        run(store)
    np.savez(store.root / 's0.npz', labels=np.array([[-1, 0, 0]]))
    with pytest.raises(ValueError, match='Invalid cluster'):
        run(store)
    store.model_entries[0]['params'] = {'delta_kind': 'delta_patch'}
    with pytest.raises(ValueError, match='full/combined'):
        run(store)


def test_ensemble_hotspots_reuse_uncertainty_and_store_coverage(store):
    # Static endpoint references are valid mean observations but no valid CIs.
    for m in store.model_entries:
        m['params'] = {'state_ids': ['s0' if m['model_id'] == 'a' else 's1' if m['model_id'] == 'b' else 's2']}
    meta = run(store, hotspots_mode='ensemble', block_length=1, bootstrap=50)
    path = store.root / 'analyses' / meta['analysis_type'] / meta['analysis_id'] / 'analysis.npz'
    with np.load(path) as data:
        assert 'hotspot_0_counts_a' in data
        assert np.isnan(data['hotspot_0_ci_low']).all()


def test_clean_json_uses_null_for_undefined_scientific_values():
    assert clean_json({'x': np.array([np.nan, np.inf, 0.])}) == {'x': [None, None, 0.]}


def test_auto_assignment_is_single_core_and_persisted_before_use(store, monkeypatch):
    import sys
    sample_entry = store.sample_entries.pop(0)
    calls = []
    def evaluate(*args, **kwargs):
        calls.append(kwargs)
        return sample_entry
    def save(*args):
        store.sample_entries.append(args[-1])
    store.save_sample_entry = save
    monkeypatch.setitem(sys.modules, 'phase.workflows.clustering', SimpleNamespace(evaluate_state_with_models=evaluate))
    result = run(store)
    assert result['sample_ids'] == ['s0', 's1', 's2', 's3']
    assert calls[0]['workers'] == 1
    assert len(store.sample_entries) == 4


def test_atomic_failure_never_leaves_a_listed_half_analysis(store, monkeypatch):
    import phase.analysis.experimental_benchmark as module
    def fail(*args):
        raise OSError('disk full')
    monkeypatch.setattr(module, '_write_json', fail)
    with pytest.raises(OSError, match='disk full'):
        run(store)
    assert not list(store.root.glob('analyses/*/*'))
    assert not list(store.root.glob('.experimental-benchmark-*'))


def test_saved_scheme_snapshot_remains_unchanged_after_edit_and_delete(store):
    _, _, states, scheme = scientific_fixture()
    saved = save_scheme(store, 'p', 's', scheme)
    meta = run_experimental_benchmark('p', 's', 'c', store=store, state_ids=states,
            contrasts=[{'model_a_id': 'a', 'model_b_id': 'b'}], scheme_id=saved['scheme_id'])
    directory = store.root / 'analyses' / meta['analysis_type'] / meta['analysis_id']
    before = (directory / 'labels.json').read_bytes()
    save_scheme(store, 'p', 's', {**scheme, 'categories': list(reversed(scheme['categories']))}, saved['scheme_id'])
    delete_scheme(store, 'p', 's', saved['scheme_id'])
    assert (directory / 'labels.json').read_bytes() == before


def test_cli_config_and_annotation_management_use_explicit_root(store, monkeypatch, capsys):
    from phase.scripts.experimental_structure_benchmark import main
    import phase.services.project_store as module
    roots = []
    def factory(base_dir):
        roots.append(base_dir)
        return store
    monkeypatch.setattr(module, 'ProjectStore', factory)
    monkeypatch.setenv('PHASE_DATA_ROOT', 'original')
    config = store.root / 'benchmark.json'
    config.write_text(json.dumps({'state_ids': list(store.states), 'label_scheme': scientific_fixture()[3],
                                 'contrasts': [{'model_a_id': 'a', 'model_b_id': 'b'}]}))
    common = ['--root', str(store.root), '--project-id', 'p', '--system-id', 's']
    assert main(common + ['--cluster-id', 'c', '--config', str(config), '--progress']) == 0
    result = json.loads(capsys.readouterr().out)
    assert result['immutable'] and roots == [store.root / 'projects']
    assert main(common + ['--list-label-schemes']) == 0
    assert json.loads(capsys.readouterr().out) == []


def test_real_project_store_sidecars_and_result_location(tmp_path):
    from phase.services.project_store import ProjectStore, DescriptorState
    store = ProjectStore(tmp_path / 'projects')
    store.create_project('Project', project_id='p')
    system = store.create_system('p', name='System', system_id='s')
    directories = store.ensure_cluster_directories('p', 's', 'c')
    cluster = directories['cluster_dir']; root = directories['system_dir']
    models, frames, ids, scheme = scientific_fixture()
    entries = []
    for mid, model in models.items():
        path = cluster / 'potts_models' / mid / 'model.npz'; path.parent.mkdir(parents=True)
        save_potts_model(model, str(path))
        entries.append({'model_id': mid, 'name': mid, 'path': str(path.relative_to(root))})
    samples = []
    for sid, frame in zip(ids, frames):
        folder = root / 'states' / sid; folder.mkdir(parents=True)
        (folder / 'structure.pdb').write_text('source PDB')
        system.states[sid] = DescriptorState(state_id=sid, name=sid, pdb_file=f'states/{sid}/structure.pdb', n_frames=1)
        path = cluster / 'samples' / sid / 'sample.npz'; path.parent.mkdir(parents=True)
        np.savez(path, labels=frame[None, :], frame_state_ids=[sid])
        samples.append({'sample_id': sid, 'state_id': sid, 'name': sid, 'type': 'md_eval', 'path': str(path.relative_to(root))})
    np.savez(cluster / 'cluster.npz', residue_keys=['res_2_SER', 'res_3_ALA', 'res_4_TYR'])
    system.metastable_clusters = [{'cluster_id': 'c', 'path': 'clusters/c/cluster.npz', 'potts_models': entries, 'samples': samples}]
    store.save_system(system)
    result = run_experimental_benchmark('p', 's', 'c', store=store, state_ids=ids, label_scheme=scheme,
                                       contrasts=[{'model_a_id': 'a', 'model_b_id': 'b'}], assign_missing=False)
    path = cluster / 'analyses' / result['analysis_type'] / result['analysis_id']
    assert (path / 'summary.json').is_file()
    assert len(store.list_samples('p', 's', 'c')) == 4
    assert set(store.get_system('p', 's').states) == set(ids)
    assert result['parameters']['distance'] == 'euclidean'


def test_raw_masked_scores_and_common_scores_have_explicit_different_support():
    models = {'a': PottsModel([np.array([0., 2.]), np.array([100., 102.])], {}, []),
              'b': PottsModel([np.zeros(2), np.zeros(2)], {}, [])}
    result, arrays, _ = compute_benchmark(models, [('a', 'b')], np.array([[0, -2], [1, 1]]), ['one', 'two'],
                                        {'kind': 'categorical', 'categories': ['label'], 'labels': {}})
    np.testing.assert_allclose(arrays['raw_scores'][:, 0], [-1, 2])
    np.testing.assert_allclose(arrays['common_scores'][:, 0], [-1, 1])
    assert np.isnan(arrays['term_contributions'][0, 0, 1])
    assert result['coverage_qc'][0]['node_counts'].tolist() == [1, 2]


def test_global_ties_can_coexist_with_distinct_local_fingerprints():
    models = {'a': PottsModel([np.array([1., -1.]), np.array([1., -1.])], {}, []),
              'b': PottsModel([np.zeros(2), np.zeros(2)], {}, [])}
    result, arrays, _ = compute_benchmark(models, [('a', 'b')], np.array([[0, 1], [1, 0]]), ['one', 'two'],
            {'kind': 'categorical', 'categories': ['label'], 'labels': {}}, global_threshold=.1, local_threshold=1)
    np.testing.assert_allclose(arrays['common_scores'], 0)
    assert result['discordance'][0]['global_rank_difference'] == 0
    assert result['discordance'][0]['local_distance'] > 1
    assert result['discordance'][0]['threshold_match']


def test_shared_delta_decomposition_precision_and_masking_are_opt_in():
    from phase.potts.analysis_run import _compute_delta_energy_components
    model = PottsModel([np.array([1.0000000001, -1.]), np.array([2., -2.])], {(0, 1): np.array([[1.0000000001, 2.], [3., 4.]])}, [(0, 1)])
    zero = PottsModel([np.zeros(2), np.zeros(2)], {}, [])
    labels = np.array([[0, 0], [0, -2]])
    old_node, old_edge, _ = _compute_delta_energy_components(labels, model, zero)
    node, edge, _ = _compute_delta_energy_components(labels, model, zero, dtype=np.float64, mask_missing=True)
    assert old_node.dtype == np.float32 and old_edge.dtype == np.float32
    assert old_node[1, 1] == old_edge[1, 0] == 0
    assert node.dtype == np.float64 and edge.dtype == np.float64
    assert node[0, 0] == edge[0, 0] == 1.0000000001
    assert np.isnan(node[1, 1]) and np.isnan(edge[1, 0])
