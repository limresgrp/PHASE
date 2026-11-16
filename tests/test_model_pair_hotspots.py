import json

import numpy as np
import pytest

from phase.potts.hotspots import (block_statistics, contrast_models, run_hotspot_analysis,
                                 separation_statistics, summarize_blocks, suggested_references)
from phase.potts.potts_model import PottsModel, save_potts_model


def models():
    a = PottsModel([np.array([-1., 1.]), np.array([0., 0.])],
                   {(0, 1): np.array([[2., -2.], [-2., 2.]])}, [(0, 1)])
    b = PottsModel([np.zeros(2), np.zeros(2)], {}, [])
    return a, b


def test_gauge_union_alphabet_normalized_and_swap():
    a, b = models()
    delta = contrast_models(a, b)
    assert delta.edges == [(0, 1)]
    assert np.sqrt(np.mean(delta.J[0, 1] ** 2)) == 2
    assert np.sqrt(np.mean(delta.h[0] ** 2)) == 1
    np.testing.assert_allclose(contrast_models(b, a).J[0, 1], -delta.J[0, 1])
    # Move a row-dependent potential from field to coupling: same Hamiltonian.
    v = np.array([3., -4.])
    shifted = PottsModel([a.h[0] - v, a.h[1] + 20], {(0, 1): a.J[0, 1] + v[:, None]}, a.edges)
    shifted_delta = contrast_models(shifted, b)
    np.testing.assert_allclose(shifted_delta.h, delta.h)
    np.testing.assert_allclose(shifted_delta.J[0, 1], delta.J[0, 1])


def test_invalid_model_rejected():
    a, _ = models()
    with pytest.raises(ValueError, match='alphabets'):
        contrast_models(a, PottsModel([np.zeros(3)], {}, []))
    a.J[0, 1][0, 0] = np.nan
    with pytest.raises(ValueError, match='coupling'):
        contrast_models(a, a)


def test_masked_values_are_not_zero_observations():
    a, b = models()
    sums, counts = block_statistics(contrast_models(a, b), np.array([[0, -2], [1, 1]]), block_length=1)
    np.testing.assert_array_equal(counts, [[1, 0, 0], [1, 1, 1]])
    summary = summarize_blocks([(sums, counts)], 0)
    np.testing.assert_allclose(summary['mean'], [0, 0, 2])
    with pytest.raises(ValueError, match='Invalid cluster'):
        block_statistics(a, np.array([[-1, 0]]))


def test_gaps_and_state_boundaries_break_blocks():
    a, _ = models()
    sums, _ = block_statistics(a, np.zeros((5, 2), dtype=int),
                               [0, 1, 10, 11, 12], ['a', 'a', 'a', 'b', 'b'], block_length=100)
    assert len(sums) == 3


def test_equal_trajectory_weight_not_frame_weight():
    a = PottsModel([np.array([-1., 1.])], {}, [])
    long = block_statistics(a, np.zeros((100, 1), dtype=int), block_length=10)
    short = block_statistics(a, np.ones((2, 1), dtype=int), block_length=1)
    assert summarize_blocks([long, short], 0)['mean'][0] == 0


def test_separation_sign_determinism_and_single_frame_no_ci():
    a = PottsModel([np.array([-1., 1.])], {}, [])
    low = block_statistics(a, np.zeros((20, 1), dtype=int), block_length=2)
    high = block_statistics(a, np.ones((20, 1), dtype=int), block_length=2)
    result = separation_statistics([low], [high], seed=4)
    assert result['separation'][0] == 2
    assert result['ci_low'][0] == 2
    assert separation_statistics([high], [low])['separation'][0] == -2
    for key in result:
        np.testing.assert_allclose(result[key], separation_statistics([low], [high], seed=4)[key])
    one = block_statistics(a, np.ones((1, 1), dtype=int))
    assert np.isnan(separation_statistics([low], [one])['ci_low'][0])
    with pytest.raises(ValueError):
        separation_statistics([low], [high], bootstrap=1)


class Store:
    def __init__(self, root):
        self.root = root
        self.model_entries = []
        self.sample_entries = []
    def _system_dir(self, *args):
        return self.root
    def get_cluster_entry(self, *args):
        return {'path': 'cluster.npz'}
    def resolve_path(self, p, s, path):
        return self.root / path
    def list_potts_models(self, *args):
        return self.model_entries
    def list_samples(self, *args):
        return self.sample_entries


@pytest.fixture
def store(tmp_path):
    store = Store(tmp_path)
    np.savez(tmp_path / 'cluster.npz', residue_keys=['res_2_SER', 'res_3_ALA'])
    for identifier, model in zip(['a', 'b'], models()):
        save_potts_model(model, str(tmp_path / f'{identifier}.npz'))
        store.model_entries.append({'model_id': identifier, 'name': identifier.upper(), 'path': f'{identifier}.npz'})
        np.savez(tmp_path / f's{identifier}.npz', labels=np.full((12, 2), identifier == 'b', dtype=int),
                 frame_indices=np.arange(12), frame_state_ids=np.full(12, identifier))
        store.sample_entries.append({'sample_id': f's{identifier}', 'name': identifier, 'type': 'md_eval',
                                     'state_id': identifier, 'path': f's{identifier}.npz'})
    return store


def test_workflow_compact_artifact_and_progress(store):
    progress = []
    result = run_hotspot_analysis('p', 's', 'c', store=store, model_a_id='a', model_b_id='b',
                                  reference_a=['sa'], reference_b=['sb'], block_length=2,
                                  progress_callback=lambda *args: progress.append(args))
    assert result['edge_basis'] == 'union'
    assert result['name'] == 'A / B'
    assert progress[-1][0] == 'Complete'
    path = store.root / result['paths']['summary_npz']
    assert json.loads((path.parent / 'analysis_metadata.json').read_text())['sample_ids'] == ['sa', 'sb']
    with np.load(path, allow_pickle=False) as data:
        assert data['node_separation'].shape == (2,)
        assert data['edge_separation'].shape == (1,)
        assert data['node_preference'].shape == (2, 2)
        np.testing.assert_allclose(data['node_preference'][0], -data['node_preference'][1])
        assert all(not key.startswith('frame_') for key in data)


def test_parameter_mode_needs_no_samples_and_raw_delta_is_rejected(store):
    store.sample_entries = []
    result = run_hotspot_analysis('p', 's', 'c', store=store, model_a_id='a', model_b_id='b', mode='parameter')
    with np.load(store.root / result['paths']['summary_npz']) as data:
        assert 'node_separation' not in data
    store.model_entries[0]['params'] = {'delta_kind': 'delta_patch'}
    with pytest.raises(ValueError, match='raw delta'):
        run_hotspot_analysis('p', 's', 'c', store=store, model_a_id='a', model_b_id='b', mode='parameter')


def test_overlapping_references_and_wrong_basis_rejected(store):
    kwargs = dict(store=store, model_a_id='a', model_b_id='b', reference_a=['sa'], reference_b=['sb'])
    np.savez(store.root / 'sb.npz', labels=np.ones((12, 2), dtype=int), frame_indices=np.arange(12), frame_state_ids=np.full(12, 'a'))
    with pytest.raises(ValueError, match='overlap'):
        run_hotspot_analysis('p', 's', 'c', **kwargs)
    store.model_entries[0]['residue_keys'] = ['wrong', 'order']
    with pytest.raises(ValueError, match='ordering'):
        run_hotspot_analysis('p', 's', 'c', **kwargs)


def test_reference_suggestions(store):
    assert suggested_references({'params': {'state_ids': ['a']}}, store.sample_entries) == ['sa']


def test_bootstrap_has_uncertainty_for_variable_blocks():
    a = PottsModel([np.array([-1., 1.])], {}, [])
    mixed = block_statistics(a, np.array([[0], [0], [1], [1]] * 5), block_length=2)
    reference = block_statistics(a, np.zeros((20, 1), dtype=int), block_length=2)
    result = separation_statistics([reference], [mixed], bootstrap=200, seed=8)
    assert result['ci_low'][0] < result['separation'][0] < result['ci_high'][0]


def test_reference_suggestions_and_validation_of_embedded_metadata(store):
    a, _ = models()
    save_potts_model(a, str(store.root / 'a.npz'), metadata={'residue_keys': ['wrong', 'order']})
    with pytest.raises(ValueError, match='ordering'):
        run_hotspot_analysis('p', 's', 'c', store=store, model_a_id='a', model_b_id='b', mode='parameter')


def test_cli_uses_explicit_root_even_if_store_module_was_preimported(store, monkeypatch, capsys):
    from phase.scripts.potts_model_pair_hotspots import main
    import phase.services.project_store as service
    roots = []
    monkeypatch.setenv("PHASE_DATA_ROOT", "original-root")
    def factory(base_dir):
        roots.append(base_dir)
        return store
    monkeypatch.setattr(service, 'ProjectStore', factory)
    assert main(['--root', str(store.root), '--project-id', 'p', '--system-id', 's', '--cluster-id', 'c',
                 '--model-a-id', 'a', '--model-b-id', 'b', '--mode', 'parameter', '--progress']) == 0
    assert roots == [store.root / 'projects']
    assert json.loads(capsys.readouterr().out)['mode'] == 'parameter'


def test_preference_only_sample_uses_one_block_even_with_random_source_frames(store, monkeypatch):
    import phase.potts.hotspots as module
    np.savez(store.root / 'extra.npz', labels=np.zeros((100, 2), dtype=int), frame_indices=np.arange(100)[::-1])
    store.sample_entries.append({'sample_id': 'extra', 'name': 'SA reads', 'type': 'potts_sampling', 'path': 'extra.npz'})
    original = module.block_statistics
    lengths = []
    def track(*args, **kwargs):
        result = original(*args, **kwargs)
        lengths.append(len(result[0]))
        return result
    monkeypatch.setattr(module, 'block_statistics', track)
    result = run_hotspot_analysis('p', 's', 'c', store=store, model_a_id='a', model_b_id='b',
                                  reference_a=['sa'], reference_b=['sb'], sample_ids=['extra'], block_length=2, bootstrap=0)
    assert lengths == [6, 6, 1]
    with np.load(store.root / result['paths']['summary_npz']) as data:
        np.testing.assert_array_equal(data['sample_node_counts'][2], [100, 100])


def test_failed_artifact_write_never_publishes_partial_analysis(store, monkeypatch):
    import phase.potts.hotspots as module
    def fail(*args):
        raise OSError('simulated disk failure')
    monkeypatch.setattr(module, '_write_json', fail)
    with pytest.raises(OSError, match='disk failure'):
        run_hotspot_analysis('p', 's', 'c', store=store, model_a_id='a', model_b_id='b', mode='parameter')
    assert not list(store.root.glob('analyses/*/*'))
    assert not list(store.root.glob('.hotspots-*'))
