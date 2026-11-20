import pytest

from phase.potts.model_pairs import model_pairs, run_model_pairs


def test_all_unique_pairs_and_orientation():
    assert model_pairs(['a', 'b', 'a', 'c']) == [('a', 'b'), ('a', 'c'), ('b', 'c')]
    assert model_pairs(model_a_id='b', model_b_id='a') == [('b', 'a')]


@pytest.mark.parametrize('kwargs', [{}, {'model_ids': ['a', 'a']},
                                  {'model_ids': ['a', 'b'], 'model_a_id': 'a'}])
def test_invalid_selection(kwargs):
    with pytest.raises(ValueError):
        model_pairs(**kwargs)


def test_failed_pair_does_not_abort_remaining_pairs():
    called, messages = [], []
    def run(a, b):
        called.append((a, b))
        if b == 'b':
            raise ValueError('bad reference')
    assert run_model_pairs(model_pairs(['a', 'b', 'c']), run, messages.append) == 1
    assert called == [('a', 'b'), ('a', 'c'), ('b', 'c')]
    assert '2/3 pairs succeeded' in messages[-1]


def test_delta_cli_batch_preserves_sample_selection_and_limits(monkeypatch):
    from phase.scripts import potts_delta_energy as cli
    calls = []
    def run(**kwargs):
        calls.append(kwargs)
        return {'metadata': {}}
    monkeypatch.setattr(cli, 'upsert_delta_energy_analysis', run)
    assert cli.main(['--project-id', 'p', '--system-id', 's', '--cluster-id', 'c',
                     '--model-ids', 'a,b,c', '--sample-ids', 'md1,md2', '--frame-limit', 'md2:12']) == 0
    assert [(c['model_a_ref'], c['model_b_ref']) for c in calls] == model_pairs(['a', 'b', 'c'])
    assert all(c['sample_ids'] == ['md1', 'md2'] and c['frame_limits'] == {'md2': 12} for c in calls)


def test_hotspot_cli_batch_uses_references_of_each_model(monkeypatch, tmp_path):
    from phase.scripts.potts_model_pair_hotspots import main
    from phase.potts import hotspots
    calls = []
    monkeypatch.setattr(hotspots, 'run_hotspot_analysis', lambda **kwargs: calls.append(kwargs) or {})
    assert main(['--root', str(tmp_path), '--project-id', 'p', '--system-id', 's', '--cluster-id', 'c',
                 '--model-ids', 'a,b,c', '--model-reference', 'a:md1', '--model-reference', 'b:md2',
                 '--model-reference', 'c:md3,md4', '--sample-ids', 'extra']) == 0
    assert [(c['model_a_id'], c['model_b_id']) for c in calls] == model_pairs(['a', 'b', 'c'])
    assert [(c['reference_a'], c['reference_b']) for c in calls] == [
        (['md1'], ['md2']), (['md1'], ['md3', 'md4']), (['md2'], ['md3', 'md4'])]
    assert all(c['sample_ids'] == ['extra'] for c in calls)
