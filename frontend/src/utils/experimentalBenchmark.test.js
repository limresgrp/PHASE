import { axisPairs, benchmarkTraces, discordancePairs } from './experimentalBenchmark';
const result = {
  states: [{ state_id: 's', name: 'Structure' }, { state_id: 't', name: 'Other' }],
  label_scheme: { categories: ['User category'] },
  scores: [
    { state_id: 's', axis: 0, label: 'User category', raw_score: -100, common_score: -2, raw_rank: 0, common_rank: .25, observed_nodes: 10, observed_edges: 3 },
    { state_id: 's', axis: 1, label: 'User category', raw_score: 200, common_score: 3, raw_rank: 1, common_rank: .75, observed_nodes: 10, observed_edges: 4 },
    { state_id: 't', axis: 0, label: null, raw_score: null, common_score: null, raw_rank: null, common_rank: null },
    { state_id: 't', axis: 1, label: null, raw_score: 10, common_score: 2, raw_rank: 0, common_rank: 0 },
  ],
  discordance: [{ state_a: 's', state_b: 't', global_rank_difference: .1, local_distance: 3 },
    { state_a: 's', state_b: 'u', global_rank_difference: null, local_distance: 5 }],
};
test('generates all unordered axis pairs', () => {
  expect(axisPairs(3)).toEqual([[0, 1], [0, 2], [1, 2]]);
  expect(axisPairs(1)).toEqual([]);
});
test('maps ranks only, includes raw/common support in hover, and excludes null coordinates', () => {
  const traces = benchmarkTraces(result, 0, 1);
  expect(traces).toHaveLength(1);
  expect(traces[0].x).toEqual([.25]);
  expect(traces[0].y).toEqual([.75]);
  expect(traces[0].customdata[0]).toContain(-100);
  expect(benchmarkTraces(result, 0, 1, 'raw')[0].x).toEqual([0]);
});
test('display discordance thresholds are numeric criteria, never impute null to zero', () => {
  expect(discordancePairs(result, .2, 2)).toHaveLength(1);
  expect(discordancePairs(result, .05, 2)).toHaveLength(0);
});
