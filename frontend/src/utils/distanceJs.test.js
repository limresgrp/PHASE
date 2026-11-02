import { binsBelowCutoff, distanceJsResidues, distanceJsRows, JS_MAX } from './distanceJs';

const data = {
  edge_global_js: [.1, .2], distance_bin_labels: ['near', 'far'],
  edge_bin_js: [[0, .2], [JS_MAX, .3]],
  edge_bin_counts_a: [[30, 30], [30, 1]], edge_bin_counts_b: [[30, 30], [30, 30]],
  edge_index: [[0, 1], [1, 2]], residue_keys: ['res_1', 'res_2', 'res_3'],
};
test('ranks changes only across two supported bins, not global-vs-bin JS', () => {
  const rows = distanceJsRows(data);
  expect(rows).toHaveLength(1); expect(rows[0].change).toBeCloseTo(1);
  expect(distanceJsResidues(rows, data).map((r) => r.key)).toEqual(['res_1', 'res_2']);
  expect(distanceJsRows(data, { bins: [0] })).toEqual([]);
});
test('null is missing, never zero; legacy single-bin files have no change', () => {
  expect(distanceJsRows({ ...data, edge_bin_js: [[null, null], [null, null]] })).toEqual([]);
  expect(distanceJsRows(data, { sort: 'global' })).toHaveLength(2);
});
test('distance cutoffs retain only complete physical bins', () => {
  expect(binsBelowCutoff({ distance_available: true, thresholds: [6, 8, 10] }, 4, '8')).toEqual([0, 1]);
  expect(binsBelowCutoff({ distance_available: false }, 1, '8')).toEqual([]);
  expect(binsBelowCutoff({}, 2, '')).toEqual([0, 1]);
});
