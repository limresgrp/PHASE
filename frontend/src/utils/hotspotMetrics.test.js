import { hotspotRows, hotspotColor } from './hotspotMetrics';

const data = { residue_keys: ['A', 'B', 'C'], edges: [[0, 1], [1, 2]],
  node_separation: [1, -2, null], edge_separation: [.5, -.8],
  node_parameter_rms: [1, 2, 3], edge_parameter_rms: [.5, .8],
  node_ci_low: [.1, -3, null], node_ci_high: [2, -.1, null], edge_ci_low: [-.1, -1], edge_ci_high: [.9, -.2],
  node_preference: [[-.5, .3, null]], edge_preference: [[.4, -.2]],
};

test('null observations stay unavailable; absolute effects are ranked and identities preserved', () => {
  const result = hotspotRows(data, { labels: { B: 'TYR193' } });
  expect(result.nodes.map((r) => r.index)).toEqual([1, 0]);
  expect(result.nodes[0].label).toBe('TYR193');
  expect(result.edges[0].indices).toEqual([1, 2]);
});
test('CI, sign, minimum, percentile and selection filters operate on individual effects', () => {
  expect(hotspotRows(data, { reliable: true }).edges.map((r) => r.index)).toEqual([1]);
  expect(hotspotRows(data, { direction: 'positive' }).nodes.map((r) => r.index)).toEqual([0]);
  expect(hotspotRows(data, { minimum: 1.5 }).nodes.map((r) => r.index)).toEqual([1]);
  expect(hotspotRows(data, { selection: [1] }).edges).toHaveLength(0);
  expect(hotspotRows(data, { selection: [2], selectionMode: 'include' }).edges).toHaveLength(1);
  expect(hotspotRows(data, { percentile: 50 }).nodes).toHaveLength(1);
});
test('blending is display only; all incident edges contribute with documented mean or sum', () => {
  const original = hotspotRows(data).nodes[0];
  const mean = hotspotRows(data, { blend: true, alpha: 1 }).nodes[0];
  const sum = hotspotRows(data, { blend: true, alpha: 1, aggregation: 'sum' }).nodes[0];
  expect(mean.effect).toBe(original.effect);
  expect(mean.color).toBeCloseTo(.65);
  expect(sum.color).toBeCloseTo(1.3);
});
test('preference is independently signed, with red positive and blue negative', () => {
  const rows = hotspotRows(data, { view: 'preference' });
  expect(rows.nodes.find((r) => r.index === 0).color).toBe(-.5);
  expect(hotspotColor(1, 1, true)).toBe(0xd22d37);
  expect(hotspotColor(-1, 1, true)).toBe(0x2d69d2);
});
