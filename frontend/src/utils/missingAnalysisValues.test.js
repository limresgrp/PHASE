import { preserveMissingAnalysisValues } from './missingAnalysisValues';

test('missing scientific values cannot become artificial zero JS or commitment', () => {
  const result = preserveMissingAnalysisValues({ metadata: { optional: null }, data: {
    js_node_a: [[0, null, 0.3]], q_residue_all: [[null, 0.5]], node_js: [null, 0], sample_ids: ['pose'],
  } });
  expect(Number.isNaN(Number(result.data.js_node_a[0][1]))).toBe(true);
  expect(result.data.js_node_a[0][0]).toBe(0);
  expect(result.data.q_residue_all[0][1]).toBe(0.5);
  expect(result.metadata.optional).toBeNull();
  expect(result.data.sample_ids).toEqual(['pose']);
  const framewise = preserveMissingAnalysisValues({ node: { sym: [[null, 0.1]] }, edge: { pol: [[null]] } });
  expect(Number.isNaN(framewise.node.sym[0][0])).toBe(true);
  expect(Number.isNaN(framewise.edge.pol[0][0])).toBe(true);
});
