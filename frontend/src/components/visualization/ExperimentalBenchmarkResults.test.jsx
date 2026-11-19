import { render, screen, fireEvent } from '@testing-library/react';
import ExperimentalBenchmarkResults from './ExperimentalBenchmarkResults';

jest.mock('react-plotly.js', () => ({ data }) => <div data-testid="plot">{JSON.stringify(data)}</div>);
const result = {
  metadata: { name: 'Benchmark', parameters: {}, models: [], warnings: [] }, disclaimer: 'Quantitative evidence only',
  label_scheme: { name: 'Categorical', kind: 'categorical', categories: ['x'], labels: { s: 'x', t: 'x' } },
  axes: [{ axis: 0, name: 'A / B' }, { axis: 1, name: 'A / C' }],
  states: [{ state_id: 's', name: 'One' }, { state_id: 't', name: 'Two' }],
  scores: ['s', 't'].flatMap((state_id, i) => [0, 1].map((axis) => ({ state_id, state_name: state_id, label: 'x', axis,
    raw_score: 100 + i, common_score: 10 + i, raw_rank: i, common_rank: i }))),
  statistics: [{ axis: 0, variant: 'common', auc: { auc: .25, n: 2, positive: 'x', negative: 'y' } }],
  axis_redundancy: { rank_spearman: { matrix: [[1, 1], [1, 1]], n: [[2, 2], [2, 2]] } },
  axis_algebra: [], pcoa: { available: false, reason: 'Undefined correlation distance' },
  distance_matrices: { euclidean: [[0, 1], [1, 0]], correlation: [[0, null], [null, 0]] },
  global_threshold: .15, local_threshold: 1,
  discordance: [{ state_a: 's', state_b: 't', global_rank_difference: .1, local_distance: 2 }],
  hotspots: [{ nodes: [{ identity: 'SER2', score: 1 }], edges: [] }], hotspot_comparisons: [], coverage_qc: [],
};
test('renders all quantitative result tabs and handles undefined PCoA without crashing', () => {
  render(<ExperimentalBenchmarkResults result={result} />);
  expect(screen.getByText(/These are not probabilities/)).toBeInTheDocument();
  fireEvent.click(screen.getByRole('button', { name: 'Label statistics' }));
  expect(screen.getByText(/monotonic ordering, not linearity/)).toBeInTheDocument();
  fireEvent.click(screen.getByRole('button', { name: 'Axis relationships' }));
  expect(screen.getByText('Identical-support axis algebra')).toBeInTheDocument();
  fireEvent.click(screen.getByRole('button', { name: 'Local fingerprints / PCoA' }));
  expect(screen.getByText('Undefined correlation distance')).toBeInTheDocument();
  fireEvent.click(screen.getByRole('button', { name: 'Global / local differences' }));
  expect(screen.getByText('Pairs matching the supplied thresholds')).toBeInTheDocument();
  fireEvent.click(screen.getByRole('button', { name: 'Hotspots' }));
  expect(screen.getByText('SER2')).toBeInTheDocument();
  fireEvent.click(screen.getByRole('button', { name: 'Coverage QC' }));
  expect(screen.getByText(/not observed zeros/)).toBeInTheDocument();
});

test('unsupported common ranks explain the missing map instead of rendering an empty plot', () => {
  render(<ExperimentalBenchmarkResults result={{ ...result, scores: result.scores.map((s) => ({ ...s, common_rank: null })) }} />);
  expect(screen.getByText(/No supported ranks for this axis pair/)).toBeInTheDocument();
});
