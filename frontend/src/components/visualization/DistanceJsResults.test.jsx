import { fireEvent, render, screen } from '@testing-library/react';
import DistanceJsResults from './DistanceJsResults';

jest.mock('react-plotly.js', () => () => <div data-testid="plot" />);
jest.mock('./DistanceJsStructure', () => () => <div>Embedded structure viewer</div>);

const result = { metadata: { distance_available: true, thresholds: [6], condition_a: { state_ids: ['a'] }, condition_b: { state_ids: ['b'] } },
  data: { edge_global_js: [.1], edge_bin_js: [[0], [.5]], edge_bin_counts_a: [[40], [40]], edge_bin_counts_b: [[40], [40]],
    edge_bin_occupancy_a: [[.5], [.5]], edge_bin_occupancy_b: [[.5], [.5]],
    edge_labels: ['res_1 - res_2'], edge_index: [[0], [1]], residue_keys: ['res_1', 'res_2'], distance_bin_labels: ['<6', '>=6'] } };

test('renders supported distance effects with names, profile and 3D tab', () => {
  render(<DistanceJsResults result={result} residueNames={{ res_1: 'SER', res_2: 'TYR' }} />);
  expect(screen.getByRole('button', { name: 'SER1 — TYR2' })).toBeInTheDocument();
  expect(screen.getByText('SER1')).toBeInTheDocument();
  fireEvent.click(screen.getByText('3D structure'));
  expect(screen.getByText('Embedded structure viewer')).toBeInTheDocument();
});

test('cutoff removes complete bins without pretending a single-bin JS range exists', () => {
  render(<DistanceJsResults result={result} />);
  fireEvent.change(screen.getByLabelText('Use distances below'), { target: { value: '6' } });
  expect(screen.getByRole('status')).toHaveTextContent('No pairs pass');
});

test('raw display preserves normalized range filters', () => {
  render(<DistanceJsResults result={result} />);
  fireEvent.change(screen.getByLabelText('Minimum JS range (0–1)'), { target: { value: '.7' } });
  expect(screen.getByRole('button', { name: 'res_1 — res_2' })).toBeInTheDocument();
  fireEvent.click(screen.getByLabelText('Show normalized JS'));
  expect(screen.getByRole('button', { name: 'res_1 — res_2' })).toBeInTheDocument();
  fireEvent.change(screen.getByLabelText('Minimum JS range (0–1)'), { target: { value: '.8' } });
  expect(screen.getByRole('status')).toHaveTextContent('No pairs pass');
});
