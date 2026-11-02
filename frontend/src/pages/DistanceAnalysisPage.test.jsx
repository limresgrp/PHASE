import { fireEvent, render, screen, waitFor, within } from '@testing-library/react';
import { MemoryRouter, Route, Routes } from 'react-router-dom';
import DistanceAnalysisPage from './DistanceAnalysisPage';
import * as api from '../api/distances';
import * as projects from '../api/projects';
import * as jobs from '../api/jobs';

jest.mock('react-plotly.js', () => () => <div data-testid="plot" />);
jest.mock('../components/common/HelpDrawer', () => () => null);
jest.mock('../components/visualization/DistanceJsResults', () => () => <div />);
jest.mock('../api/distances');
jest.mock('../api/projects');
jest.mock('../api/jobs');

beforeEach(() => {
  jest.clearAllMocks();
  projects.fetchSystem.mockResolvedValue({ metastable_clusters: [{ cluster_id: 'c', name: 'Toy cluster', path: 'cluster.npz' }] });
  projects.fetchClusterAnalyses.mockResolvedValue({ analyses: [] });
  api.fetchDistanceStatus.mockResolvedValue({ distance_available: false, pair_count: 1, metadata: {}, states: [
    { state_id: 'a', name: 'Active', frame_count: 20 }, { state_id: 'b', name: 'Inactive', frame_count: 30 },
  ], metastable_ids: [] });
  api.fetchDistancePresets.mockResolvedValue({ presets: [] });
  api.fetchDistancePairs.mockResolvedValue({ pairs: [], total: 0 });
  api.submitDistanceJs.mockResolvedValue({ job_id: 'queued-js' });
  jobs.fetchJobStatus.mockImplementation(() => new Promise(() => {}));
});

function page() {
  return render(<MemoryRouter initialEntries={['/projects/p/systems/s/sampling/distance?cluster_id=c']} future={{ v7_startTransition: true, v7_relativeSplatPath: true }}>
    <Routes><Route path="/projects/:projectId/systems/:systemId/sampling/distance" element={<DistanceAnalysisPage />} /></Routes>
  </MemoryRouter>);
}

test('requires a cache build when neither candidate nor contact edges exist', async () => {
  const status = await api.fetchDistanceStatus();
  api.fetchDistanceStatus.mockResolvedValue({ ...status, pair_count: 0 });
  page();
  await screen.findByText(/No residue pairs are available\. Click/);
  expect(screen.getByText('New conditional-JS analysis')).toBeDisabled();
  expect(screen.getByRole('button', { name: 'Build / refresh cache' })).toBeEnabled();
  expect(api.submitDistanceJs).not.toHaveBeenCalled();
});

test('starts inspection on one state and keeps draft trajectory changes until refresh', async () => {
  page();
  await screen.findByText(/No distance cache: each pair uses one/);
  await waitFor(() => expect(api.fetchDistancePairs).toHaveBeenCalled());
  expect(api.fetchDistancePairs.mock.calls.at(-1)[3]).toMatchObject({ state_ids: 'a' });
  fireEvent.click(screen.getByText('Select all trajectories'));
  expect(api.fetchDistancePairs.mock.calls.at(-1)[3]).toMatchObject({ state_ids: 'a' });
  fireEvent.click(screen.getByText('Refresh inspection'));
  await waitFor(() => expect(api.fetchDistancePairs.mock.calls.at(-1)[3]).toMatchObject({ state_ids: 'a,b' }));
});

test('submits conditions through the background job dialog even without a geometry cache', async () => {
  page();
  await screen.findByText(/No distance cache: each pair uses one/);
  await waitFor(() => expect(screen.getByText('New conditional-JS analysis')).toBeEnabled());
  fireEvent.click(screen.getByText('New conditional-JS analysis'));
  const dialog = screen.getByRole('dialog');
  expect(within(dialog).getByText(/requested bins collapse to All distances/)).toBeInTheDocument();
  fireEvent.click(within(dialog).getByText('Run background job'));
  await waitFor(() => expect(api.submitDistanceJs).toHaveBeenCalledWith('p', 's', 'c', expect.objectContaining({
    condition_a: { state_ids: ['a'] }, condition_b: { state_ids: ['b'] }, thresholds: [6, 8, 10], min_bin_frames: 1,
  })));
  await waitFor(() => expect(screen.queryByRole('dialog')).not.toBeInTheDocument());
});

test('applies the physical contact filter only on explicit inspection refresh', async () => {
  const status = await api.fetchDistanceStatus();
  api.fetchDistanceStatus.mockResolvedValue({ ...status, distance_available: true, metadata: { max_distance_cutoff: 12 } });
  page();
  const filter = await screen.findByLabelText('Only pairs contacting below the cutoff');
  await waitFor(() => expect(filter).toBeEnabled());
  fireEvent.click(filter);
  expect(api.fetchDistancePairs.mock.calls.at(-1)[3].min_contact_fraction).toBeUndefined();
  fireEvent.click(screen.getByText('Refresh inspection'));
  await waitFor(() => expect(api.fetchDistancePairs.mock.calls.at(-1)[3]).toMatchObject({ contact_cutoff: 6, min_contact_fraction: .1 }));
});
