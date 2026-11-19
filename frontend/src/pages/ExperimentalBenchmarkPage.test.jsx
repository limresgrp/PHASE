import { fireEvent, render, screen, waitFor, within } from '@testing-library/react';
import { MemoryRouter, Route, Routes } from 'react-router-dom';
import ExperimentalBenchmarkPage from './ExperimentalBenchmarkPage';
import * as projects from '../api/projects';
import * as api from '../api/client';
import * as jobs from '../api/jobs';

jest.mock('../api/projects');
jest.mock('../api/client');
jest.mock('../api/jobs');
jest.mock('../components/common/HelpDrawer', () => () => null);
jest.mock('../components/visualization/ExperimentalBenchmarkResults', () => ({ result }) => <div data-testid="results">{result.metadata.name}</div>);

const choices = {
  states: ['s1', 's2'].map((state_id) => ({ state_id, name: state_id, n_frames: 1, descriptor_ready: true, assigned_samples: [{ sample_id: state_id, name: state_id }] })),
  models: [{ model_id: 'a', name: 'Model A' }, { model_id: 'b', name: 'Model B' }], md_samples: [],
  label_schemes: [{ scheme_id: 'labels', name: 'User labels', kind: 'ordinal', revision: 1, categories: ['arbitrary-low', 'arbitrary-high'], labels: { s1: 'arbitrary-low', s2: 'arbitrary-high' } }],
};

beforeEach(() => {
  jest.clearAllMocks();
  projects.fetchSystem.mockResolvedValue({ metastable_clusters: [{ cluster_id: 'c', path: 'cluster.npz', name: 'Cluster' }] });
  projects.fetchClusterAnalyses.mockResolvedValue({ analyses: [{ analysis_id: 'saved', name: 'Saved benchmark', state_ids: ['s1', 's2'], parameters: { contrasts: [{}] } }] });
  api.requestJSON.mockImplementation((url) => Promise.resolve(url.includes('/context') ? choices : url.includes('/analyses/saved') ? { metadata: { name: 'Immutable benchmark' } } : { job_id: 'queued' }));
  jobs.fetchJobStatus.mockImplementation(() => new Promise(() => {}));
});
function page(path = '/projects/p/systems/s/experimental_benchmarks?cluster_id=c') {
  return render(<MemoryRouter initialEntries={[path]} future={{ v7_startTransition: true, v7_relativeSplatPath: true }}><Routes>
    <Route path="/projects/:projectId/systems/:systemId/experimental_benchmarks" element={<ExperimentalBenchmarkPage />} />
    <Route path="/projects/:projectId/systems/:systemId/experimental_benchmarks/:analysisId" element={<ExperimentalBenchmarkPage />} />
  </Routes></MemoryRouter>);
}
test('lists existing benchmarks and opens dedicated result route with download', async () => {
  page(); fireEvent.click(await screen.findByText('Saved benchmark'));
  expect(await screen.findByTestId('results')).toHaveTextContent('Immutable benchmark');
  expect(screen.getByRole('link', { name: /Download CSV/ })).toHaveAttribute('href', expect.stringMatching(/analyses\/saved\/download/));
});
test('run hyperparameters stay in the popup and submit a background benchmark', async () => {
  page(); await screen.findByText('Saved benchmark');
  await waitFor(() => expect(api.requestJSON).toHaveBeenCalledWith(expect.stringMatching(/context/)));
  fireEvent.click(screen.getByRole('button', { name: /New benchmark/ }));
  const dialog = screen.getByRole('dialog');
  expect(within(dialog).getByText('Model contrasts')).toBeInTheDocument();
  fireEvent.click(within(dialog).getByRole('button', { name: 'Run benchmark' }));
  await waitFor(() => expect(api.requestJSON).toHaveBeenCalledWith(expect.stringMatching(/experimental-benchmarks\/analyses$/), expect.objectContaining({ method: 'POST', body: expect.objectContaining({ state_ids: ['s1', 's2'], scheme_id: 'labels', fingerprint_mode: 'standardized' }) })));
});
test('annotation editor uses explicit ordinal ordering and saves labels separately', async () => {
  page(); await screen.findByText('Saved benchmark');
  fireEvent.click(screen.getByRole('button', { name: /Label schemes/ }));
  const dialog = screen.getByRole('dialog');
  fireEvent.change(within(dialog).getByLabelText('Load scheme'), { target: { value: 'labels' } });
  fireEvent.click(within(dialog).getByRole('button', { name: 'Save scheme' }));
  await waitFor(() => expect(api.requestJSON).toHaveBeenCalledWith(expect.stringMatching(/label-schemes\/labels$/), expect.objectContaining({ method: 'PUT', body: expect.objectContaining({ categories: ['arbitrary-low', 'arbitrary-high'], labels: { s1: 'arbitrary-low', s2: 'arbitrary-high' } }) })));
});
