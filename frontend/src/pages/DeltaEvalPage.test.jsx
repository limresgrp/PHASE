import { fireEvent, render, screen, waitFor } from '@testing-library/react';
import { MemoryRouter, Route, Routes } from 'react-router-dom';
import DeltaEvalPage from './DeltaEvalPage';
import * as projects from '../api/projects';
import * as jobs from '../api/jobs';

jest.mock('../api/projects');
jest.mock('../api/jobs');
jest.mock('react-plotly.js', () => () => null);
jest.mock('../components/common/HelpDrawer', () => () => null);
jest.mock('../components/common/EnergyDistributionPlot', () => ({
  __esModule: true, default: () => null, EnergySeriesSelectorButton: () => null,
  buildEnergyDistributionPlot: () => ({ data: [], layout: {} }), pickEnergyColor: () => '#000000',
  useEnergySeriesSelection: () => ({ selectedSeries: [], selectedIds: [], selectedSeriesIds: [] }),
}));

test('delta energy advanced dialog queues all unique pairs with the same sample and frame limits', async () => {
  const models = ['a', 'b', 'c'].map((id) => ({ model_id: id, name: `Model ${id}`, params: { fit_mode: 'delta' } }));
  projects.fetchSystem.mockResolvedValue({ name: 'System', metastable_clusters: [{ cluster_id: 'c', path: 'cluster.npz', potts_models: models,
    samples: [{ sample_id: 'md1', name: 'MD one', type: 'md_eval', n_frames: 100 }] }] });
  projects.fetchClusterAnalyses.mockResolvedValue({ analyses: [] });
  projects.fetchClusterUiSetups.mockResolvedValue({ setups: [] });
  projects.fetchPottsClusterInfo.mockResolvedValue({});
  jobs.submitDeltaEnergyJob.mockImplementation(async (p) => ({ job_id: `${p.model_a_id}-${p.model_b_id}` }));
  jobs.fetchJobStatus.mockImplementation(() => new Promise(() => {}));
  render(<MemoryRouter initialEntries={['/projects/p/systems/s/sampling/delta_eval?cluster_id=c']} future={{ v7_startTransition: true, v7_relativeSplatPath: true }}>
    <Routes><Route path="/projects/:projectId/systems/:systemId/sampling/delta_eval" element={<DeltaEvalPage />} /></Routes>
  </MemoryRouter>);
  fireEvent.click(await screen.findByRole('button', { name: /New analysis/ }));
  fireEvent.click(screen.getByText('Advanced: compare multiple models'));
  fireEvent.click(screen.getByLabelText('Run all unique model pairs'));
  for (const id of ['a', 'b', 'c']) fireEvent.click(screen.getByLabelText(`Model ${id}`));
  expect(screen.getByText('3 independent analyses')).toBeInTheDocument();
  // The old sidebar run form is hidden; submit the visible modal button.
  fireEvent.click(screen.getAllByRole('button', { name: 'Run analysis' })[0]);
  await waitFor(() => expect(jobs.submitDeltaEnergyJob).toHaveBeenCalledTimes(3));
  const bodies = jobs.submitDeltaEnergyJob.mock.calls.map(([body]) => body);
  expect(bodies.map((p) => [p.model_a_id, p.model_b_id])).toEqual([['a', 'b'], ['a', 'c'], ['b', 'c']]);
  expect(bodies.every((p) => p.sample_ids.join(',') === 'md1' && Object.keys(p.frame_limits).length === 0)).toBe(true);
});
