import { render, screen } from '@testing-library/react';
import { MemoryRouter, Route, Routes } from 'react-router-dom';
import SimulationResultPage from './SimulationResultPage';
import { fetchResult } from '../api/jobs';

jest.mock('../api/jobs');

function renderResult() {
  return render(<MemoryRouter initialEntries={['/result/job']} future={{ v7_startTransition: true, v7_relativeSplatPath: true }}>
    <Routes>
      <Route path="/result/:jobId" element={<SimulationResultPage />} />
      <Route path="/projects/p/systems/s/sampling/visualize" element={<p>Sampling destination</p>} />
    </Routes>
  </MemoryRouter>);
}

test('redirects completed sampling jobs to their system explorer', async () => {
  fetchResult.mockResolvedValue({ analysis_type: 'simulation', system_reference: {
    project_id: 'p', system_id: 's', cluster_id: 'c', sample_id: 'sample',
  } });
  renderResult();
  expect(await screen.findByText('Sampling destination')).toBeInTheDocument();
});

test('handles results without a system reference without navigating or crashing', async () => {
  fetchResult.mockResolvedValue({ analysis_type: 'simulation' });
  renderResult();
  expect(await screen.findByRole('button', { name: 'Open Sampling Explorer' })).toBeDisabled();
});
