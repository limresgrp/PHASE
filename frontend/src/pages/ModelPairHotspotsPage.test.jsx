import { fireEvent, render, screen, waitFor } from '@testing-library/react';
import { MemoryRouter, Route, Routes } from 'react-router-dom';
import ModelPairHotspotsPage from './ModelPairHotspotsPage';
import * as projects from '../api/projects';
import * as client from '../api/client';
import * as jobs from '../api/jobs';

jest.mock('../api/projects');
jest.mock('../api/client');
jest.mock('../api/jobs');
jest.mock('../components/common/HelpDrawer', () => () => null);
jest.mock('../components/visualization/HotspotStructure', () => (props) => <div data-testid="network">{props.nodes.length} residues, {props.edges.length} edges</div>);

beforeEach(() => {
  jest.clearAllMocks();
  projects.fetchSystem.mockResolvedValue({ metastable_clusters: [{ cluster_id: 'c', name: 'Cluster', path: 'cluster.npz' }],
    reference_residue_records: [{ residue_key: 'res_2_SER', resname: 'SER', resid: 2 }] });
  projects.fetchClusterUiSetups.mockResolvedValue({ setups: [] });
  projects.fetchClusterAnalyses.mockResolvedValue({ analyses: [{ analysis_id: 'result', name: 'Active / inactive', mode: 'ensemble' }] });
  projects.fetchClusterAnalysisData.mockResolvedValue({ metadata: { name: 'Active / inactive', mode: 'ensemble', samples: [] }, data: {
    residue_keys: ['res_2_SER'], edges: [], node_separation: [2], node_parameter_rms: [1], node_ci_low: [.2], node_ci_high: [3], node_counts_a: [20], node_counts_b: [30], edge_separation: [],
  } });
  client.requestJSON.mockImplementation((url) => Promise.resolve(url.endsWith('/context') ? {
    models: [{ model_id: 'a', name: 'Active', suggested_references: ['sa'] }, { model_id: 'b', name: 'Inactive', suggested_references: ['sb'] }],
    samples: [{ sample_id: 'sa', name: 'Active MD', type: 'md_eval' }, { sample_id: 'sb', name: 'Inactive MD', type: 'md_eval' }],
  } : { job_id: 'queued' }));
  jobs.fetchJobStatus.mockImplementation(() => new Promise(() => {}));
});
const page = () => render(<MemoryRouter initialEntries={['/projects/p/systems/s/sampling/model_pair_hotspots?cluster_id=c']} future={{ v7_startTransition: true, v7_relativeSplatPath: true }}>
  <Routes><Route path="/projects/:projectId/systems/:systemId/sampling/model_pair_hotspots" element={<ModelPairHotspotsPage />} /></Routes></MemoryRouter>);

test('loads all analyses and rank table, links residue to embedded 3D', async () => {
  page();
  const residue = await screen.findByRole('button', { name: 'SER2' });
  expect(screen.getByText('Residue terms')).toBeInTheDocument();
  expect(projects.fetchClusterAnalyses).toHaveBeenCalledWith('p', 's', 'c', { analysisType: 'model_pair_hotspots' });
  fireEvent.click(residue);
  expect(await screen.findByTestId('network')).toHaveTextContent('1 residues, 0 edges');
});

test('new-analysis dialog submits independent background job with suggested editable references', async () => {
  page(); await screen.findByRole('button', { name: 'SER2' });
  fireEvent.click(screen.getByRole('button', { name: /New analysis/ }));
  expect(screen.getByRole('heading', { name: 'New discriminative hotspot analysis' })).toBeInTheDocument();
  fireEvent.click(screen.getByRole('button', { name: 'Run analysis' }));
  await waitFor(() => expect(client.requestJSON).toHaveBeenCalledWith(expect.stringMatching(/hotspots\/analyses$/), expect.objectContaining({
    method: 'POST', body: expect.objectContaining({ model_a_id: 'a', model_b_id: 'b', reference_a: ['sa'], reference_b: ['sb'], block_length: 50 }),
  })));
});

test('metadata opens separately and analysis deletion is wired', async () => {
  window.confirm = jest.fn(() => true);
  projects.deleteClusterAnalysis.mockResolvedValue({ status: 'deleted' });
  page(); await screen.findByRole('button', { name: 'SER2' });
  fireEvent.click(screen.getByTitle('Analysis metadata'));
  expect(screen.getByLabelText('Close metadata')).toBeInTheDocument();
  fireEvent.click(screen.getByLabelText('Close metadata'));
  fireEvent.click(screen.getByTitle('Delete analysis'));
  await waitFor(() => expect(projects.deleteClusterAnalysis).toHaveBeenCalledWith('p', 's', 'c', 'model_pair_hotspots', 'result'));
});

test('advanced selection submits all three pairs with model-specific references', async () => {
  client.requestJSON.mockImplementation((url) => Promise.resolve(url.endsWith('/context') ? {
    models: ['a', 'b', 'c'].map((id) => ({ model_id: id, name: `Model ${id}`, suggested_references: [`md-${id}`] })),
    samples: ['a', 'b', 'c'].map((id) => ({ sample_id: `md-${id}`, name: `MD ${id}`, type: 'md_eval' })),
  } : { job_id: 'queued' }));
  page(); await screen.findByRole('button', { name: 'SER2' });
  fireEvent.click(screen.getByRole('button', { name: /New analysis/ }));
  fireEvent.click(screen.getByText('Advanced: compare multiple models'));
  fireEvent.click(screen.getByLabelText('Run all unique model pairs'));
  for (const id of ['a', 'b', 'c']) fireEvent.click(screen.getByLabelText(`Model ${id}`));
  expect(screen.getByText('3 independent analyses')).toBeInTheDocument();
  fireEvent.click(screen.getByRole('button', { name: 'Run analysis' }));
  await waitFor(() => expect(client.requestJSON.mock.calls.filter(([url]) => url.endsWith('/analyses'))).toHaveLength(3));
  const bodies = client.requestJSON.mock.calls.filter(([url]) => url.endsWith('/analyses')).map(([, options]) => options.body);
  expect(bodies.map((p) => [p.model_a_id, p.model_b_id, p.reference_a, p.reference_b])).toEqual([
    ['a', 'b', ['md-a'], ['md-b']], ['a', 'c', ['md-a'], ['md-c']], ['b', 'c', ['md-b'], ['md-c']],
  ]);
});
