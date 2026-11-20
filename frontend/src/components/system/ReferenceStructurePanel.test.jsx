import { fireEvent, render, screen, waitFor } from '@testing-library/react';
import ReferenceStructurePanel, { AlignmentDetails } from './ReferenceStructurePanel';
import { requestJSON } from '../../api/client';

jest.mock('../../api/client');
jest.mock('../common/HelpDrawer', () => () => null);
const system = { reference_pdb_file: 'reference/reference.pdb', reference_chain_id: 'A', reference_sequence: 'ACDE', reference_residue_keys: ['res_1'], states: {} };

beforeEach(() => jest.clearAllMocks());

test('new uploaded structures require an alignment preview before state creation', async () => {
  requestJSON.mockResolvedValue({ identity: 1, coverage_reference: 1, quality: 'excellent', alignment_status: 'aligned', correspondence: [] });
  render(<ReferenceStructurePanel projectId="p" systemId="s" system={system} onChanged={jest.fn()} />);
  fireEvent.click(screen.getByText('+ Add state / replica'));
  fireEvent.change(screen.getByLabelText('State name'), { target: { value: 'Target' } });
  fireEvent.change(screen.getByLabelText('Topology'), { target: { value: 'upload' } });
  fireEvent.change(screen.getByLabelText('PDB file'), { target: { files: [new File(['ATOM'], 'target.pdb')] } });
  expect(screen.getByText('Create selected state(s)')).toBeDisabled();
  fireEvent.click(screen.getByText('Preview alignment(s)'));
  await waitFor(() => expect(screen.getByText('Create selected state(s)')).toBeEnabled());
  expect(requestJSON.mock.calls[0][0]).toMatch(/reference\/preview$/);
  expect(requestJSON).toHaveBeenCalledTimes(1);
  fireEvent.change(screen.getByLabelText(/Target chain/), { target: { value: 'B' } });
  expect(screen.getByText('Create selected state(s)')).toBeDisabled();
});

test('reference is immutable once states exist and partial details name missing residues', () => {
  render(<><ReferenceStructurePanel projectId="p" systemId="s" system={{ ...system, states: { a: { state_id: 'a' } } }} onChanged={jest.fn()} />
    <AlignmentDetails status="partial" alignment={{ missing_selected_residues: ['res_45A'], quality: 'warning' }} /></>);
  expect(screen.getByText('Change reference / selection')).toBeDisabled();
  expect(screen.getByText(/Missing coordinates \(masked in descriptors and energies\): res_45A/)).toBeInTheDocument();
});

test('a PDB in the reference trajectory field switches to aligned structure upload', async () => {
  requestJSON.mockResolvedValueOnce({ identity: 1, alignment_status: 'aligned', correspondence: [] })
    .mockResolvedValueOnce({ states: {} });
  render(<ReferenceStructurePanel projectId="p" systemId="s" system={system} onChanged={jest.fn()} />);
  fireEvent.click(screen.getByText('+ Add state / replica'));
  fireEvent.change(screen.getByLabelText('State name'), { target: { value: 'Crystal' } });
  const file = new File(['ATOM'], 'crystal.pdb');
  fireEvent.change(screen.getByLabelText('Optional trajectory (must match this topology)'), { target: { files: [file] } });
  expect(screen.getByLabelText('Topology')).toHaveValue('upload');
  expect(screen.getByText('Create selected state(s)')).toBeDisabled();
  fireEvent.click(screen.getByText('Preview alignment(s)'));
  await waitFor(() => expect(screen.getByText('Create selected state(s)')).toBeEnabled());
  expect(requestJSON.mock.calls[0][1].body.get('pdb')).toBe(file);
  fireEvent.click(screen.getByText('Create selected state(s)'));
  await waitFor(() => expect(requestJSON).toHaveBeenCalledTimes(2));
  const payload = requestJSON.mock.calls[1][1].body;
  expect(payload.get('pdb')).toBe(file);
  expect(payload.has('trajectory')).toBe(false);
  expect(payload.has('use_reference')).toBe(false);
});

test('creation builds by default and queues assignment only for newly created states', async () => {
  const onChanged = jest.fn();
  requestJSON.mockResolvedValueOnce({ states: { old: { state_id: 'old', descriptor_file: 'old.npz' }, fresh: { state_id: 'fresh', descriptor_file: 'fresh.npz' } } })
    .mockResolvedValueOnce({ job_id: 'assignment-job' });
  render(<ReferenceStructurePanel projectId="p" systemId="s" system={{ ...system, states: { old: { state_id: 'old' } }, metastable_clusters: [{ cluster_id: 'c', name: 'Cluster One', path: 'cluster.npz' }] }} onChanged={onChanged} />);
  fireEvent.click(screen.getByText('+ Add state / replica'));
  fireEvent.change(screen.getByLabelText('State name'), { target: { value: 'New pose' } });
  expect(screen.getByLabelText('Build descriptors now (PDB-only = one frame)')).toBeChecked();
  fireEvent.click(screen.getByLabelText('Cluster One'));
  fireEvent.click(screen.getByText('Create selected state(s)'));
  await waitFor(() => expect(onChanged).toHaveBeenCalled());
  expect(requestJSON.mock.calls[0][1].body.get('build_descriptors_after_upload')).toBe('true');
  expect(requestJSON.mock.calls[1]).toEqual(['/submit/md_samples_refresh', {
    method: 'POST', body: { project_id: 'p', system_id: 's', cluster_id: 'c', state_ids: ['fresh'], overwrite: true, cleanup: true },
  }]);
  expect(screen.getByText('View assignment progress')).toHaveAttribute('href', '/jobs/assignment-job');
});

test('existing descriptor-ready states can be assigned from System page', async () => {
  requestJSON.mockResolvedValue({ job_id: 'existing-job' });
  render(<ReferenceStructurePanel projectId="p" systemId="s" system={{ ...system,
    states: { pose: { state_id: 'pose', name: 'Crystal pose', descriptor_file: 'pose.npz' }, unbuilt: { state_id: 'unbuilt', name: 'Not ready' } },
    metastable_clusters: [{ cluster_id: 'c', name: 'Cluster One', path: 'cluster.npz' }],
  }} onChanged={jest.fn()} />);
  fireEvent.click(screen.getByText('Assign states to clusters'));
  expect(screen.queryByText('Not ready')).not.toBeInTheDocument();
  fireEvent.click(screen.getByLabelText('Crystal pose'));
  fireEvent.click(screen.getByLabelText('Cluster One'));
  fireEvent.click(screen.getByText('Run assignment'));
  expect(await screen.findByText('View assignment progress')).toBeInTheDocument();
  expect(requestJSON.mock.calls[0][1].body.state_ids).toEqual(['pose']);
});

test('partial PDB alignments allow automatic descriptor building with a missing mask', async () => {
  requestJSON.mockResolvedValueOnce({ quality: 'warning', alignment_status: 'partial', missing_selected_residues: ['res_45A'], correspondence: [] });
  render(<ReferenceStructurePanel projectId="p" systemId="s" system={system} onChanged={jest.fn()} />);
  fireEvent.click(screen.getByText('+ Add state / replica'));
  fireEvent.change(screen.getByLabelText('State name'), { target: { value: 'Partial pose' } });
  fireEvent.change(screen.getByLabelText('Topology'), { target: { value: 'upload' } });
  fireEvent.change(screen.getByLabelText('PDB file'), { target: { files: [new File(['ATOM'], 'partial.pdb')] } });
  fireEvent.click(screen.getByText('Preview alignment(s)'));
  await waitFor(() => expect(screen.getByText('Create selected state(s)')).toBeEnabled());
  expect(screen.getByText(/Available residues will be built and assigned/)).toBeInTheDocument();
  expect(screen.getByLabelText('Build descriptors now (PDB-only = one frame)')).toBeChecked();
});
