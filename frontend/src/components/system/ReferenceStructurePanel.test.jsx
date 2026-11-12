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
  expect(screen.getByText(/Descriptors cannot be built: missing res_45A/)).toBeInTheDocument();
});
