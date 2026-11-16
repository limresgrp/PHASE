import { fireEvent, render, screen } from '@testing-library/react';
import { StateCard } from './StateCards';

jest.mock('./ReferenceStructurePanel', () => ({ AlignmentDetails: () => null }));

test('Build always builds, even after a trajectory-only upload choice', () => {
  const onUpload = jest.fn();
  render(<StateCard state={{ state_id: 'pose', name: 'Pose', pdb_file: 'pose.pdb' }} onUpload={onUpload} />);
  const fileInput = screen.getByLabelText('Upload/save trajectory file');
  fireEvent.change(fileInput, { target: { files: [new File(['traj'], 'traj.xtc')] } });
  fireEvent.click(screen.getByLabelText(/Rebuild descriptors after upload/));
  expect(screen.getByText('Upload only')).toBeInTheDocument();
  fireEvent.change(fileInput, { target: { files: [] } });
  expect(screen.getByLabelText(/Rebuild descriptors after upload/)).toBeChecked();
  fireEvent.click(screen.getByText('Build'));
  expect(onUpload.mock.calls[0][5]).toBe(true);
});
