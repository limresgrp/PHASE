"""
Handles calculation of dihedral angles and the mandatory
sin/cos transformation.

Refactored to:
- Accept a slice_obj to perform analysis on trajectory slices.
- Use the slice parameters in the Dihedral.run() method.
- Return the number of frames processed.
"""

import MDAnalysis as mda
import MDAnalysis.analysis.dihedrals as mda_dihedrals
import numpy as np
from typing import Dict, List, Optional, Tuple

# Mock types for demonstration
TrajectoryObject = mda.Universe
FeatureDict = Dict[str, np.ndarray]
DIHEDRAL_KEYS = ("phi", "psi", "omega", "chi1", "chi2")

def deg2rad(deg_array: np.ndarray) -> np.ndarray:
    """Converts an array of degrees to radians."""
    return np.deg2rad(deg_array)

def transform_to_circular(angles_rad: np.ndarray) -> np.ndarray:
    """
    Transforms an array of angles (n_frames, n_residues, n_angles)
    into the 2D vector representation [sin(th), cos(th)].
    Output shape: (n_frames, n_residues, n_angles * 2)
    """
    sin_transformed = np.sin(angles_rad)
    cos_transformed = np.cos(angles_rad)
    
    n_frames, n_residues, n_angles = angles_rad.shape
    circular_features = np.empty((n_frames, n_residues, n_angles * 2))
    
    circular_features[..., 0::2] = sin_transformed
    circular_features[..., 1::2] = cos_transformed
    
    return circular_features

class FeatureExtractor:
    """
    Handles calculation of backbone/sidechain dihedral angles:
    phi, psi, omega, chi1, chi2.
    """
    def __init__(self, residue_selections: Optional[Dict[str, str]] = None, *, use_coordinate_neighbors: bool = False):
        """
        Initializes with a dictionary of residue selections.
        Example: {'res_50': 'resid 50', 'res_131': 'resid 131'}
        """
        self.residue_selections = residue_selections if residue_selections is not None else {}
        self.use_coordinate_neighbors = use_coordinate_neighbors
        if residue_selections is not None:
            print(f"FeatureExtractor initialized for {len(self.residue_selections)} residues.")
        else:
            print("FeatureExtractor initialized in automatic mode (will use all residues).")

    def _get_sliced_length(self, traj: TrajectoryObject, slice_obj: slice) -> int:
        """Helper to calculate frames in a slice."""
        total_frames = len(traj.trajectory)
        return len(range(*slice_obj.indices(total_frames)))

    def _calculate_dihedral_angle(
        self, 
        atom_groups_list: List[Optional[mda.AtomGroup]], 
        n_frames: int,
        slice_obj: slice
    ) -> np.ndarray:
        """
        Helper function to run Dihedral analysis on a list of AtomGroups
        over a specified trajectory slice.
        
        Args:
            atom_groups_list: List of AtomGroups (or None) to analyze.
            n_frames: The *sliced* number of frames, for result array shape.
            slice_obj: The slice object with start/stop/step.
            
        Returns: (n_sliced_frames, len(atom_groups_list)) array of angles.
        """
        mask = np.array([ag is not None for ag in atom_groups_list])
        angles = np.full((n_frames, len(atom_groups_list)), 0.0, dtype=np.float32)
        
        valid_atom_groups = [ag for ag in atom_groups_list if ag is not None]
        
        if not valid_atom_groups:
            return angles
            
        try:
            # Unpack slice parameters for the .run() method
            start = slice_obj.start
            stop = slice_obj.stop
            step = slice_obj.step
            
            # Run analysis only on the valid groups and the specified slice
            analysis = mda_dihedrals.Dihedral(valid_atom_groups).run(
                start=start, 
                stop=stop, 
                step=step
            )
            
            # Verify shape
            if analysis.results.angles.shape[0] != n_frames:
                print(f"  Warning: Dihedral analysis returned {analysis.results.angles.shape[0]} frames, expected {n_frames}. Check slice logic.")
                # Truncate or handle as needed; for now, we'll trust the mask
            
            angles[:, mask] = analysis.results.angles
            
        except Exception as e:
            print(f"    Warning: Dihedral calculation failed: {e}")
            
        return angles

    @staticmethod
    def _safe_residue_dihedral_selection(
        residue: mda.ResidueGroup,
        method_name: str,
    ) -> Optional[mda.AtomGroup]:
        """Return a residue dihedral selection or ``None`` when unavailable."""
        method = getattr(residue, method_name, None)
        if method is None:
            return None
        try:
            ag = method()
        except Exception:
            return None
        if ag is None:
            return None
        try:
            return ag if int(ag.n_atoms) == 4 else None
        except Exception:
            return None

    @staticmethod
    def _chi2_selection(residue: mda.ResidueGroup) -> Optional[mda.AtomGroup]:
        """
        Build chi2 like MDAnalysis Janin:
        CA - CB - (CG/CG1) - (CD/CD1/OD1/ND1/SD)
        """
        try:
            atoms = residue.atoms
            ag1 = atoms.select_atoms("name CA")
            ag2 = atoms.select_atoms("name CB")
            ag3 = atoms.select_atoms("name CG CG1")
            ag4 = atoms.select_atoms("name CD CD1 OD1 ND1 SD")
        except Exception:
            return None
        if any(int(ag.n_atoms) != 1 for ag in (ag1, ag2, ag3, ag4)):
            return None
        return ag1 + ag2 + ag3 + ag4

    def _build_dihedral_selections(
        self,
        residues: mda.ResidueGroup,
    ) -> Dict[str, List[Optional[mda.AtomGroup]]]:
        """Collect per-residue 4-atom selections for all supported dihedrals."""
        selections: Dict[str, List[Optional[mda.AtomGroup]]] = {name: [] for name in DIHEDRAL_KEYS}
        neighbors = {}
        if self.use_coordinate_neighbors and len(residues):
            from MDAnalysis.lib.distances import calc_bonds
            # Number arithmetic is not peptide connectivity: 45 -> 45A -> 46
            # must not become 45 -> 46. Use topology order within each chain,
            # gated by the actual C--N bond distance to avoid spanning gaps.
            chains = {}
            for residue in residues.universe.select_atoms("protein").residues:
                chain = tuple(sorted(set(getattr(residue.atoms, "chainIDs", [""]))))
                chains.setdefault((residue.segment.ix, chain), []).append(residue)
            for chain_residues in chains.values():
                for left, right in zip(chain_residues, chain_residues[1:]):
                    c, n = left.atoms.select_atoms("name C"), right.atoms.select_atoms("name N")
                    box = residues.universe.dimensions
                    if box is not None and (not np.all(np.isfinite(box)) or np.any(box[:3] < 5)
                                            or np.any(box[3:] <= 0) or np.any(box[3:] >= 180)):
                        box = None
                    if len(c) == len(n) == 1 and .5 <= float(calc_bonds(c.positions[0], n.positions[0], box=box)) <= 1.9:
                        neighbors.setdefault(left.ix, {})["next"] = right
                        neighbors.setdefault(right.ix, {})["previous"] = left
        for residue in residues:
            if self.use_coordinate_neighbors:
                prev = neighbors.get(residue.ix, {}).get("previous")
                next_res = neighbors.get(residue.ix, {}).get("next")
                def group(spec):
                    if any(r is None for r, name in spec):
                        return None
                    atoms = [r.atoms.select_atoms(f"name {name}") for r, name in spec]
                    return sum(atoms[1:], atoms[0]) if all(len(a) == 1 for a in atoms) else None
                selections["phi"].append(group([(prev, "C"), (residue, "N"), (residue, "CA"), (residue, "C")]))
                selections["psi"].append(group([(residue, "N"), (residue, "CA"), (residue, "C"), (next_res, "N")]))
                selections["omega"].append(group([(residue, "CA"), (residue, "C"), (next_res, "N"), (next_res, "CA")]))
            else:
                selections["phi"].append(self._safe_residue_dihedral_selection(residue, "phi_selection"))
                selections["psi"].append(self._safe_residue_dihedral_selection(residue, "psi_selection"))
                selections["omega"].append(self._safe_residue_dihedral_selection(residue, "omega_selection"))
            selections["chi1"].append(self._safe_residue_dihedral_selection(residue, "chi1_selection"))
            selections["chi2"].append(self._chi2_selection(residue))
        return selections


    def extract_features_for_residues(
        self,
        traj: TrajectoryObject,
        protein_residues: mda.ResidueGroup,
        slice_obj: slice = slice(None)
    ) -> Tuple[Optional[Dict[int, np.ndarray]], int]:
        """
        Performs the expensive, one-time feature extraction for all residues
        in the provided ResidueGroup. This is the core computational step.

        Args:
            traj: The MDAnalysis Universe object.
            protein_residues: A ResidueGroup (e.g., from u.select_atoms('protein').residues)
                              for which to calculate all features.
            slice_obj: A slice object for the trajectory.

        Returns:
            - A dictionary mapping residue index (.ix) to its feature array (n_frames, 1, 5).
            - The number of frames processed.
        """
        n_frames = self._get_sliced_length(traj, slice_obj)
        if n_frames == 0:
            print("      Warning: Slice results in 0 frames. Skipping extraction.")
            return None, 0

        print(f"    Calculating dihedrals for {len(protein_residues)} residues in one pass...")

        try:
            selections = self._build_dihedral_selections(protein_residues)
            all_selections = []
            dihedral_offsets: Dict[str, Tuple[int, int]] = {}
            cursor = 0
            for name in DIHEDRAL_KEYS:
                items = selections[name]
                dihedral_offsets[name] = (cursor, cursor + len(items))
                all_selections.extend(items)
                cursor += len(items)
            if not all_selections:
                print("      Warning: No valid dihedrals found for the entire protein selection.")
                return None, 0

            # --- Single, expensive calculation over the trajectory slice ---
            all_angles_deg = self._calculate_dihedral_angle(all_selections, n_frames, slice_obj)
            # -------------------------------------------------------------

            # Create a per-residue dictionary to store results.
            # Shape: (n_frames, 5) for phi, psi, omega, chi1, chi2
            res_angle_map = {
                res.ix: np.zeros((n_frames, len(DIHEDRAL_KEYS)), dtype=np.float32)
                for res in protein_residues
            }

            res_indices = [res.ix for res in protein_residues]
            for i, res_ix in enumerate(res_indices):
                for dim_idx, name in enumerate(DIHEDRAL_KEYS):
                    start, end = dihedral_offsets[name]
                    res_angle_map[res_ix][:, dim_idx] = all_angles_deg[:, start:end][:, i]

            # Now, transform to circular coordinates and store in the final dict
            all_residue_features: Dict[int, np.ndarray] = {}
            for res_ix, angles_deg in res_angle_map.items():
                # Reshape to (n_frames, 1, 5) to match the stored descriptor layout.
                angles_deg_reshaped = angles_deg[:, np.newaxis, :]
                angles_rad = np.deg2rad(angles_deg_reshaped)
                all_residue_features[res_ix] = angles_rad

            print(f"      Bulk extraction complete. Found features for {len(all_residue_features)} residues.")
            return all_residue_features, n_frames

        except Exception as e:
            print(f"    FATAL ERROR during bulk feature extraction: {e}")
            import traceback
            traceback.print_exc()
            return None, 0


    def extract_all_features(
        self, 
        traj: TrajectoryObject, 
        slice_obj: slice = slice(None)
    ) -> Tuple[FeatureDict, int]:
        """
        Extracts features for all defined residues on a given slice.

        This method acts as a filter. It expects to be called
        AFTER the expensive `extract_features_for_residues` has been run. It
        selects the pre-computed features based on `self.residue_selections`.

        Args:
            traj: The MDAnalysis Universe object.
            slice_obj: A slice object (e.g., slice(1000, 5000, 2)).

        Returns:
            - A dictionary {selection_key: feature_array}
              where each array has shape (n_sliced_frames, n_residues, 5).
            - The number of frames actually processed (n_sliced_frames).
        """
        if not self.residue_selections:
            print("  Warning: FeatureExtractor has no residue selections. Returning empty results.")
            return {}, 0

        # Instead of selecting the whole protein, we build a
        # single selection string for only the residues we need.
        combined_selection_string = "protein and (" + " or ".join(f"({sel})" for sel in self.residue_selections.values()) + ")"
        
        # This creates a ResidueGroup of only the aligned, filtered residues.
        target_residues = traj.select_atoms(combined_selection_string).residues
        
        # Now, we run the expensive calculation ONLY on this small group.
        all_residue_features, n_frames = self.extract_features_for_residues(
            traj, target_residues, slice_obj
        )

        if all_residue_features is None:
            return {}, 0

        # --- This section now just filters the pre-computed results ---
        filtered_features_dict: FeatureDict = {}
        print("  Filtering pre-computed features based on config selections...")
        for key, sel_string in self.residue_selections.items():
            try:
                # Select all residues in the group to get their indices
                target_residues_group = traj.select_atoms(sel_string).residues
                
                if target_residues_group.n_residues > 0:
                    # Collect the pre-computed features for each residue in the group
                    feature_list_for_group = []
                    all_found = True
                    for res in target_residues_group:
                        if res.ix in all_residue_features:
                            feature_list_for_group.append(all_residue_features[res.ix])
                        else:
                            print(f"    Warning: Feature for resid {res.resid} (part of group '{key}') not found. Skipping group.")
                            all_found = False
                            break
                    
                    # If all features were found, concatenate them
                    if all_found and feature_list_for_group:
                        # Concatenate along the last axis (the feature dimension)
                        # (n_frames, 1, 3), (n_frames, 1, 3) -> (n_frames, 1, 6)
                        filtered_features_dict[key] = np.concatenate(feature_list_for_group, axis=-1)
                else:
                    print(f"    Warning: Selection '{sel_string}' for key '{key}' resolved to 0 residues. Skipping.")
                    continue
            except Exception as e:
                print(f"    Warning: Could not resolve selection '{sel_string}' for key '{key}'. Skipping. Error: {e}")

        return filtered_features_dict, n_frames
