from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Tuple

import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals

from phase.features.extraction import DIHEDRAL_KEYS, FeatureExtractor
from phase.services.project_store import ProjectStore


_CHI1_ATOM = {
    "ARG": "CG",
    "ASN": "CG",
    "ASP": "CG",
    "CYS": "SG",
    "GLN": "CG",
    "GLU": "CG",
    "HIS": "CG",
    "ILE": "CG1",
    "LEU": "CG",
    "LYS": "CG",
    "MET": "CG",
    "PHE": "CG",
    "PRO": "CG",
    "SER": "OG",
    "THR": "OG1",
    "TRP": "CG",
    "TYR": "CG",
    "VAL": "CG1",
}


def _extract_resid_from_key(key: str) -> int | None:
    match = re.search(r"(\d+)$", key)
    if match:
        return int(match.group(1))
    return None


def _build_residue_index_map(residue_keys: List[str]) -> Dict[int, int]:
    mapping: Dict[int, int] = {}
    for idx, key in enumerate(residue_keys):
        resid = _extract_resid_from_key(str(key))
        if resid is None:
            continue
        mapping[resid] = idx
    return mapping


def _sorted_residue_keys(residue_keys: List[str]) -> Tuple[List[str], np.ndarray, np.ndarray]:
    indexed: List[Tuple[int, int, str]] = []
    for idx, key in enumerate(residue_keys):
        resid = _extract_resid_from_key(str(key))
        resid_val = resid if resid is not None else 1_000_000 + idx
        indexed.append((resid_val, idx, str(key)))
    indexed.sort(key=lambda item: item[0])
    ordered = [key for _, _, key in indexed]
    order = np.array([idx for _, idx, _ in indexed], dtype=int)
    resids = np.array([_extract_resid_from_key(k) or -1 for k in ordered], dtype=int)
    return ordered, order, resids


def _resolve_selection(state_selection: str | None) -> str:
    if state_selection and state_selection.strip():
        return f"protein and ({state_selection.strip()})"
    return "protein"


def _normalize_label(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def _candidate_labels(raw: str) -> List[str]:
    raw = raw.strip()
    candidates = [raw]
    if "-" in raw:
        candidates.append(raw.split("-")[-1].strip())
    if "+" in raw:
        candidates.append(raw.split("+")[-1].strip())
    return candidates


def _map_frame_state_ids(frame_state_ids: np.ndarray, system_states: Dict[str, object]) -> np.ndarray:
    if not system_states:
        return frame_state_ids
    normalized_lookup = {
        sid: (_normalize_label(getattr(state, "name", sid) or sid), _normalize_label(sid))
        for sid, state in system_states.items()
    }
    mapping: Dict[str, str] = {}
    for raw_id in sorted({str(v) for v in frame_state_ids.tolist()}):
        if raw_id in system_states:
            mapping[raw_id] = raw_id
            continue
        label_candidates = _candidate_labels(raw_id)
        label_candidates = [
            c.replace("_descriptors", "").replace("descriptors", "").strip()
            for c in label_candidates
            if c
        ]
        normalized_candidates = {_normalize_label(c) for c in label_candidates if c}
        matched = None
        for sid, (norm_name, norm_id) in normalized_lookup.items():
            if norm_name in normalized_candidates or norm_id in normalized_candidates:
                matched = sid
                break
        if matched is not None:
            mapping[raw_id] = matched
    if not mapping:
        return frame_state_ids
    mapped = np.array([mapping.get(str(v), str(v)) for v in frame_state_ids], dtype=str)
    return mapped


def _build_dihedral_indices(u: mda.Universe, residue_index_map: Dict[int, int]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    phi_list: List[Tuple[int, int]] = []
    psi_list: List[Tuple[int, int]] = []
    chi_list: List[Tuple[int, int]] = []

    phi_atoms: List[Tuple[int, int, int, int]] = []
    psi_atoms: List[Tuple[int, int, int, int]] = []
    chi_atoms: List[Tuple[int, int, int, int]] = []

    residues = list(u.residues)
    for idx, res in enumerate(residues):
        resid = int(res.resid)
        if resid not in residue_index_map:
            continue
        r_index = residue_index_map[resid]

        prev_res = residues[idx - 1] if idx > 0 else None
        next_res = residues[idx + 1] if idx < len(residues) - 1 else None

        n_atom = res.atoms.select_atoms("name N")
        ca_atom = res.atoms.select_atoms("name CA")
        c_atom = res.atoms.select_atoms("name C")

        if prev_res is not None:
            prev_c = prev_res.atoms.select_atoms("name C")
            if prev_c.n_atoms and n_atom.n_atoms and ca_atom.n_atoms and c_atom.n_atoms:
                phi_atoms.append((prev_c[0].index, n_atom[0].index, ca_atom[0].index, c_atom[0].index))
                phi_list.append((r_index, resid))

        if next_res is not None:
            next_n = next_res.atoms.select_atoms("name N")
            if n_atom.n_atoms and ca_atom.n_atoms and c_atom.n_atoms and next_n.n_atoms:
                psi_atoms.append((n_atom[0].index, ca_atom[0].index, c_atom[0].index, next_n[0].index))
                psi_list.append((r_index, resid))

        chi_atom_name = _CHI1_ATOM.get(res.resname)
        if chi_atom_name:
            cb_atom = res.atoms.select_atoms("name CB")
            chi_atom = res.atoms.select_atoms(f"name {chi_atom_name}")
            if n_atom.n_atoms and ca_atom.n_atoms and cb_atom.n_atoms and chi_atom.n_atoms:
                chi_atoms.append((n_atom[0].index, ca_atom[0].index, cb_atom[0].index, chi_atom[0].index))
                chi_list.append((r_index, resid))

    def _to_array(items: List[Tuple[int, int, int, int]]) -> np.ndarray:
        if not items:
            return np.zeros((0, 4), dtype=int)
        return np.array(items, dtype=int)

    phi_idx = _to_array(phi_atoms)
    psi_idx = _to_array(psi_atoms)
    chi_idx = _to_array(chi_atoms)
    return phi_idx, psi_idx, chi_idx, np.array(phi_list), np.array(psi_list), np.array(chi_list)


def _compute_dihedrals_for_frame(positions: np.ndarray, idxs: np.ndarray) -> np.ndarray:
    if idxs.size == 0:
        return np.array([], dtype=float)
    a = positions[idxs[:, 0]]
    b = positions[idxs[:, 1]]
    c = positions[idxs[:, 2]]
    d = positions[idxs[:, 3]]
    return calc_dihedrals(a, b, c, d)


def _build_named_dihedral_indices(
    residues: mda.core.groups.ResidueGroup,
    residue_index_map: Dict[int, int],
    *, by_resindex: bool = False,
) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    extractor = FeatureExtractor(use_coordinate_neighbors=by_resindex)
    selections = extractor._build_dihedral_selections(residues)
    out: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for name in DIHEDRAL_KEYS:
        idx_rows: List[Tuple[int, int, int, int]] = []
        residue_rows: List[int] = []
        for residue, atom_group in zip(residues, selections.get(name) or []):
            if atom_group is None:
                continue
            try:
                atom_count = int(atom_group.n_atoms)
            except Exception:
                continue
            if atom_count != 4:
                continue
            residue_col = residue_index_map.get(int(residue.ix if by_resindex else residue.resid))
            if residue_col is None:
                continue
            idx_rows.append(tuple(int(atom.index) for atom in atom_group))
            residue_rows.append(int(residue_col))
        out[name] = (
            np.asarray(idx_rows, dtype=int) if idx_rows else np.zeros((0, 4), dtype=int),
            np.asarray(residue_rows, dtype=int) if residue_rows else np.zeros((0,), dtype=int),
        )
    return out


def _build_dense_dihedral_atom_indices(
    dihedral_indices: Dict[str, Tuple[np.ndarray, np.ndarray]],
    *,
    n_residues: int,
    atom_index_to_local: Dict[int, int],
) -> tuple[np.ndarray, np.ndarray]:
    dense = np.full((n_residues, len(DIHEDRAL_KEYS), 4), -1, dtype=np.int32)
    mask = np.zeros((n_residues, len(DIHEDRAL_KEYS)), dtype=bool)
    for dim_idx, name in enumerate(DIHEDRAL_KEYS):
        idxs, residue_rows = dihedral_indices.get(name, (np.zeros((0, 4), dtype=int), np.zeros((0,), dtype=int)))
        if idxs.size == 0 or residue_rows.size == 0:
            continue
        for row_atoms, residue_row in zip(np.asarray(idxs, dtype=int), np.asarray(residue_rows, dtype=int)):
            local = [atom_index_to_local.get(int(atom_idx), -1) for atom_idx in np.asarray(row_atoms, dtype=int).tolist()]
            dense[int(residue_row), dim_idx, :] = np.asarray(local, dtype=np.int32)
            mask[int(residue_row), dim_idx] = all(v >= 0 for v in local)
    return dense, mask


def _load_cluster_label_metadata(cluster_path: Path) -> tuple[np.ndarray, np.ndarray]:
    with np.load(cluster_path, allow_pickle=True) as cluster_npz:
        residue_keys = np.asarray(cluster_npz["residue_keys"]).astype(str)
        if "cluster_counts" in cluster_npz:
            cluster_counts = np.asarray(cluster_npz["cluster_counts"], dtype=np.int32)
        else:
            cluster_counts = np.asarray(cluster_npz["merged__cluster_counts"], dtype=np.int32)
    return residue_keys, cluster_counts


def build_sample_backmapping_dataset(
    project_id: str,
    system_id: str,
    cluster_id: str,
    sample_id: str,
    trajectory_path: Path,
    output_path: Path,
    progress_callback: Callable[[int, int], None] | None = None,
) -> Dict[str, Any]:
    store = ProjectStore()
    system = store.get_system(project_id, system_id)
    sample = store.get_sample_entry(project_id, system_id, cluster_id, sample_id)
    if str(sample.get("type") or "") != "md_eval":
        raise ValueError("Backmapping dataset is only supported for md_eval samples.")
    state_id = str(sample.get("state_id") or "").strip()
    if not state_id:
        raise ValueError("Backmapping dataset currently requires an md_eval sample tied to a single state.")
    state = system.states.get(state_id)
    if not state:
        raise ValueError(f"State '{state_id}' referenced by sample '{sample_id}' is missing.")
    if not state.pdb_file:
        raise ValueError(f"State '{state_id}' is missing its PDB file.")

    cluster_dir = store.ensure_cluster_directories(project_id, system_id, cluster_id)["cluster_dir"]
    cluster_path = cluster_dir / "cluster.npz"
    if not cluster_path.exists():
        raise FileNotFoundError(f"Cluster NPZ not found for cluster '{cluster_id}'.")
    residue_keys, cluster_counts = _load_cluster_label_metadata(cluster_path)

    raw_sample_npz_path = str(sample.get("path") or (sample.get("paths") or {}).get("summary_npz") or "").strip()
    if not raw_sample_npz_path:
        raise FileNotFoundError(f"Sample '{sample_id}' is missing its NPZ path.")
    sample_npz_path = Path(raw_sample_npz_path)
    if not sample_npz_path.is_absolute():
        sample_npz_path = store.resolve_path(project_id, system_id, str(sample_npz_path))
    if not sample_npz_path.exists():
        raise FileNotFoundError(f"Sample NPZ is missing on disk: {sample_npz_path}")

    with np.load(sample_npz_path, allow_pickle=True) as sample_npz:
        labels = np.asarray(sample_npz["labels"], dtype=np.int32)
        frame_indices = np.asarray(sample_npz["frame_indices"], dtype=np.int64)
        frame_state_ids = (
            np.asarray(sample_npz["frame_state_ids"]).astype(str)
            if "frame_state_ids" in sample_npz
            else np.full(labels.shape[0], state_id, dtype=str)
        )

    if labels.ndim != 2:
        raise ValueError("Sample labels must have shape (n_frames, n_residues).")
    if labels.shape[1] != residue_keys.shape[0]:
        raise ValueError(
            f"Sample residue count mismatch: expected {residue_keys.shape[0]}, got {labels.shape[1]}."
        )
    if frame_indices.shape[0] != labels.shape[0]:
        raise ValueError("Sample frame_indices length does not match label rows.")
    normalized_frame_state_ids = np.asarray([str(v) for v in frame_state_ids.tolist()], dtype=object)
    unique_frame_ids = {str(v) for v in normalized_frame_state_ids.tolist()}
    if unique_frame_ids == {state_id[:1]}:
        normalized_frame_state_ids = np.asarray([state_id] * int(labels.shape[0]), dtype=object)
        unique_frame_ids = {state_id}
    if unique_frame_ids != {state_id}:
        raise ValueError(
            f"Backmapping dataset only supports samples with frames from one state; found frame_state_ids={sorted(unique_frame_ids)!r} for state_id={state_id!r}."
        )
    frame_state_ids = np.asarray([state_id] * int(labels.shape[0]), dtype=f"U{max(1, len(state_id))}")

    pdb_path = store.resolve_path(project_id, system_id, str(state.pdb_file))
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file missing on disk: {pdb_path}")
    if not trajectory_path.exists():
        raise FileNotFoundError(f"Uploaded trajectory is missing on disk: {trajectory_path}")

    selection = _resolve_selection(state.residue_selection)
    universe = mda.Universe(str(pdb_path), str(trajectory_path))
    atoms = universe.select_atoms(selection)
    if atoms.n_atoms == 0:
        raise ValueError(f"Selection '{selection}' yielded no atoms for state '{state_id}'.")

    residue_index_map = _build_residue_index_map([str(key) for key in residue_keys.tolist()])
    reference_aligned = state.alignment_status in {"reference", "aligned", "partial"}
    if reference_aligned:
        from phase.workflows.reference_states import canonical_mapping
        correspondence = canonical_mapping(store, system, state)
        residue_index_map = {}
        for column, key in enumerate(residue_keys.tolist()):
            selected = universe.select_atoms(correspondence[str(key)]).residues
            if len(selected) != 1:
                raise ValueError(f"Canonical backmapping residue {key} does not resolve uniquely.")
            residue_index_map[int(selected[0].ix)] = column
    atom_residue_index = np.full(atoms.n_atoms, -1, dtype=np.int32)
    atom_resids = np.asarray([int(atom.resid) for atom in atoms], dtype=np.int32)
    atom_names = np.asarray([str(atom.name) for atom in atoms], dtype=str)
    atom_index_to_local = {int(atom.index): int(idx_atom) for idx_atom, atom in enumerate(atoms)}
    for idx_atom, atom in enumerate(atoms):
        mapped = residue_index_map.get(int(atom.resindex if reference_aligned else atom.resid))
        if mapped is not None:
            atom_residue_index[idx_atom] = int(mapped)

    n_rows, n_residues = labels.shape
    coordinates = np.zeros((n_rows, atoms.n_atoms, 3), dtype=np.float32)
    dihedrals = np.full((n_rows, n_residues, len(DIHEDRAL_KEYS)), np.nan, dtype=np.float32)

    dihedral_indices = _build_named_dihedral_indices(atoms.residues, residue_index_map, by_resindex=reference_aligned)
    dihedral_atom_indices, dihedral_mask = _build_dense_dihedral_atom_indices(
        dihedral_indices,
        n_residues=n_residues,
        atom_index_to_local=atom_index_to_local,
    )

    frame_to_rows: Dict[int, List[int]] = {}
    for row_idx, frame_idx in enumerate(frame_indices.tolist()):
        frame_to_rows.setdefault(int(frame_idx), []).append(int(row_idx))

    total_rows = int(n_rows)
    processed_rows = 0
    seen_rows: set[int] = set()
    for ts in universe.trajectory:
        frame_idx = int(ts.frame)
        rows = frame_to_rows.get(frame_idx)
        if not rows:
            continue
        positions = universe.atoms.positions.astype(np.float32)
        selected_positions = atoms.positions.astype(np.float32)
        angle_values = {
            name: _compute_dihedrals_for_frame(positions, idxs)
            for name, (idxs, _) in dihedral_indices.items()
        }
        for row_idx in rows:
            coordinates[row_idx] = selected_positions
            for dim_idx, name in enumerate(DIHEDRAL_KEYS):
                _, residue_rows = dihedral_indices[name]
                values = angle_values[name]
                if residue_rows.size and values.size:
                    dihedrals[row_idx, residue_rows, dim_idx] = values.astype(np.float32)
            seen_rows.add(int(row_idx))
        processed_rows += len(rows)
        if progress_callback:
            progress_callback(processed_rows, total_rows)

    if len(seen_rows) != total_rows:
        missing = sorted(set(range(total_rows)) - seen_rows)
        raise ValueError(
            f"Trajectory does not provide all requested sample frames. Missing frame rows: {missing[:10]}"
            + ("..." if len(missing) > 10 else "")
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "trajectory": coordinates,
        "atom_resids": atom_resids,
        "atom_names": atom_names,
        "atom_residue_index": atom_residue_index,
        "residue_keys": residue_keys,
        "residue_cluster_ids": labels,
        "residue_cluster_counts": cluster_counts,
        "frame_indices": frame_indices,
        "frame_state_ids": frame_state_ids,
        "state_id": np.asarray(state_id, dtype=str),
        "sample_id": np.asarray(sample_id, dtype=str),
        "dihedrals": dihedrals,
        "dihedral_keys": np.asarray(DIHEDRAL_KEYS, dtype=str),
        "dihedral_atom_indices": dihedral_atom_indices,
        "dihedral_mask": dihedral_mask,
    }
    np.savez_compressed(output_path, **payload)
    return {
        "path": str(output_path),
        "n_frames": int(coordinates.shape[0]),
        "n_atoms": int(coordinates.shape[1]),
        "n_residues": int(labels.shape[1]),
        "dihedral_keys": list(DIHEDRAL_KEYS),
        "sample_id": str(sample_id),
        "state_id": str(state_id),
    }


def build_backmapping_npz(
    project_id: str,
    system_id: str,
    cluster_path: Path,
    output_path: Path,
    progress_callback: Callable[[int, int], None] | None = None,
    trajectory_overrides: Dict[str, Path] | None = None,
) -> Path:
    store = ProjectStore()
    system = store.get_system(project_id, system_id)

    with np.load(cluster_path, allow_pickle=True) as cluster_npz:
        residue_keys = [str(v) for v in cluster_npz["residue_keys"].tolist()]
        labels = np.asarray(cluster_npz["merged__labels_assigned"], dtype=np.int32)
        if "cluster_counts" in cluster_npz:
            cluster_counts = np.asarray(cluster_npz["cluster_counts"], dtype=np.int32)
        else:
            cluster_counts = np.asarray(cluster_npz["merged__cluster_counts"], dtype=np.int32)
        frame_state_ids = np.asarray(cluster_npz["merged__frame_state_ids"]).astype(str)
        frame_indices = np.asarray(cluster_npz["merged__frame_indices"], dtype=np.int64)

    n_frames, n_residues = labels.shape
    ordered_keys, order, ordered_resids = _sorted_residue_keys(residue_keys)
    residue_index_map = _build_residue_index_map(ordered_keys)

    frame_state_ids = _map_frame_state_ids(frame_state_ids, system.states or {})
    valid_state_ids = {str(sid) for sid in (system.states or {}).keys()}
    if valid_state_ids:
        keep_mask = np.array([sid in valid_state_ids for sid in frame_state_ids], dtype=bool)
        if not np.any(keep_mask):
            raise ValueError("No frames match available system states for backmapping.")
        labels = labels[keep_mask]
        frame_state_ids = frame_state_ids[keep_mask]
        frame_indices = frame_indices[keep_mask]
        n_frames = labels.shape[0]

    if order.size and order.shape[0] == labels.shape[1]:
        labels = labels[:, order]
        cluster_counts = cluster_counts[order]

    output_path.parent.mkdir(parents=True, exist_ok=True)

    xyz = None
    atom_residue_index = None
    phi = np.full((n_frames, n_residues), np.nan, dtype=np.float32)
    psi = np.full((n_frames, n_residues), np.nan, dtype=np.float32)
    chi1 = np.full((n_frames, n_residues), np.nan, dtype=np.float32)
    total_rows = n_frames
    processed_rows = 0

    state_to_rows: Dict[str, List[int]] = {}
    for idx, state_id in enumerate(frame_state_ids):
        state_to_rows.setdefault(state_id, []).append(idx)

    for state_id, rows in state_to_rows.items():
        state = system.states.get(state_id)
        if not state:
            continue
        override_path = None
        if trajectory_overrides:
            override_path = trajectory_overrides.get(str(state_id))
        if not state.trajectory_file and override_path is None:
            raise ValueError(f"State '{state_id}' is missing trajectory; please upload trajectory first.")
        if not state.pdb_file:
            raise ValueError(f"State '{state_id}' is missing PDB; please upload PDB first.")

        traj_path = override_path
        if traj_path is None:
            traj_path = store.resolve_path(project_id, system_id, state.trajectory_file)
        pdb_path = store.resolve_path(project_id, system_id, state.pdb_file)
        if not traj_path.exists():
            raise ValueError(f"Trajectory file missing on disk: {traj_path}")
        if not pdb_path.exists():
            raise ValueError(f"PDB file missing on disk: {pdb_path}")

        selection = _resolve_selection(state.residue_selection)
        u = mda.Universe(str(pdb_path), str(traj_path))
        atoms = u.select_atoms(selection)
        if atoms.n_atoms == 0:
            raise ValueError(f"Selection '{selection}' yielded no atoms for state '{state_id}'.")

        if xyz is None:
            xyz = np.zeros((n_frames, atoms.n_atoms, 3), dtype=np.float32)
            atom_residue_index = np.full(atoms.n_atoms, -1, dtype=np.int32)
            for idx_atom, atom in enumerate(atoms):
                resid = int(atom.resid)
                if resid in residue_index_map:
                    atom_residue_index[idx_atom] = residue_index_map[resid]
        else:
            if atoms.n_atoms != xyz.shape[1]:
                raise ValueError("Atom selection size mismatch across states; cannot combine trajectories.")

        phi_idx, psi_idx, chi_idx, phi_map, psi_map, chi_map = _build_dihedral_indices(u, residue_index_map)

        needed_frame_indices = frame_indices[rows]
        frame_to_rows: Dict[int, List[int]] = {}
        for row_idx, frame_idx in zip(rows, needed_frame_indices):
            frame_to_rows.setdefault(int(frame_idx), []).append(row_idx)

        for ts in u.trajectory:
            frame_idx = int(ts.frame)
            if frame_idx not in frame_to_rows:
                continue
            positions = u.atoms.positions.astype(np.float32)
            sel_positions = atoms.positions.astype(np.float32)
            phi_vals = _compute_dihedrals_for_frame(positions, phi_idx)
            psi_vals = _compute_dihedrals_for_frame(positions, psi_idx)
            chi_vals = _compute_dihedrals_for_frame(positions, chi_idx)

            for row_idx in frame_to_rows[frame_idx]:
                xyz[row_idx] = sel_positions
                for val_idx, (res_idx, _) in enumerate(phi_map):
                    phi[row_idx, res_idx] = float(phi_vals[val_idx])
                for val_idx, (res_idx, _) in enumerate(psi_map):
                    psi[row_idx, res_idx] = float(psi_vals[val_idx])
                for val_idx, (res_idx, _) in enumerate(chi_map):
                    chi1[row_idx, res_idx] = float(chi_vals[val_idx])
            processed_rows += len(frame_to_rows[frame_idx])
            if progress_callback:
                progress_callback(processed_rows, total_rows)

    if xyz is None or atom_residue_index is None:
        raise ValueError("No trajectory data was collected.")

    payload = {
        "xyz": xyz,
        "atom_residue_index": atom_residue_index,
        "residue_keys": np.array(ordered_keys, dtype=str),
        "residue_resids": ordered_resids,
        "residue_cluster_labels": labels,
        "residue_cluster_counts": cluster_counts,
        "phi": phi,
        "psi": psi,
        "chi1": chi1,
        "frame_state_ids": frame_state_ids,
        "frame_indices": frame_indices,
    }

    np.savez_compressed(output_path, **payload)
    return output_path
