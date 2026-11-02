"""Optional C-alpha distance data aligned to a cluster's unfiltered merged rows.

Distances remain independent of labels, bin thresholds, and Potts parameters.
The only authoritative distance metadata lives inside distance_cache.npz.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import tempfile
from typing import Callable, Sequence

import numpy as np

from phase.common.slice_utils import descriptor_frame_to_source_frame
from phase.services.project_store import DescriptorState, ProjectStore


DISTANCE_CACHE_FILENAME = "distance_cache.npz"
SCHEMA_VERSION = 1


@dataclass(frozen=True)
class ClusterFrames:
    residue_keys: np.ndarray
    labels: np.ndarray
    frame_state_ids: np.ndarray
    frame_metastable_ids: np.ndarray
    frame_indices: np.ndarray
    contact_edge_index: np.ndarray
    metadata: dict


@dataclass(frozen=True)
class DistanceCache:
    residue_keys: np.ndarray
    edge_index: np.ndarray
    frame_state_ids: np.ndarray
    frame_metastable_ids: np.ndarray
    frame_indices: np.ndarray
    distances: np.ndarray | None
    metadata: dict

    @property
    def available(self) -> bool:
        return self.distances is not None


def _metadata(data) -> dict:
    return json.loads(str(data["metadata_json"].item())) if "metadata_json" in data else {}


def _edges(array, n_residues: int) -> np.ndarray:
    edges = np.asarray(array)
    if not np.issubdtype(edges.dtype, np.integer):
        raise ValueError("edge_index must contain integers.")
    if edges.ndim != 2 or edges.shape[0] != 2:
        raise ValueError("edge_index must have shape (2, E).")
    if edges.size and (np.any(edges < 0) or np.any(edges >= n_residues) or np.any(edges[0] >= edges[1])):
        raise ValueError("Distance edges must be in bounds, canonical (r < s), and unique.")
    if len(set(map(tuple, edges.T))) != edges.shape[1]:
        raise ValueError("Distance edges must be unique.")
    return edges


def read_cluster_frames(cluster_path: Path | str) -> ClusterFrames:
    """Read all merged rows, including halo/invalid labels and fitting holdouts."""
    with np.load(cluster_path, allow_pickle=False) as data:
        meta = _metadata(data)
        keys = (meta.get("merged") or {}).get("npz_keys") or {}
        def key(name, default):
            return str(keys.get(name) or default)
        assigned_key = key("labels_assigned", "merged__labels_assigned")
        label_key = assigned_key if assigned_key in data else key("labels", "merged__labels")
        labels = np.asarray(data[label_key], dtype=np.int32)
        residues = np.asarray(data["residue_keys"], dtype=str)
        states = np.asarray(data[key("frame_state_ids", "merged__frame_state_ids")], dtype=str)
        frames = np.asarray(data[key("frame_indices", "merged__frame_indices")])
        meta_key = key("frame_metastable_ids", "merged__frame_metastable_ids")
        metastable = np.asarray(data[meta_key], dtype=str) if meta_key in data else np.full(len(states), "", dtype="U1")
        if labels.ndim != 2 or labels.shape != (len(states), len(residues)):
            raise ValueError("Merged labels, residue keys, and frame identities are misaligned.")
        if residues.ndim != 1 or len(set(residues)) != len(residues):
            raise ValueError("Cluster residue keys must be a unique 1D array.")
        if states.ndim != 1 or metastable.shape != states.shape or frames.shape != states.shape:
            raise ValueError("Merged frame identity arrays must be aligned 1D arrays.")
        if not np.issubdtype(frames.dtype, np.integer) or np.any(frames < 0):
            raise ValueError("Merged frame_indices must contain nonnegative integers.")
        contact = _edges(data["contact_edge_index"], len(residues)) if "contact_edge_index" in data else np.empty((2, 0), dtype=np.int64)
        return ClusterFrames(residues, labels, states, metastable, frames.astype(np.int64), contact, meta)


def validate_distance_alignment(cluster: ClusterFrames, cache: DistanceCache) -> None:
    """Equal row counts are insufficient: compare every frame and residue identity."""
    if not np.issubdtype(cache.frame_indices.dtype, np.integer):
        raise ValueError("Distance cache frame_indices must contain integers.")
    for name in ("residue_keys", "frame_state_ids", "frame_metastable_ids", "frame_indices"):
        if not np.array_equal(getattr(cluster, name), getattr(cache, name)):
            raise ValueError(f"Distance cache alignment mismatch: {name}. Rebuild the distance cache.")
    _edges(cache.edge_index, len(cluster.residue_keys))
    if cache.available and cache.distances.shape != (len(cluster.frame_indices), cache.edge_index.shape[1]):
        raise ValueError("Distance matrix shape does not match frame/edge identities.")


def load_distance_cache(cluster_path: Path | str, cache_path: Path | str | None = None) -> DistanceCache:
    """Missing cache means one unconditional bin; malformed/stale caches raise."""
    cluster_path = Path(cluster_path)
    cluster = read_cluster_frames(cluster_path)
    path = Path(cache_path) if cache_path is not None else cluster_path.with_name(DISTANCE_CACHE_FILENAME)
    if not path.exists():
        return DistanceCache(
            cluster.residue_keys, cluster.contact_edge_index, cluster.frame_state_ids,
            cluster.frame_metastable_ids, cluster.frame_indices, None,
            {"schema_version": SCHEMA_VERSION, "distance_available": False,
             "bin_labels": ["All distances"], "fallback": "single_unconditional_bin"},
        )
    with np.load(path, allow_pickle=False) as data:
        meta = _metadata(data)
        if meta.get("schema_version") != SCHEMA_VERSION or meta.get("distance_type") != "ca_ca":
            raise ValueError("Unsupported distance cache schema or distance_type.")
        if meta.get("distance_unit") != "angstrom":
            raise ValueError("Distance cache must use angstrom units.")
        cache = DistanceCache(
            np.asarray(data["residue_keys"], dtype=str), np.asarray(data["edge_index"]),
            np.asarray(data["frame_state_ids"], dtype=str), np.asarray(data["frame_metastable_ids"], dtype=str),
            np.asarray(data["frame_indices"]), np.asarray(data["distances"], dtype=np.float32), meta,
        )
    validate_distance_alignment(cluster, cache)
    return cache


def distance_bin_indices(distances: np.ndarray | None, thresholds: Sequence[float] = (), *, frame_count: int | None = None) -> np.ndarray:
    """[min, max) bins; exact thresholds belong to the bin on their right.

    Without geometry, return a single bin regardless of requested thresholds.
    Missing values in an existing cache are excluded (-1), never invented.
    """
    if distances is None:
        if frame_count is None or frame_count < 0:
            raise ValueError("frame_count is required for the single-bin fallback.")
        return np.zeros(frame_count, dtype=np.int32)
    boundaries = np.asarray(thresholds, dtype=float)
    if boundaries.ndim != 1 or not np.all(np.isfinite(boundaries)) or np.any(boundaries <= 0) or np.any(np.diff(boundaries) <= 0):
        raise ValueError("Distance thresholds must be finite, positive, and strictly increasing.")
    values = np.asarray(distances)
    result = np.searchsorted(boundaries, values, side="right").astype(np.int32)
    result[~np.isfinite(values)] = -1
    return result


def validate_distance_preset(cache: DistanceCache, thresholds: Sequence[float], distance_type: str = "ca_ca") -> None:
    if distance_type != "ca_ca":
        raise ValueError("Only ca_ca distance presets are supported.")
    distance_bin_indices(np.empty(0), thresholds)
    if cache.available and len(thresholds) and max(thresholds) > float(cache.metadata["max_distance_cutoff"]):
        raise ValueError("Preset thresholds exceed the candidate-edge cutoff; rebuild/extend the cache.")


def _frame_mask(cache: DistanceCache, state_ids=None, metastable_ids=None) -> np.ndarray:
    mask = np.ones(len(cache.frame_indices), dtype=bool)
    for requested, available in ((state_ids, cache.frame_state_ids), (metastable_ids, cache.frame_metastable_ids)):
        if requested is not None:
            requested = list(requested)
            unknown = set(requested) - set(available)
            if unknown:
                raise ValueError(f"Unknown frame filter IDs: {sorted(unknown)}")
            mask &= np.isin(available, requested)
    return mask


def pair_distance_histogram(
    cluster_path: Path | str, r: int, s: int, *, cache_path=None,
    state_ids=None, metastable_ids=None, thresholds: Sequence[float] = (), bins: int = 60,
) -> dict:
    """Return compact pair data and occupancies, never a full T x E payload."""
    cache = load_distance_cache(cluster_path, cache_path)
    r, s = sorted((int(r), int(s)))
    if r < 0 or s >= len(cache.residue_keys) or r == s:
        raise ValueError("Pair indices must identify two distinct cluster residues.")
    if bins < 1 or bins > 10000:
        raise ValueError("Histogram bins must be between 1 and 10000.")
    validate_distance_preset(cache, thresholds)
    mask = _frame_mask(cache, state_ids, metastable_ids)
    n_frames = int(mask.sum())
    values = None
    if cache.available:
        found = np.flatnonzero((cache.edge_index[0] == r) & (cache.edge_index[1] == s))
        if not len(found):
            raise ValueError(f"Pair ({r}, {s}) is outside the stored candidate-edge set.")
        values = cache.distances[mask, int(found[0])]
        if not np.all(np.isfinite(values)):
            raise ValueError("Selected pair contains missing distance values; rebuild the cache.")
    indices = distance_bin_indices(values, thresholds, frame_count=n_frames)
    labels = ["All distances"]
    if cache.available and len(thresholds):
        boundaries = list(thresholds)
        labels = [f"< {boundaries[0]:g} A"] + [
            f"{a:g} <= d < {b:g} A" for a, b in zip(boundaries, boundaries[1:])
        ] + [f">= {boundaries[-1]:g} A"]
    counts = np.bincount(indices, minlength=len(labels))
    histogram, edges = np.histogram(values, bins=bins) if values is not None and len(values) else (np.empty(0, dtype=int), np.empty(0))
    summary = None
    if values is not None and len(values):
        q10, median, q90 = np.quantile(values, [0.1, 0.5, 0.9])
        summary = {"min": float(values.min()), "max": float(values.max()), "median": float(median), "q10": float(q10), "q90": float(q90)}
    return {
        "pair": [r, s], "residue_keys": cache.residue_keys[[r, s]].tolist(),
        "distance_available": cache.available, "frame_count": n_frames,
        "state_count": len(set(cache.frame_state_ids[mask])), "distance_unit": "angstrom",
        "bin_labels": labels, "bin_counts": counts.tolist(),
        "bin_fractions": (counts / max(1, n_frames)).tolist(),
        "histogram": {"counts": histogram.tolist(), "edges": edges.tolist()}, "summary": summary,
    }


def _ca_atoms(universe, state: DescriptorState, residue_keys: np.ndarray, cluster_mapping: dict):
    selection = str(state.residue_selection or "protein")
    if selection == "system_selections":
        selection = "protein"
    base = universe.select_atoms(f"protein and ({selection})")
    indices = []
    for key in residue_keys:
        raw_selection = (state.residue_mapping or {}).get(str(key))
        if not raw_selection:
            match = re.fullmatch(r"res_(-?\d+)(?:_[A-Za-z]+)?", str(key))
            if match:
                raw_selection = f"resid {int(match.group(1)) - int(state.resid_shift or 0)}"
            elif not state.residue_mapping and not state.resid_shift:
                raw_selection = cluster_mapping.get(str(key))
        if not raw_selection:
            raise ValueError(f"State '{state.state_id}': no single-residue mapping for {key}.")
        atoms = base.select_atoms(f"({raw_selection}) and name CA")
        if len(atoms) != 1:
            raise ValueError(f"State '{state.state_id}': {key} maps to {len(atoms)} CA atoms (expected exactly one).")
        indices.append(int(atoms.indices[0]))
    if len(set(indices)) != len(indices):
        raise ValueError(f"State '{state.state_id}': residue keys map to duplicate CA atoms.")
    return universe.atoms[indices]


def _valid_box(dimensions, periodic: bool):
    if not periodic or dimensions is None:
        return None
    box = np.asarray(dimensions, dtype=float)
    if box.shape != (6,) or not np.all(np.isfinite(box)) or np.any(box[:3] <= 0) or np.any(box[3:] <= 0) or np.any(box[3:] >= 180):
        return None
    return box


def build_distance_cache(
    project_id: str, system_id: str, cluster_id: str, *, store: ProjectStore | None = None,
    max_distance_cutoff: float = 12.0, periodic: bool = True,
    trajectory_overrides: dict[str, Path | str] | None = None,
    output_path: Path | str | None = None,
    progress_callback: Callable[[str, int, int], None] | None = None,
) -> dict:
    """Two passes over exact merged frames: candidate union, then all distances.

    Only the disk-backed T x E matrix is assembled. No T x N x N or full
    trajectory coordinate tensor is allocated. Atomic replacement preserves an
    existing cache if any source/mapping/frame validation fails.
    """
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import calc_bonds, self_capped_distance

    cutoff = float(max_distance_cutoff)
    if not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("max_distance_cutoff must be finite and positive.")
    store = store if store is not None else ProjectStore()
    system = store.get_system(project_id, system_id)
    entry = store.get_cluster_entry(project_id, system_id, cluster_id)
    cluster_path = store.resolve_path(project_id, system_id, str(entry["path"]))
    cluster = read_cluster_frames(cluster_path)
    if not len(cluster.frame_indices):
        raise ValueError("Cannot build distances for an empty cluster.")
    overrides = trajectory_overrides or {}
    if set(overrides) - set(cluster.frame_state_ids):
        raise ValueError("Trajectory overrides contain states absent from the merged cluster.")
    path = Path(output_path) if output_path is not None else cluster_path.with_name(DISTANCE_CACHE_FILENAME)
    if path.resolve() == cluster_path.resolve():
        raise ValueError("Distance cache cannot overwrite cluster.npz.")
    path.parent.mkdir(parents=True, exist_ok=True)
    candidates: set[tuple[int, int]] = set()
    sources = {}
    n_frames = len(cluster.frame_indices)

    def frames(stage):
        completed = 0
        for state_id in dict.fromkeys(cluster.frame_state_ids.tolist()):
            state = system.states.get(state_id)
            if state is None or not state.pdb_file:
                raise ValueError(f"Missing topology/state metadata for '{state_id}'.")
            top = store.resolve_path(project_id, system_id, state.pdb_file)
            traj_ref = overrides.get(state_id) or state.trajectory_file or state.source_traj or state.pdb_file
            traj = store.resolve_path(project_id, system_id, str(traj_ref))
            if not top.is_file() or not traj.is_file():
                raise FileNotFoundError(f"State '{state_id}': topology/trajectory unavailable: {top}, {traj}")
            universe = mda.Universe(str(top), str(traj))
            try:
                atoms = _ca_atoms(universe, state, cluster.residue_keys, cluster.metadata.get("residue_mapping") or {})
                slice_spec = state.slice_spec or (f"::{state.stride}" if int(state.stride or 1) > 1 else None)
                rows = np.flatnonzero(cluster.frame_state_ids == state_id)
                # Repeated source frames (e.g. overlapping macro/metastable
                # selections) are read once and written to every merged row.
                groups: dict[int, list[int]] = {}
                for row in rows:
                    source_frame = descriptor_frame_to_source_frame(int(cluster.frame_indices[row]), slice_spec, len(universe.trajectory))
                    groups.setdefault(source_frame, []).append(int(row))
                provenance = {
                    "topology": str(top), "trajectory": str(traj), "slice_spec": slice_spec,
                    "stride": int(state.stride or 1), "resid_shift": int(state.resid_shift or 0),
                    "ca_atom_indices": atoms.indices.tolist(), "source_frame_count": len(universe.trajectory),
                    "source_size_bytes": traj.stat().st_size, "source_mtime_ns": traj.stat().st_mtime_ns,
                    "topology_size_bytes": top.stat().st_size, "topology_mtime_ns": top.stat().st_mtime_ns,
                }
                if state_id in sources and sources[state_id] != provenance:
                    raise ValueError(f"State '{state_id}': trajectory/topology changed between distance passes.")
                sources[state_id] = provenance
                for source_frame in sorted(groups):
                    ts = universe.trajectory[source_frame]
                    positions = atoms.positions
                    if not np.all(np.isfinite(positions)):
                        raise ValueError(f"State '{state_id}', frame {source_frame}: non-finite CA coordinates.")
                    yield groups[source_frame], positions, _valid_box(ts.dimensions, periodic)
                    completed += len(groups[source_frame])
                    if progress_callback:
                        progress_callback(stage, completed, n_frames)
            finally:
                universe.trajectory.close()

    for _, positions, box in frames("Finding candidate pairs"):
        pairs = self_capped_distance(positions, max_cutoff=cutoff, box=box, return_distances=False)
        candidates.update((min(int(a), int(b)), max(int(a), int(b))) for a, b in pairs)
    edge_index = np.asarray(sorted(candidates), dtype=np.int64).reshape(-1, 2).T
    with tempfile.TemporaryDirectory(prefix=".distance-cache-", dir=path.parent) as work:
        work = Path(work)
        shape = (n_frames, len(candidates))
        distances = np.lib.format.open_memmap(work / "distances.npy", mode="w+", dtype=np.float32, shape=shape) if candidates else np.empty(shape, dtype=np.float32)
        for rows, positions, box in frames("Computing aligned distances"):
            distances[rows] = calc_bonds(positions[edge_index[0]], positions[edge_index[1]], box=box).astype(np.float32)
        contacts = set(map(tuple, cluster.contact_edge_index.T))
        is_contact = np.asarray([tuple(edge) in contacts for edge in edge_index.T], dtype=bool)
        meta = {
            "schema_version": SCHEMA_VERSION, "distance_type": "ca_ca", "distance_unit": "angstrom",
            "coordinate_source": "atomistic_trajectory", "distance_available": True,
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "project_id": project_id, "system_id": system_id, "cluster_id": cluster_id,
            "max_distance_cutoff": cutoff, "residue_count": len(cluster.residue_keys),
            "edge_count": len(candidates), "frame_count": n_frames,
            "selected_state_ids": list(dict.fromkeys(cluster.frame_state_ids.tolist())),
            "selected_metastable_ids": list(dict.fromkeys(cluster.frame_metastable_ids.tolist())),
            "frame_alignment": {"state_ids": "merged__frame_state_ids", "metastable_ids": "merged__frame_metastable_ids", "frame_indices": "merged__frame_indices"},
            "frame_index_semantics": "descriptor row; source frame resolved through state slice/stride",
            "candidate_rule": "union of pairs at or below cutoff over every merged frame",
            "distance_coverage": "all merged frames, including distances above candidate cutoff",
            "periodic_policy": "minimum_image_when_valid_box" if periodic else "nonperiodic",
            "contact_edges_count": cluster.contact_edge_index.shape[1], "stored_contact_edges_count": int(is_contact.sum()),
            "contact_cutoff": cluster.metadata.get("contact_cutoff", entry.get("contact_cutoff")),
            "distance_candidate_edges_count": len(candidates), "sources": sources,
            "subsampling": None,
        }
        temporary = work / DISTANCE_CACHE_FILENAME
        np.savez_compressed(
            temporary, residue_keys=cluster.residue_keys, edge_index=edge_index,
            edge_residue_keys=cluster.residue_keys[edge_index.T], is_contact_edge=is_contact,
            frame_state_ids=cluster.frame_state_ids, frame_metastable_ids=cluster.frame_metastable_ids,
            frame_indices=cluster.frame_indices, distances=distances,
            metadata_json=np.asarray(json.dumps(meta)),
        )
        # Re-read the cluster to detect frame/residue changes during extraction.
        load_distance_cache(cluster_path, temporary)
        temporary.replace(path)
        if isinstance(distances, np.memmap):
            distances._mmap.close()
    return {"path": str(path), "size_bytes": path.stat().st_size, "metadata": meta}
