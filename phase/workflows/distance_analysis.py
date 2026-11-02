"""Distance inspection and separately persisted distance-conditioned pair JS."""

from __future__ import annotations

from datetime import datetime, timezone
from functools import lru_cache
import json
from pathlib import Path
import tempfile
from threading import RLock
import uuid

import numpy as np

from phase.potts.metrics import js_divergence
from phase.services.project_store import ProjectStore
from phase.workflows.distances import (
    DISTANCE_CACHE_FILENAME, _frame_mask, distance_bin_indices, load_distance_cache,
    read_cluster_frames, validate_distance_alignment, validate_distance_preset,
)


_context_lock = RLock()


def _signature(path: Path):
    return (path.stat().st_mtime_ns, path.stat().st_size) if path.exists() else None


@lru_cache(maxsize=1)
def _read_context(path: str, cluster_signature, cache_signature):
    cluster = read_cluster_frames(path)
    cache = load_distance_cache(path)
    validate_distance_alignment(cluster, cache)
    with np.load(path, allow_pickle=False) as data:
        count_key = ((cluster.metadata.get("merged") or {}).get("npz_keys") or {}).get("cluster_counts", "cluster_counts")
        counts = data.get(count_key, data.get("merged__cluster_counts"))
        if counts is None:
            counts = np.maximum(cluster.labels.max(axis=0, initial=-1) + 1, 1)
        counts = np.asarray(counts, dtype=np.int64)
    if counts.shape != cluster.residue_keys.shape or np.any(counts < 0):
        raise ValueError("Cluster alphabet sizes are invalid.")
    return cluster, cache, counts


def distance_context(cluster_path):
    path = Path(cluster_path).resolve()
    # Serialize the initial large decompression and retain only the last cluster.
    # File signatures invalidate this read cache after rebuilds/reassignments.
    with _context_lock:
        return _read_context(str(path), _signature(path), _signature(path.with_name(DISTANCE_CACHE_FILENAME)))


def bin_labels(cache, thresholds):
    if not cache.available or not len(thresholds):
        return ["All distances"]
    return [f"< {thresholds[0]:g} A"] + [
        f"{a:g} <= d < {b:g} A" for a, b in zip(thresholds, thresholds[1:])
    ] + [f">= {thresholds[-1]:g} A"]


def _bins(cache, index, thresholds):
    return distance_bin_indices(
        cache.distances[:, index] if cache.available else None,
        thresholds, frame_count=len(cache.frame_indices),
    )


def _mask(cache, condition):
    condition = condition or {}
    return _frame_mask(cache, condition.get("state_ids"), condition.get("metastable_ids"))


def _joint(labels, mask, r, s, counts):
    valid = mask & (labels[:, r] >= 0) & (labels[:, s] >= 0)
    valid &= (labels[:, r] < counts[r]) & (labels[:, s] < counts[s])
    histogram = np.bincount(labels[valid, r] * counts[s] + labels[valid, s], minlength=int(counts[r] * counts[s]))
    return histogram.reshape(int(counts[r]), int(counts[s])), int(valid.sum())


def pair_detail(cluster_path, r: int, s: int, *, condition=None, thresholds=(), histogram_bins=60):
    cluster, cache, counts = distance_context(cluster_path)
    validate_distance_preset(cache, thresholds)
    r, s = sorted((int(r), int(s)))
    if r < 0 or s >= len(cluster.residue_keys) or r == s:
        raise ValueError("Select two distinct residue-column indices.")
    if not 1 <= histogram_bins <= 500:
        raise ValueError("Histogram bins must be between 1 and 500.")
    matches = np.flatnonzero((cache.edge_index[0] == r) & (cache.edge_index[1] == s))
    if cache.available and not len(matches):
        raise ValueError("Pair is not in the stored candidate set.")
    index = int(matches[0]) if len(matches) else 0
    mask = _mask(cache, condition)
    bins = _bins(cache, index, thresholds)
    labels = bin_labels(cache, thresholds)
    values = cache.distances[mask, index] if cache.available else None
    if values is not None and not np.all(np.isfinite(values)):
        raise ValueError("Missing distances in this pair; rebuild the cache.")
    hist, edges = np.histogram(values, bins=histogram_bins) if values is not None and len(values) else ([], [])
    def distribution(selected, label):
        joint, n = _joint(cluster.labels, selected, r, s, counts)
        return {"label": label, "count": n, "probabilities": (joint / max(1, n)).tolist()}
    panels = [distribution(mask, "All distances")]
    if cache.available and len(thresholds):
        panels.extend(distribution(mask & (bins == i), label) for i, label in enumerate(labels))
    occupancy = [int(np.count_nonzero(mask & (bins == i))) for i in range(len(labels))]
    overlays = []
    if values is not None and len(values):
        for state_id in dict.fromkeys(cache.frame_state_ids[mask].tolist()):
            state_hist, _ = np.histogram(cache.distances[mask & (cache.frame_state_ids == state_id), index], bins=edges)
            overlays.append({"state_id": state_id, "counts": state_hist.tolist()})
    return {
        "pair": [r, s], "residue_keys": cluster.residue_keys[[r, s]].tolist(),
        "distance_available": cache.available, "frame_count": int(mask.sum()),
        "bin_labels": labels, "bin_counts": occupancy,
        "histogram": {"counts": np.asarray(hist).tolist(), "edges": np.asarray(edges).tolist(), "by_state": overlays},
        "panels": panels, "cluster_labels_r": list(range(int(counts[r]))), "cluster_labels_s": list(range(int(counts[s]))),
    }


@lru_cache(maxsize=4)
def _summaries(path, cluster_signature, cache_signature, state_ids, metastable_ids, cutoff):
    _, cache, _ = distance_context(path)
    mask = _mask(cache, {"state_ids": state_ids, "metastable_ids": metastable_ids})
    n_frames = int(mask.sum())
    n_states = len(set(cache.frame_state_ids[mask]))
    rows = []
    for index, (r, s) in enumerate(cache.edge_index.T):
        values = cache.distances[mask, index] if cache.available else None
        summary = {"r": int(r), "s": int(s), "residue_a": str(cache.residue_keys[r]), "residue_b": str(cache.residue_keys[s]),
                   "frame_count": n_frames, "state_count": n_states}
        for key in ("min", "median", "q10", "q90", "max", "range", "contact_occupancy"):
            summary[key] = None
        if values is not None and len(values):
            finite = values[np.isfinite(values)]
            if len(finite) != len(values):
                raise ValueError("Cache contains missing distances; rebuild it.")
            q10, median, q90 = np.quantile(finite, [.1, .5, .9])
            summary.update(min=float(finite.min()), median=float(median), q10=float(q10), q90=float(q90), max=float(finite.max()),
                           range=float(np.ptp(finite)), contact_occupancy=float(np.mean(finite < cutoff)))
        rows.append(summary)
    return rows


def distance_pair_table(cluster_path, *, state_ids=None, metastable_ids=None, contact_cutoff=6.0,
                        sort="range", descending=True, search="", offset=0, limit=50,
                        min_contact_fraction=None):
    if not np.isfinite(contact_cutoff) or contact_cutoff <= 0:
        raise ValueError("Display contact cutoff must be finite and positive.")
    path = Path(cluster_path).resolve()
    rows = _summaries(str(path), _signature(path), _signature(path.with_name(DISTANCE_CACHE_FILENAME)),
                      tuple(state_ids) if state_ids is not None else None,
                      tuple(metastable_ids) if metastable_ids is not None else None, float(contact_cutoff))
    if min_contact_fraction is not None:
        if not np.isfinite(min_contact_fraction) or not 0 <= min_contact_fraction <= 1:
            raise ValueError("Minimum contact fraction must be between 0 and 1.")
        _, cache, _ = distance_context(path)
        if not cache.available:
            raise ValueError("A physical distance filter requires a distance cache.")
        rows = [row for row in rows if row["contact_occupancy"] is not None
                and row["contact_occupancy"] > 0 and row["contact_occupancy"] >= min_contact_fraction]
    if search:
        rows = [row for row in rows if search.lower() in f"{row['residue_a']} {row['residue_b']}".lower()]
    allowed = {"r", "s", "residue_a", "residue_b", "frame_count", "state_count", "min", "median", "q10", "q90", "max", "range", "contact_occupancy"}
    if sort not in allowed:
        raise ValueError("Unknown pair sort field.")
    ordered = sorted(rows, key=lambda row: (row[sort] is not None, row[sort] if row[sort] is not None else 0), reverse=descending)
    return {"total": len(rows), "pairs": ordered[max(0, offset):max(0, offset) + min(200, max(1, limit))]}


def run_distance_js_analysis(project_id, system_id, cluster_id, *, condition_a, condition_b, thresholds=(),
                             name=None, min_bin_frames=1, store=None, progress_callback=None):
    store = store if store is not None else ProjectStore()
    entry = store.get_cluster_entry(project_id, system_id, cluster_id)
    cluster_path = store.resolve_path(project_id, system_id, entry["path"])
    signature = (_signature(cluster_path), _signature(cluster_path.with_name(DISTANCE_CACHE_FILENAME)))
    cluster, cache, counts = distance_context(cluster_path)
    validate_distance_preset(cache, thresholds)
    if min_bin_frames < 1:
        raise ValueError("min_bin_frames must be positive.")
    mask_a, mask_b = _mask(cache, condition_a), _mask(cache, condition_b)
    if not mask_a.any() or not mask_b.any():
        raise ValueError("Both conditions must select at least one merged frame.")
    labels = bin_labels(cache, thresholds)
    n_bins, n_edges = len(labels), cache.edge_index.shape[1]
    if not n_edges:
        raise ValueError("No candidate/contact pairs are available for distance JS. Build the distance cache first to compute candidate pairs from the trajectories, or rebuild an empty cache with a larger candidate cutoff. The single-bin fallback requires existing contact edges.")
    js = np.full((n_bins, n_edges), np.nan)
    global_js = np.full(n_edges, np.nan)
    counts_a = np.zeros((n_bins, n_edges), dtype=np.int64)
    counts_b = np.zeros_like(counts_a)
    valid_a, valid_b = np.zeros(n_edges, dtype=np.int64), np.zeros(n_edges, dtype=np.int64)
    for e, (r, s) in enumerate(cache.edge_index.T):
        all_a, valid_a[e] = _joint(cluster.labels, mask_a, r, s, counts)
        all_b, valid_b[e] = _joint(cluster.labels, mask_b, r, s, counts)
        if valid_a[e] >= min_bin_frames and valid_b[e] >= min_bin_frames:
            global_js[e] = js_divergence((all_a / valid_a[e]).ravel(), (all_b / valid_b[e]).ravel())
        assignments = _bins(cache, e, thresholds)
        if cache.available and np.any(assignments < 0):
            raise ValueError("Missing distances in cache; rebuild before conditional JS.")
        for b in range(n_bins):
            joint_a, na = _joint(cluster.labels, mask_a & (assignments == b), r, s, counts)
            joint_b, nb = _joint(cluster.labels, mask_b & (assignments == b), r, s, counts)
            counts_a[b, e], counts_b[b, e] = na, nb
            if na >= min_bin_frames and nb >= min_bin_frames:
                js[b, e] = js_divergence((joint_a / na).ravel(), (joint_b / nb).ravel())
        if progress_callback:
            progress_callback("Computing conditional pair JS", e + 1, n_edges)
    if signature != (_signature(cluster_path), _signature(cluster_path.with_name(DISTANCE_CACHE_FILENAME))):
        raise ValueError("Cluster/distance cache changed during analysis; run again.")
    analysis_id = str(uuid.uuid4())
    root = cluster_path.parent / "analyses/distance_js"
    root.mkdir(parents=True, exist_ok=True)
    directory = root / analysis_id
    now = datetime.now(timezone.utc).isoformat()
    name_a = condition_a.get("label") or ", ".join(condition_a.get("state_ids") or condition_a.get("metastable_ids") or ["All frames"])
    name_b = condition_b.get("label") or ", ".join(condition_b.get("state_ids") or condition_b.get("metastable_ids") or ["All frames"])
    metadata = {
        "schema_version": 1, "analysis_type": "distance_js", "analysis_id": analysis_id,
        "name": name or f"{name_a} vs {name_b}", "created_at": now, "updated_at": now,
        "project_id": project_id, "system_id": system_id, "cluster_id": cluster_id,
        "condition_a": condition_a, "condition_b": condition_b,
        "state_ids": sorted(set(cache.frame_state_ids[mask_a | mask_b])),
        "metastable_ids": sorted(set(cache.frame_metastable_ids[mask_a | mask_b]) - {""}),
        "distance_available": cache.available, "distance_type": "ca_ca", "distance_unit": "angstrom",
        "requested_thresholds": list(thresholds), "thresholds": list(thresholds) if cache.available else [],
        "bin_labels": labels, "min_bin_frames": min_bin_frames, "label_mode": "assigned",
        "js_units": "nats", "js_max": float(np.log(2)), "occupancy_denominator": "pair-valid frames within each condition",
        "invalid_policy": "exclude invalid labels for this pair only", "empty_bin_policy": "NaN/null, never zero JS",
        "cache_generated_at": cache.metadata.get("generated_at"),
        "summary": {"n_edges": n_edges, "n_bins": n_bins, "frames_a": int(mask_a.sum()), "frames_b": int(mask_b.sum()),
                    "overlap_frames": int((mask_a & mask_b).sum())},
        "paths": {"analysis_npz": str((directory / "analysis.npz").relative_to(store.resolve_path(project_id, system_id, ".")))},
    }
    with tempfile.TemporaryDirectory(prefix=".distance-js-", dir=root) as work:
        work = Path(work)
        np.savez_compressed(work / "analysis.npz", edge_index=cache.edge_index,
            edge_labels=np.asarray([f"{cache.residue_keys[r]} - {cache.residue_keys[s]}" for r, s in cache.edge_index.T]),
            distance_bin_labels=np.asarray(labels), edge_bin_js=js, edge_global_js=global_js,
            edge_bin_counts_a=counts_a, edge_bin_counts_b=counts_b, edge_valid_counts_a=valid_a, edge_valid_counts_b=valid_b,
            edge_bin_occupancy_a=counts_a / np.maximum(1, valid_a)[None, :],
            edge_bin_occupancy_b=counts_b / np.maximum(1, valid_b)[None, :], residue_keys=cluster.residue_keys)
        (work / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        # Unique directory publication keeps incomplete results out of listings.
        work.rename(directory)
        work.mkdir()
    return {"metadata": metadata, "analysis_npz": str(directory / "analysis.npz")}
