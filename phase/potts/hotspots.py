"""Gauge-fixed, mask-aware endpoint discrimination; no mutation-causality claim."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import uuid
import tempfile
import warnings

import numpy as np

from phase.potts.masking import observed_labels
from phase.potts.potts_model import PottsModel, load_potts_model, load_potts_model_metadata, zero_sum_gauge_model
from phase.potts.sample_io import load_sample_npz
from phase.services.project_store import _write_json

ANALYSIS_TYPE = "model_pair_hotspots"


def contrast_models(a: PottsModel, b: PottsModel) -> PottsModel:
    """RMS norms of this contrast are alphabet-normalized, not energy signs."""
    if a.K_list() != b.K_list() or not a.h:
        raise ValueError("Endpoint models must have identical nonempty residue alphabets.")
    for model in (a, b):
        for h in model.h:
            if h.ndim != 1 or not h.size or not np.isfinite(h).all():
                raise ValueError("Invalid model fields.")
        if set(model.edges) != set(model.J) or len(model.edges) != len(set(model.edges)):
            raise ValueError("Model edges and couplings disagree.")
        for (r, s), matrix in model.J.items():
            if not 0 <= r < s < len(model.h) or matrix.shape != (len(model.h[r]), len(model.h[s])) or not np.isfinite(matrix).all():
                raise ValueError("Invalid model coupling.")
    a, b = zero_sum_gauge_model(a), zero_sum_gauge_model(b)
    edges = sorted(set(a.edges) | set(b.edges))
    return PottsModel([x - y for x, y in zip(a.h, b.h)],
                      {edge: a.J.get(edge, 0) - b.J.get(edge, 0) for edge in edges}, edges)


def finite_mean(values, axis=0):
    values = np.asarray(values, dtype=float)
    count = np.isfinite(values).sum(axis=axis)
    return np.divide(np.nansum(values, axis=axis), count, out=np.full(np.shape(count), np.nan), where=count > 0)


def block_statistics(model, labels, frame_indices=None, frame_state_ids=None, block_length=50):
    """Nonoverlapping blocks never cross state boundaries or irregular frame gaps.

    Each block holds term sums/counts, not frame-by-edge energy arrays. Missing
    coordinates are absent observations, not zero-valued observations.
    """
    x, observed = observed_labels(labels, model.K_list())
    n = len(x)
    if not n or block_length < 1:
        raise ValueError("A nonempty sample and positive block length are required.")
    frames = np.arange(n) if frame_indices is None else np.asarray(frame_indices)
    states = np.full(n, "") if frame_state_ids is None else np.asarray(frame_state_ids)
    if frames.shape != (n,) or states.shape != (n,):
        raise ValueError("Frame provenance must match labels.")
    diffs = np.diff(frames)
    positive = diffs[diffs > 0]
    stride = int(np.min(positive)) if positive.size else 1
    breaks = np.flatnonzero((diffs != stride) | (states[1:] != states[:-1])) + 1
    boundaries = np.r_[0, breaks, n]
    slices = [(i, min(i + block_length, end)) for start, end in zip(boundaries[:-1], boundaries[1:])
              for i in range(start, end, block_length)]
    n_terms = len(model.h) + len(model.edges)
    sums = np.zeros((len(slices), n_terms))
    counts = np.zeros_like(sums, dtype=np.int64)
    starts = np.array([start for start, _ in slices], dtype=int)
    # One vectorized frame-to-block reduction per term; never allocate T x E.
    for r, h in enumerate(model.h):
        mask = observed[:, r]
        sums[:, r] = np.add.reduceat(np.where(mask, h[x[:, r]], 0.), starts)
        counts[:, r] = np.add.reduceat(mask.astype(np.int64), starts)
    for j, (r, s) in enumerate(model.edges, len(model.h)):
        mask = observed[:, r] & observed[:, s]
        sums[:, j] = np.add.reduceat(np.where(mask, model.J[r, s][x[:, r], x[:, s]], 0.), starts)
        counts[:, j] = np.add.reduceat(mask.astype(np.int64), starts)
    return sums, counts


def summarize_blocks(blocks, bootstrap=200, seed=0):
    """Equal trajectory weight; resample blocks independently within trajectories.

    Chunking bootstrap matrix products bounds peak memory even for large graphs.
    CIs require >=4 observed blocks in total and >=2 in every contributing
    trajectory. A single structure therefore never gets a spurious narrow CI.
    """
    if not blocks:
        raise ValueError("Select at least one reference MD sample per endpoint.")
    rng = np.random.default_rng(seed)
    means, draws, coverage, support, reliable = [], [], [], [], []
    for sums, counts in blocks:
        total = counts.sum(axis=0)
        means.append(np.divide(sums.sum(axis=0), total, out=np.full(len(total), np.nan), where=total > 0))
        coverage.append(total)
        supported = (counts > 0).sum(axis=0)
        support.append(supported)
        reliable.append((total == 0) | (supported >= 2))
        result = np.full((bootstrap, len(total)), np.nan)
        for start in range(0, bootstrap, 32):
            weights = rng.multinomial(len(sums), np.full(len(sums), 1 / len(sums)), size=min(32, bootstrap - start))
            denominators = weights @ counts
            np.divide(weights @ sums, denominators, out=result[start:start + len(weights)], where=denominators > 0)
        draws.append(result)
    return {"mean": finite_mean(means), "draws": finite_mean(draws),
            "counts": np.sum(coverage, axis=0), "blocks": np.sum(support, axis=0),
            "reliable": np.all(reliable, axis=0) & (np.sum(support, axis=0) >= 4)}


def separation_statistics(a, b, bootstrap=200, seed=0):
    if bootstrap < 0 or (bootstrap and bootstrap < 50):
        raise ValueError("Use zero (no confidence intervals) or at least 50 bootstrap replicates.")
    sa = summarize_blocks(a, bootstrap, seed)
    sb = summarize_blocks(b, bootstrap, seed + 1)
    separation = sb["mean"] - sa["mean"]
    low, high = np.full_like(separation, np.nan), np.full_like(separation, np.nan)
    reliable = sa["reliable"] & sb["reliable"] & np.isfinite(separation)
    if bootstrap and reliable.any():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            intervals = np.nanquantile(sb["draws"] - sa["draws"], [.025, .975], axis=0)
        low[reliable], high[reliable] = intervals[:, reliable]
    return {"mean_a": sa["mean"], "mean_b": sb["mean"], "separation": separation,
            "ci_low": low, "ci_high": high, "counts_a": sa["counts"], "counts_b": sb["counts"],
            "blocks_a": sa["blocks"], "blocks_b": sb["blocks"]}


def is_full_model(entry):
    params = entry.get("params") or {}
    return params.get("delta_kind") != "delta_patch" and not (params.get("fit_mode") == "delta" and not params.get("delta_kind"))


def suggested_references(model, samples):
    params = model.get("params") or {}
    state_ids = set(map(str, params.get("state_ids") or []))
    return [s["sample_id"] for s in samples if s.get("type") == "md_eval"
            and set(map(str, (s.get("params") or {}).get("state_ids") or [s.get("state_id", "")])) & state_ids]


def run_hotspot_analysis(project_id, system_id, cluster_id, *, store, model_a_id, model_b_id,
                         mode="ensemble", reference_a=None, reference_b=None, sample_ids=None,
                         block_length=50, bootstrap=200, seed=0, progress_callback=None):
    if mode not in ("parameter", "ensemble") or model_a_id == model_b_id:
        raise ValueError("Select different endpoint models and parameter or ensemble mode.")
    if not isinstance(block_length, int) or block_length < 1:
        raise ValueError("Block length must be a positive integer.")
    if not isinstance(bootstrap, int) or bootstrap < 0 or (0 < bootstrap < 50):
        raise ValueError("Bootstrap replicates must be zero or >=50.")
    if progress_callback:
        progress_callback("Validating endpoint models", 0, 1)
    entry = store.get_cluster_entry(project_id, system_id, cluster_id)
    cluster_path = store.resolve_path(project_id, system_id, entry["path"])
    if cluster_path.is_dir():
        cluster_path = cluster_path / "cluster.npz"
    cluster_dir = cluster_path.parent
    with np.load(cluster_path, allow_pickle=False) as data:
        keys = data["residue_keys"].astype(str)
    models = {m["model_id"]: m for m in store.list_potts_models(project_id, system_id, cluster_id)}
    resolved, provenance = [], []
    for identifier in (model_a_id, model_b_id):
        meta = models.get(identifier)
        if not meta or not is_full_model(meta):
            raise ValueError(f"Select a full/combined endpoint model, not a raw delta patch: {identifier}")
        if meta.get("residue_keys") is not None and list(meta["residue_keys"]) != keys.tolist():
            raise ValueError("Model canonical residue ordering differs from this cluster.")
        path = store.resolve_path(project_id, system_id, meta["path"])
        model_bytes = path.read_bytes()
        embedded = load_potts_model_metadata(io.BytesIO(model_bytes)) or {}
        for source in (embedded, embedded.get("params") or {}, meta.get("params") or {}):
            if source.get("residue_keys") is not None and list(source["residue_keys"]) != keys.tolist():
                raise ValueError("Model canonical residue ordering differs from this cluster.")
            if source.get("cluster_id") and str(source["cluster_id"]) != cluster_id:
                raise ValueError("Model belongs to a different cluster categorical basis.")
        if not is_full_model({"params": embedded.get("params") or embedded}):
            raise ValueError("The model NPZ is a raw delta patch, not a full endpoint.")
        resolved.append(load_potts_model(io.BytesIO(model_bytes)))
        provenance.append({"model_id": identifier, "name": meta.get("name") or identifier,
                           "path": meta["path"], "sha256": hashlib.sha256(model_bytes).hexdigest()})
    delta = contrast_models(*resolved)
    n = len(keys)
    if len(delta.h) != n:
        raise ValueError("Model residue count differs from the selected cluster.")
    arrays = {"residue_keys": keys, "edges": np.asarray(delta.edges, dtype=int).reshape(-1, 2),
              "node_parameter_rms": np.array([np.sqrt(np.mean(h ** 2)) for h in delta.h]),
              "edge_parameter_rms": np.array([np.sqrt(np.mean(delta.J[e] ** 2)) for e in delta.edges])}
    samples = {s["sample_id"]: s for s in store.list_samples(project_id, system_id, cluster_id)}
    refs_a = list(dict.fromkeys(reference_a if reference_a is not None else suggested_references(models[model_a_id], samples.values())))
    refs_b = list(dict.fromkeys(reference_b if reference_b is not None else suggested_references(models[model_b_id], samples.values())))
    selected = list(dict.fromkeys((refs_a + refs_b + (sample_ids or [])) if mode == "ensemble" else []))
    if mode == "ensemble" and (not refs_a or not refs_b or set(refs_a) & set(refs_b)):
        raise ValueError("Select nonempty, disjoint A and B reference MD samples explicitly (inference can be ambiguous).")
    stats, sample_meta, seen_frames = {}, [], {}
    means_by_sample, counts_by_sample = {}, {}
    for i, identifier in enumerate(selected):
        if progress_callback:
            progress_callback("Evaluating trajectory blocks", i, len(selected) + 1)
        sample_entry = samples.get(identifier)
        if sample_entry is None:
            raise ValueError(f"Sample not found in this cluster: {identifier}")
        is_reference = identifier in refs_a or identifier in refs_b
        if is_reference and sample_entry.get("type") != "md_eval":
            raise ValueError("Endpoint reference samples must be assigned MD trajectories.")
        paths = sample_entry.get("paths") or {}
        rel = paths.get("summary_npz") or sample_entry.get("path")
        if not rel:
            raise ValueError(f"Sample has no NPZ path: {identifier}")
        path = store.resolve_path(project_id, system_id, rel)
        if not path.exists() and (cluster_dir / rel).is_file():
            path = cluster_dir / rel
        before = path.stat()
        sample = load_sample_npz(path)
        keep = np.ones(len(sample.labels), dtype=bool)
        if sample.invalid_mask is not None:
            keep &= ~np.asarray(sample.invalid_mask, dtype=bool)
        # Ordinary halo labels are not missing coordinates. Exclude such frames,
        # explicitly reporting their count; -2 remains a masked observation.
        keep &= ~np.any(sample.labels == -1, axis=1)
        frames = np.arange(len(keep)) if sample.frame_indices is None else sample.frame_indices
        states = sample.frame_state_ids
        sid = sample_entry.get("state_id") or (sample_entry.get("params") or {}).get("state_id") or identifier
        if states is None:
            states = np.full(len(keep), sid)
        else:
            states = np.where(np.asarray(states).astype(str) == "", sid, states)
        if is_reference:
            if len(np.unique(states[keep])) != 1:
                raise ValueError("Each reference MD sample must belong to one source state; assign/select separate state samples.")
            group = "a" if identifier in refs_a else "b"
            for state, frame in zip(states[keep], frames[keep]):
                identity = (str(state), int(frame))
                if identity in seen_frames:
                    raise ValueError("Reference MD samples overlap in source state/frame; select independent, non-overlapping references.")
                seen_frames[identity] = group
        if is_reference:
            stats[identifier] = block_statistics(delta, sample.labels[keep], frames[keep], states[keep], block_length)
        else:
            # Preference-only ensembles need a mean, not bootstrap blocks. In
            # particular, random SA source-frame IDs must not create T x E blocks.
            stats[identifier] = block_statistics(delta, sample.labels[keep], block_length=max(1, int(keep.sum())))
        sums, counts = stats[identifier]
        totals = counts.sum(axis=0)
        means_by_sample[identifier] = np.divide(sums.sum(axis=0), totals, out=np.full(len(totals), np.nan), where=totals > 0)
        counts_by_sample[identifier] = totals
        if not is_reference:
            del stats[identifier]
        sample_hash = hashlib.sha256(path.read_bytes()).hexdigest()
        after = path.stat()
        if (before.st_ino, before.st_size, before.st_mtime_ns) != (after.st_ino, after.st_size, after.st_mtime_ns):
            raise ValueError(f"Sample changed during analysis; wait for assignment/sampling to finish and retry: {identifier}")
        sample_meta.append({"sample_id": identifier, "name": sample_entry.get("name") or identifier,
                            "path": str(rel), "frames": int(keep.sum()), "excluded_frames": int((~keep).sum()),
                            "sha256": sample_hash})
    if mode == "ensemble":
        if progress_callback:
            progress_callback("Block bootstrap uncertainty", len(selected), len(selected) + 1)
        summary = separation_statistics([stats[s] for s in refs_a], [stats[s] for s in refs_b], bootstrap, seed)
        for name, values in summary.items():
            arrays[f"node_{name}"], arrays[f"edge_{name}"] = values[:n], values[n:]
        sample_means = np.array([means_by_sample[s] for s in selected])
        arrays["sample_ids"] = np.array(selected, dtype=str)
        arrays["node_preference"] = .5 * (summary["mean_a"][:n] + summary["mean_b"][:n]) - sample_means[:, :n]
        arrays["edge_preference"] = .5 * (summary["mean_a"][n:] + summary["mean_b"][n:]) - sample_means[:, n:]
        arrays["sample_node_counts"] = np.array([counts_by_sample[s][:n] for s in selected])
        arrays["sample_edge_counts"] = np.array([counts_by_sample[s][n:] for s in selected])
    identifier = str(uuid.uuid4())
    destination = cluster_dir / "analyses" / ANALYSIS_TYPE / identifier
    metadata = {"analysis_id": identifier, "analysis_type": ANALYSIS_TYPE, "version": 1,
                "created_at": datetime.now(timezone.utc).isoformat(), "cluster_id": cluster_id,
                "model_a_id": model_a_id, "model_b_id": model_b_id, "model_ids": [model_a_id, model_b_id], "models": provenance,
                "name": f"{provenance[0]['name']} / {provenance[1]['name']}", "mode": mode,
                "reference_sample_ids_a": refs_a if mode == "ensemble" else [], "reference_sample_ids_b": refs_b if mode == "ensemble" else [],
                "sample_ids": selected, "samples": sample_meta, "gauge": "zero_sum", "edge_basis": "union",
                "weighting": "equal_trajectory", "block_length": block_length, "bootstrap": bootstrap, "seed": seed,
                "uncertainty": "95% within-trajectory nonoverlapping block bootstrap; fixed model uncertainty only",
                "warnings": ["Verify reference independence and block-length sensitivity. CIs exclude model-fitting uncertainty.",
                             "Legacy model files without residue keys assume the selected cluster categorical basis.",
                             "Scores prioritize discriminative candidates; they do not predict mutation causality."],
                "paths": {"summary_npz": str((destination / "analysis.npz").relative_to(store._system_dir(project_id, system_id)))}}
    # Commit the complete artifact atomically; health cleanup never sees a
    # half-written analysis and failed writes leave no listed analysis behind.
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".hotspots-", dir=cluster_dir) as temporary:
        staging = Path(temporary) / identifier
        staging.mkdir()
        np.savez_compressed(staging / "analysis.npz", **arrays)
        _write_json(staging / "analysis_metadata.json", metadata)
        staging.replace(destination)
    if progress_callback:
        progress_callback("Complete", len(selected) + 1, len(selected) + 1)
    return metadata
