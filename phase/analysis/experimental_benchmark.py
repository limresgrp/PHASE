"""Quantitative static-structure benchmark, with no biological classifications."""
from __future__ import annotations

from datetime import datetime, timezone
import csv
import hashlib
import io
import itertools
import json
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path
import subprocess
import tempfile
import uuid
import warnings

import numpy as np
from scipy.stats import rankdata, spearmanr, pearsonr
from scipy.spatial.distance import pdist, squareform

from phase.potts.hotspots import contrast_models, is_full_model, block_statistics, separation_statistics, suggested_references
from phase.potts.masking import observed_labels
from phase.potts.potts_model import PottsModel, load_potts_model, load_potts_model_metadata
from phase.potts.sample_io import load_sample_npz
from phase.services.analysis_annotations import get_scheme, validate_scheme, list_schemes
from phase.services.project_store import _write_json

ANALYSIS_TYPE = 'experimental_structure_benchmark'
DISCLAIMER = 'Quantitative comparisons only. User labels are annotations, not inferred biological meaning. Ranks are not probabilities; hotspots are not causal or mutational predictions.'


def clean_json(value):
    if isinstance(value, dict):
        return {str(k): clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean_json(v) for v in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    return value


def percentiles(scores):
    scores = np.asarray(scores, dtype=float)
    valid = np.isfinite(scores)
    result = np.full(scores.shape, np.nan)
    n = valid.sum()
    if n:
        result[valid] = (rankdata(scores[valid], method='average') - 1) / (n - 1) if n > 1 else .5
    return result


def correlation(x, y, method='spearman'):
    x, y = np.asarray(x, float), np.asarray(y, float)
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    n = len(x)
    if n < 2 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return {'coefficient': np.nan, 'p': np.nan, 'n': n}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        coefficient, p = (spearmanr if method == 'spearman' else pearsonr)(x, y)
    return {'coefficient': float(coefficient), 'p': float(p) if n > 2 else np.nan, 'n': n}


def bh_adjust(p_values):
    p = np.asarray(p_values, dtype=float)
    q = np.full(p.shape, np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    order = valid[np.argsort(p[valid])]
    if len(order):
        values = p[order] * len(order) / np.arange(1, len(order) + 1)
        q[order] = np.minimum(1, np.minimum.accumulate(values[::-1])[::-1])
    return q


def signed_auc(scores, labels, positive, negative, best_oriented=False):
    if positive == negative:
        raise ValueError('Positive and negative AUC labels must differ.')
    scores = np.asarray(scores, float)
    labels = np.asarray(labels, object)
    keep = np.isfinite(scores) & np.isin(labels, [positive, negative])
    values, targets = scores[keep], labels[keep] == positive
    n_pos, n_neg = int(targets.sum()), int((~targets).sum())
    auc = np.nan
    if n_pos and n_neg:
        auc = float((rankdata(values)[targets].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))
    return {'auc': auc, 'n_positive': n_pos, 'n_negative': n_neg, 'n': n_pos + n_neg,
            'best_oriented_auc': max(auc, 1 - auc) if best_oriented and np.isfinite(auc) else None,
            'positive': positive, 'negative': negative, 'direction': 'higher E_A - E_B is positive; no automatic sign flip'}


def pcoa(distances):
    d = np.asarray(distances, float)
    n = len(d)
    if d.shape != (n, n) or not np.isfinite(d).all():
        return {'available': False, 'reason': 'Distance matrix contains undefined entries; no imputation was performed.'}
    if np.any(d < -1e-12) or not np.allclose(d, d.T) or not np.allclose(np.diag(d), 0):
        raise ValueError('PCoA requires a symmetric nonnegative distance matrix with zero diagonal.')
    center = np.eye(n) - np.ones((n, n)) / n
    gram = -.5 * center @ (d ** 2) @ center
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]
    tolerance = max(1., np.max(np.abs(eigenvalues), initial=0)) * 1e-12
    positive = eigenvalues > tolerance
    total = eigenvalues[positive].sum()
    coords = np.zeros((n, 3))
    fractions = np.zeros(3)
    n_axes = min(3, int(positive.sum()))
    if n_axes:
        coords[:, :n_axes] = eigenvectors[:, :n_axes] * np.sqrt(eigenvalues[:n_axes])
        # Stabilize arbitrary eigenvector signs, not their scientific meaning.
        for i in range(n_axes):
            if coords[np.argmax(np.abs(coords[:, i])), i] < 0:
                coords[:, i] *= -1
        fractions[:n_axes] = eigenvalues[:n_axes] / total
    absolute_mass = np.abs(eigenvalues).sum()
    return {'available': True, 'coordinates': coords, 'positive_eigenvalue_fraction': fractions,
            'eigenvalues': eigenvalues, 'negative_eigenvalue_fraction': float(np.abs(eigenvalues[eigenvalues < -tolerance]).sum() / absolute_mass) if absolute_mass else 0.,
            'degenerate': not bool(n_axes)}


def term_values(delta, labels, all_edges):
    """Float64 version of the shared masked delta-energy decomposition.

    Missing coordinates remain NaN. Absent model couplings are mathematical
    zeros, which are distinct from absent structural observations.
    """
    from phase.potts.analysis_run import _compute_delta_energy_components
    _, observed = observed_labels(labels, delta.K_list())
    zero = PottsModel([np.zeros_like(h) for h in delta.h], {}, [])
    node, local_values, local_edges = _compute_delta_energy_components(
        labels, delta, zero, dtype=np.float64, mask_missing=True)
    positions = {tuple(edge): i for i, edge in enumerate(local_edges.tolist())}
    edge = np.empty((len(labels), len(all_edges)))
    for j, (r, s) in enumerate(all_edges):
        if (r, s) in positions:
            edge[:, j] = local_values[:, positions[r, s]]
        else:
            edge[:, j] = np.where(observed[:, r] & observed[:, s], 0., np.nan)
    return np.concatenate([node, edge], axis=1)


def score_terms(terms, support):
    selected = np.asarray(terms)[:, np.asarray(support, bool)]
    if selected.shape[1] == 0:
        return np.full(len(terms), np.nan)
    result = np.nansum(selected, axis=1)
    result[~np.isfinite(selected).any(axis=1)] = np.nan
    return result


def compute_benchmark(models, contrasts, labels, state_ids, scheme, *, fingerprint_mode='standardized',
                      distance='euclidean', positive_label=None, negative_label=None, best_oriented=False,
                      bh=True, global_threshold=.15, local_threshold=None):
    """Pure computation, independent from ProjectStore and uploaded structures."""
    if len(state_ids) < 2 or len(set(state_ids)) != len(state_ids):
        raise ValueError('Select at least two distinct experimental/static states.')
    if fingerprint_mode not in ('standardized', 'raw') or distance not in ('euclidean', 'correlation'):
        raise ValueError('Choose standardized/raw fingerprints and euclidean/correlation distance.')
    if not 0 <= global_threshold <= 1 or (local_threshold is not None and (not np.isfinite(local_threshold) or local_threshold < 0)):
        raise ValueError('Discordance thresholds must be finite and nonnegative; global threshold is in [0,1].')
    if not contrasts or len(set(tuple(c) for c in contrasts)) != len(contrasts):
        raise ValueError('Select distinct ordered model contrasts.')
    deltas = []
    for a, b in contrasts:
        if a == b or a not in models or b not in models:
            raise ValueError('Each contrast requires two different selected endpoint models.')
        deltas.append(contrast_models(models[a], models[b]))
    k = deltas[0].K_list()
    if any(d.K_list() != k for d in deltas):
        raise ValueError('All axes must share one cluster categorical basis.')
    x, observed = observed_labels(labels, k)
    if len(x) != len(state_ids):
        raise ValueError('Exactly one assigned frame is required per selected state.')
    edges = sorted(set().union(*(set(d.edges) for d in deltas)))
    edge_observed = np.column_stack([observed[:, r] & observed[:, s] for r, s in edges]) if edges else np.zeros((len(x), 0), bool)
    observed_terms = np.column_stack([observed, edge_observed])
    common = observed_terms.all(axis=0)
    terms = np.stack([term_values(d, labels, edges) for d in deltas])
    n = len(k)
    axis_support = np.array([np.r_[np.ones(n, bool), [e in d.J for e in edges]] for d in deltas], dtype=bool)
    raw = np.stack([score_terms(t, support) for t, support in zip(terms, axis_support)], axis=1)
    common_scores = np.stack([score_terms(t, common & support) for t, support in zip(terms, axis_support)], axis=1)
    raw_ranks = np.column_stack([percentiles(v) for v in raw.T])
    common_ranks = np.column_stack([percentiles(v) for v in common_scores.T])
    annotation = [scheme['labels'].get(s) for s in state_ids]
    ordinal = {label: i for i, label in enumerate(scheme['categories'])} if scheme['kind'] == 'ordinal' else None
    target = np.array([ordinal.get(label, np.nan) for label in annotation]) if ordinal else None
    statistics, qc = [], []
    for i, delta in enumerate(deltas):
        for variant, matrix in [('raw', raw), ('common', common_scores)]:
            row = {'axis': i, 'variant': variant}
            if target is not None:
                row['spearman'] = correlation(matrix[:, i], target)
            if positive_label is not None or negative_label is not None:
                if positive_label not in scheme['categories'] or negative_label not in scheme['categories']:
                    raise ValueError('Select declared positive and negative AUC categories.')
                row['auc'] = signed_auc(matrix[:, i], annotation, positive_label, negative_label, best_oriented)
            statistics.append(row)
        edge_counts = edge_observed[:, [e in delta.J for e in edges]].sum(axis=1)
        node_counts = observed.sum(axis=1)
        combined_coverage = (node_counts + edge_counts) / (n + len(delta.edges))
        coverage_correlation = correlation(raw[:, i], combined_coverage)
        qc.append({'axis': i, 'node_counts': node_counts, 'edge_counts': edge_counts,
                   'node_fraction': node_counts / n, 'edge_fraction': edge_counts / len(delta.edges) if delta.edges else np.full(len(x), np.nan),
                   'common_nodes': int(common[:n].sum()), 'common_edges': int((common[n:] & axis_support[i, n:]).sum()),
                   'raw_vs_coverage': coverage_correlation,
                   'raw_vs_node_coverage': correlation(raw[:, i], node_counts),
                   'raw_vs_edge_coverage': correlation(raw[:, i], edge_counts),
                   'raw_vs_common_ranks': correlation(raw_ranks[:, i], common_ranks[:, i]),
                   'coverage_warning': bool(abs(coverage_correlation['coefficient']) >= .7)})
    # Two predeclared BH families: all raw axes, and all common-support axes.
    if target is not None:
        for variant in ('raw', 'common'):
            family = [row for row in statistics if row['variant'] == variant]
            adjusted = bh_adjust([row['spearman']['p'] for row in family])
            for row, q in zip(family, adjusted):
                row['spearman']['q'] = q if bh else None
    redundancy = {}
    for name, matrix, method in [('rank_spearman', common_ranks, 'spearman'), ('common_score_pearson', common_scores, 'pearson')]:
        details = [[correlation(a, b, method) for b in matrix.T] for a in matrix.T]
        redundancy[name] = {'matrix': [[row['coefficient'] for row in detail] for detail in details],
                            'n': [[row['n'] for row in detail] for detail in details]}
    algebra = []
    # All axes use the same structural support; absent couplings are zero.
    global_common_scores = np.stack([score_terms(t, common) for t in terms], axis=1)
    lookup_table, neighbors = {}, {}
    for i, (a, b) in enumerate(contrasts):
        lookup_table.setdefault((a, b), (i, 1))
        lookup_table.setdefault((b, a), (i, -1))
        neighbors.setdefault(a, set()).add(b)
        neighbors.setdefault(b, set()).add(a)
    def lookup(a, b):
        return lookup_table.get((a, b))
    # Enumerate actual graph triangles, not all O(number_of_models**3) triples.
    triangles = ((a, b, c) for a in sorted(neighbors) for b in sorted(neighbors[a]) if b > a
                 for c in sorted(neighbors[a] & neighbors[b]) if c > b)
    for a, b, c in triangles:
        selections = [lookup(a, b), lookup(b, c), lookup(a, c)]
        if any(s is None for s in selections):
            continue
        arrays = [sign * global_common_scores[:, i] for i, sign in selections]
        residual = arrays[0] + arrays[1] - arrays[2]
        finite = residual[np.isfinite(residual)]
        max_residual = np.max(np.abs(finite)) if len(finite) else np.nan
        scale = max(1., *(np.max(np.abs(v), initial=0) for v in arrays))
        algebra.append({'models': [a, b, c], 'identity': '(A-B)+(B-C)=(A-C)', 'support': 'global common node/edge terms',
                        'n': len(finite), 'max_residual': max_residual, 'mean_residual': np.mean(np.abs(finite)) if len(finite) else np.nan,
                        'tolerance': 1e-10 * scale, 'passed': bool(max_residual <= 1e-10 * scale) if len(finite) else None})
    # Each axis includes its own terms on the globally fixed common support.
    features, feature_axis, feature_term = [], [], []
    for i, support in enumerate(axis_support):
        selected = np.flatnonzero(common & support)
        features.append(terms[i][:, selected])
        feature_axis.extend([i] * len(selected)); feature_term.extend(selected.tolist())
    raw_fp = np.concatenate(features, axis=1)
    mean, std = raw_fp.mean(axis=0), raw_fp.std(axis=0)
    standardized = np.divide(raw_fp - mean, std, out=np.zeros_like(raw_fp), where=std > 1e-12)
    chosen = standardized if fingerprint_mode == 'standardized' else raw_fp
    euclidean = squareform(pdist(standardized, 'euclidean')) if raw_fp.shape[1] else np.full((len(x), len(x)), np.nan)
    np.fill_diagonal(euclidean, 0)
    corr_dist = np.full((len(x), len(x)), np.nan)
    np.fill_diagonal(corr_dist, 0)
    for i, j in itertools.combinations(range(len(x)), 2):
        if chosen.shape[1] and np.ptp(chosen[i]) and np.ptp(chosen[j]):
            corr_dist[i, j] = corr_dist[j, i] = np.clip(1 - np.corrcoef(chosen[i], chosen[j])[0, 1], 0, 2)
    selected_distance = euclidean if distance == 'euclidean' else corr_dist
    embedding = pcoa(selected_distance)
    finite_local = selected_distance[np.triu_indices(len(x), 1)]
    finite_local = finite_local[np.isfinite(finite_local)]
    local_cutoff = local_threshold if local_threshold is not None else (float(np.quantile(finite_local, .75)) if finite_local.size else np.nan)
    discordance = []
    for i, j in itertools.combinations(range(len(x)), 2):
        differences = common_ranks[i] - common_ranks[j]
        global_difference = np.sqrt(np.mean(differences ** 2)) if np.isfinite(differences).all() else np.nan
        local = selected_distance[i, j]
        discordance.append({'state_a': state_ids[i], 'state_b': state_ids[j], 'global_rank_difference': global_difference,
                            'local_distance': local, 'threshold_match': bool(global_difference <= global_threshold and local >= local_cutoff) if np.isfinite(local) and np.isfinite(global_difference) else None})
    result = {'state_ids': state_ids, 'labels': annotation, 'statistics': statistics, 'coverage_qc': qc,
              'axis_redundancy': redundancy, 'axis_algebra': algebra, 'pcoa': embedding,
              'discordance': discordance, 'global_threshold': global_threshold, 'local_threshold': local_cutoff,
              'disclaimer': DISCLAIMER}
    arrays = {'state_ids': np.array(state_ids, str), 'edges': np.asarray(edges, int).reshape(-1, 2),
              'raw_scores': raw, 'common_scores': common_scores, 'raw_ranks': raw_ranks, 'common_ranks': common_ranks,
              'node_observed': observed, 'edge_observed': edge_observed, 'common_support': common,
              'axis_support': axis_support, 'term_contributions': terms, 'fingerprints_raw': raw_fp,
              'fingerprints_standardized': standardized, 'fingerprint_mean': mean, 'fingerprint_std': std,
              'feature_axis': np.array(feature_axis), 'feature_term': np.array(feature_term),
              'euclidean_distance': euclidean, 'correlation_distance': corr_dist,
              'rank_spearman_matrix': np.asarray(redundancy['rank_spearman']['matrix']),
              'common_score_pearson_matrix': np.asarray(redundancy['common_score_pearson']['matrix'])}
    if embedding['available']:
        arrays['pcoa_coordinates'] = embedding['coordinates']
        arrays['pcoa_eigenvalues'] = embedding['eigenvalues']
    return result, arrays, deltas


def _sample_path(store, project, system, entry):
    rel = (entry.get('paths') or {}).get('summary_npz') or entry.get('path')
    if not rel:
        raise ValueError(f"Assigned sample has no NPZ: {entry.get('sample_id')}")
    return store.resolve_path(project, system, rel)


def _read_sample(store, project, system, entry):
    path = _sample_path(store, project, system, entry)
    before = path.stat()
    sample = load_sample_npz(path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    after = path.stat()
    if (before.st_ino, before.st_mtime_ns, before.st_size) != (after.st_ino, after.st_mtime_ns, after.st_size):
        raise ValueError('Assignment changed during analysis. Wait for assignment to finish and retry.')
    try:
        recorded_path = str(path.relative_to(store._system_dir(project, system)))
    except ValueError:
        recorded_path = str(path)
    return sample, {'sample_id': entry['sample_id'], 'path': recorded_path, 'sha256': digest}


def hotspot_comparisons(items, top_k):
    result = []
    for left, right in itertools.combinations(items, 2):
        for kind in ('nodes', 'edges'):
            a = {row['identity']: row['score'] for row in left[kind] if np.isfinite(row['score'])}
            b = {row['identity']: row['score'] for row in right[kind] if np.isfinite(row['score'])}
            def top(values):
                cutoff = sorted((abs(v) for v in values.values()), reverse=True)[min(top_k, len(values)) - 1] if values else np.inf
                return {key for key, value in values.items() if abs(value) >= cutoff}
            ta, tb = top(a), top(b)
            shared = sorted(set(a) & set(b))
            result.append({'axis_a': left['axis'], 'axis_b': right['axis'], 'kind': kind, 'top_k': top_k,
                           'top_a_count': len(ta), 'top_b_count': len(tb), 'overlap': len(ta & tb),
                           'jaccard': len(ta & tb) / len(ta | tb) if ta | tb else np.nan,
                           'rank_correlation': correlation([abs(a[k]) for k in shared], [abs(b[k]) for k in shared]),
                           'ties': 'All ties at the kth cutoff included; actual set size may exceed k.'})
    return result


def run_experimental_benchmark(project_id, system_id, cluster_id, *, store, state_ids, contrasts,
                               scheme_id=None, label_scheme=None, state_samples=None,
                               hotspots_mode='parameter', top_k=20, block_length=50, bootstrap=200, seed=0,
                               assign_missing=True, name=None, progress_callback=None, **options):
    if hotspots_mode not in ('parameter', 'ensemble') or top_k < 1 or seed < 0 or block_length < 1 or bootstrap < 0 or (0 < bootstrap < 50):
        raise ValueError('Invalid hotspot mode, top-k, seed, block length or bootstrap count (0 or >=50).')
    defaults = {'fingerprint_mode': 'standardized', 'distance': 'euclidean', 'positive_label': None,
                'negative_label': None, 'best_oriented': False, 'bh': True,
                'global_threshold': .15, 'local_threshold': None}
    unknown = set(options) - set(defaults)
    if unknown:
        raise ValueError(f"Unknown benchmark options: {sorted(unknown)}")
    options = {**defaults, **options}
    system = store.get_system(project_id, system_id)
    if len(state_ids) < 2 or len(set(state_ids)) != len(state_ids) or not set(state_ids) <= set(system.states):
        raise ValueError('Select at least two distinct existing static states.')
    if bool(scheme_id) == bool(label_scheme):
        raise ValueError('Supply one saved label scheme ID or one explicit label scheme snapshot.')
    if scheme_id:
        scheme = get_scheme(store, project_id, system_id, scheme_id)
    else:
        validation_payload = label_scheme
        if label_scheme.get('scheme_id'):
            validation_payload = {**label_scheme, 'labels': {s: v for s, v in label_scheme['labels'].items() if s in system.states}}
        scheme = validate_scheme(validation_payload, system.states)
        if label_scheme.get('scheme_id'):
            scheme['labels'] = dict(label_scheme['labels'])
        scheme.update({k: label_scheme[k] for k in ('scheme_id', 'revision', 'updated_at') if k in label_scheme})
    # Deleted annotations on unselected states do not invalidate a historical scheme.
    if any(scheme['labels'].get(s) not in scheme['categories'] for s in state_ids if s in scheme['labels']):
        raise ValueError('Invalid label assignment in scheme.')
    def progress(message, step):
        if progress_callback:
            progress_callback(message, step, len(state_ids) + len(contrasts) + 3)
    progress('Validating models and label snapshot', 0)
    cluster = store.get_cluster_entry(project_id, system_id, cluster_id)
    cluster_path = store.resolve_path(project_id, system_id, cluster['path'])
    cluster_path = cluster_path / 'cluster.npz' if cluster_path.is_dir() else cluster_path
    with np.load(cluster_path, allow_pickle=False) as archive:
        residue_keys = archive['residue_keys'].astype(str).tolist()
    model_entries = {m['model_id']: m for m in store.list_potts_models(project_id, system_id, cluster_id)}
    pairs = [(c['model_a_id'], c['model_b_id']) for c in contrasts]
    models, model_provenance = {}, []
    for identifier in dict.fromkeys(itertools.chain.from_iterable(pairs)):
        entry = model_entries.get(identifier)
        if not entry or not is_full_model(entry):
            raise ValueError(f'Endpoint must be a full/combined model in this cluster: {identifier}')
        raw = store.resolve_path(project_id, system_id, entry['path']).read_bytes()
        embedded = load_potts_model_metadata(io.BytesIO(raw)) or {}
        for source in (entry, entry.get('params') or {}, embedded, embedded.get('params') or {}):
            if source.get('residue_keys') is not None and list(source['residue_keys']) != residue_keys:
                raise ValueError('Model canonical residue ordering differs from this cluster.')
            if source.get('cluster_id') and source['cluster_id'] != cluster_id:
                raise ValueError('Model categorical basis belongs to another cluster.')
        if not is_full_model({'params': embedded.get('params') or embedded}):
            raise ValueError('Raw delta patches cannot be endpoint models.')
        models[identifier] = load_potts_model(io.BytesIO(raw))
        if len(models[identifier].h) != len(residue_keys):
            raise ValueError('Model residue count differs from cluster.')
        model_provenance.append({'model_id': identifier, 'name': entry.get('name') or identifier, 'sha256': hashlib.sha256(raw).hexdigest(), 'path': entry['path']})
    samples = {s['sample_id']: s for s in store.list_samples(project_id, system_id, cluster_id)}
    selected_samples, labels, sample_provenance, state_provenance = {}, [], [], []
    for i, sid in enumerate(state_ids):
        state = system.states[sid]
        progress(f"Reading static assignment: {state.name}", i + 1)
        candidates = [s for s in samples.values() if s.get('type') == 'md_eval' and s.get('state_id') == sid]
        explicit = (state_samples or {}).get(sid)
        if explicit:
            candidates = [s for s in candidates if s['sample_id'] == explicit]
        if not candidates and assign_missing and not explicit and getattr(state, 'descriptor_file', None) and getattr(state, 'n_frames', 0) == 1:
            from phase.workflows.clustering import evaluate_state_with_models
            generated = evaluate_state_with_models(project_id, system_id, cluster_id, sid, store=store, workers=1)
            store.save_sample_entry(project_id, system_id, cluster_id, generated['sample_id'], generated)
            samples = {s['sample_id']: s for s in store.list_samples(project_id, system_id, cluster_id)}
            candidates = [s for s in samples.values() if s.get('type') == 'md_eval' and s.get('state_id') == sid]
        if len(candidates) != 1:
            raise ValueError(f"{state.name}: select one assigned single-frame MD sample; build descriptors and Assign MD from the System page if missing.")
        entry = candidates[0]
        sample, provenance = _read_sample(store, project_id, system_id, entry)
        if len(sample.labels) != 1 or sample.labels.ndim != 2 or (sample.invalid_mask is not None and sample.invalid_mask.any()):
            raise ValueError(f"{state.name}: benchmark requires exactly one valid assigned frame, not a trajectory. Create a PDB-only state.")
        if sample.frame_state_ids is not None and any(str(v) not in ('', sid) for v in sample.frame_state_ids):
            raise ValueError(f'{state.name}: sample provenance belongs to another state.')
        observed_labels(sample.labels, next(iter(models.values())).K_list())
        labels.append(sample.labels[0]); selected_samples[sid] = entry['sample_id']; sample_provenance.append(provenance)
        pdb = store.resolve_path(project_id, system_id, state.pdb_file) if getattr(state, 'pdb_file', None) else None
        state_provenance.append({'state_id': sid, 'name': state.name, 'pdb_sha256': hashlib.sha256(pdb.read_bytes()).hexdigest() if pdb and pdb.is_file() else None})
    progress('Computing scores, support, fingerprints and statistics', len(state_ids) + 1)
    summary, arrays, deltas = compute_benchmark(models, pairs, np.asarray(labels), state_ids, scheme, **options)
    arrays['residue_keys'] = np.array(residue_keys, str)
    arrays['K_list'] = np.array(next(iter(models.values())).K_list(), int)
    arrays['assigned_labels'] = np.asarray(labels)
    hotspot_items, hotspot_arrays = [], {}
    all_edges = [tuple(e) for e in arrays['edges'].tolist()]
    edge_position = {edge: j for j, edge in enumerate(all_edges)}
    reference_ids = set()
    for i, (contrast, delta) in enumerate(zip(contrasts, deltas)):
        progress(f'Computing hotspots for axis {i + 1}', len(state_ids) + 2 + i)
        node = np.array([np.sqrt(np.mean(h ** 2)) for h in delta.h])
        edge = np.array([np.sqrt(np.mean(delta.J[e] ** 2)) for e in delta.edges])
        if hotspots_mode == 'ensemble':
            groups, used = [], set()
            for side in ('a', 'b'):
                mid = contrast[f'model_{side}_id']
                ids = contrast.get(f'reference_{side}')
                ids = ids if ids is not None else suggested_references(model_entries[mid], samples.values())
                if not ids or len(ids) != len(set(ids)):
                    raise ValueError('Ensemble hotspots require nonempty, nonduplicated assigned MD references for A and B.')
                blocks = []
                for sample_id in ids:
                    entry = samples.get(sample_id)
                    if not entry or entry.get('type') != 'md_eval':
                        raise ValueError('Hotspot references must be assigned MD samples in this cluster.')
                    sample, provenance = _read_sample(store, project_id, system_id, entry)
                    keep = ~np.any(sample.labels == -1, axis=1)
                    if sample.invalid_mask is not None:
                        keep &= ~sample.invalid_mask
                    frames = sample.frame_indices if sample.frame_indices is not None else np.arange(len(keep))
                    source = sample.frame_state_ids if sample.frame_state_ids is not None else np.full(len(keep), entry.get('state_id') or sample_id)
                    identities = [(str(s) or entry.get('state_id') or sample_id, int(f)) for s, f in zip(source[keep], frames[keep])]
                    if used.intersection(identities) or len(set(identities)) != len(identities):
                        raise ValueError('Hotspot reference MD samples overlap in source state/frame.')
                    used.update(identities)
                    blocks.append(block_statistics(delta, sample.labels[keep], frames[keep], source[keep], block_length))
                    reference_ids.add(sample_id)
                    if not any(p['sample_id'] == sample_id for p in sample_provenance):
                        sample_provenance.append(provenance)
                groups.append(blocks)
            stats = separation_statistics(*groups, bootstrap=bootstrap, seed=seed + i)
            node, edge = stats['separation'][:len(delta.h)], stats['separation'][len(delta.h):]
            for key, value in stats.items():
                hotspot_arrays[f'hotspot_{i}_{key}'] = value
        hotspot_arrays[f'hotspot_{i}_node'] = node
        hotspot_arrays[f'hotspot_{i}_edge'] = edge
        hotspot_arrays[f'hotspot_{i}_edge_index'] = np.array([edge_position[e] for e in delta.edges], int)
        nodes = [{'identity': residue_keys[j], 'index': j, 'score': v} for j, v in enumerate(node)]
        edges = [{'identity': f'{residue_keys[r]}--{residue_keys[s]}', 'index': edge_position[r, s], 'score': v} for (r, s), v in zip(delta.edges, edge)]
        nodes.sort(key=lambda row: (-abs(row['score']) if np.isfinite(row['score']) else np.inf, row['index']))
        edges.sort(key=lambda row: (-abs(row['score']) if np.isfinite(row['score']) else np.inf, row['index']))
        hotspot_items.append({'axis': i, 'mode': hotspots_mode, 'nodes': nodes, 'edges': edges})
    arrays.update(hotspot_arrays)
    axes = [{'axis': i, 'model_a_id': a, 'model_b_id': b,
             'name': contrasts[i].get('name') or f"{model_entries[a].get('name') or a} / {model_entries[b].get('name') or b}"} for i, (a, b) in enumerate(pairs)]
    score_rows = []
    for row, sid in enumerate(state_ids):
        for i, axis in enumerate(axes):
            qc = summary['coverage_qc'][i]
            score_rows.append({'state_id': sid, 'state_name': system.states[sid].name, 'label': scheme['labels'].get(sid),
                               'axis': i, 'axis_name': axis['name'], 'raw_score': arrays['raw_scores'][row, i],
                               'common_score': arrays['common_scores'][row, i], 'raw_rank': arrays['raw_ranks'][row, i],
                               'common_rank': arrays['common_ranks'][row, i], 'observed_nodes': qc['node_counts'][row],
                               'observed_edges': qc['edge_counts'][row], 'common_nodes': qc['common_nodes'], 'common_edges': qc['common_edges']})
    summary.update(axes=axes, states=state_provenance, scores=score_rows,
                   distance_matrices={'euclidean': arrays['euclidean_distance'], 'correlation': arrays['correlation_distance']},
                   hotspots=[{**h, 'nodes': h['nodes'][:100], 'edges': h['edges'][:100]} for h in hotspot_items],
                   hotspot_comparisons=hotspot_comparisons(hotspot_items, top_k))
    try:
        software_version = version('phase')
    except PackageNotFoundError:
        software_version = '0.1.0+source'
    try:
        repo = Path(__file__).resolve().parents[2]
        source_revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, stderr=subprocess.DEVNULL, text=True).strip()
        source_dirty = bool(subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=repo, text=True).strip())
    except (OSError, subprocess.SubprocessError):
        source_revision, source_dirty = None, None
    identifier = str(uuid.uuid4())
    metadata = {'analysis_id': identifier, 'analysis_type': ANALYSIS_TYPE, 'schema_version': 1, 'name': name or 'Experimental structure benchmark',
                'created_at': datetime.now(timezone.utc).isoformat(), 'cluster_id': cluster_id, 'state_ids': state_ids,
                'sample_ids': sorted(set(selected_samples.values()) | reference_ids), 'model_ids': list(models),
                'label_scheme_snapshot': scheme, 'models': model_provenance, 'samples': sample_provenance, 'states': state_provenance,
                'software_version': software_version, 'software_versions': {'phase': software_version, 'numpy': np.__version__, 'scipy': version('scipy')}, 'source_revision': source_revision, 'source_dirty': source_dirty,
                'analyzer_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'parameters': {'contrasts': contrasts, 'state_samples': selected_samples, 'hotspots_mode': hotspots_mode,
                               'top_k': top_k, 'block_length': block_length, 'bootstrap': bootstrap, 'seed': seed, 'assign_missing': assign_missing, **options},
                'fingerprint_normalization': {'ddof': 0, 'constant_std_tolerance': 1e-12},
                'hotspot_support': 'Full model terms (parameter) or reference-MD observations (ensemble); not restricted to static common support.',
                'gauge': 'zero_sum', 'support_definition': 'global intersection of observed nodes and union-of-model edges across selected static assignments',
                'warnings': [DISCLAIMER, 'Raw masked totals have different structural support; compare with common support and coverage QC.',
                             'Small experimental sets and related structures limit statistical inference; p/q are exploratory and do not establish biological meaning.',
                             'Legacy models without residue-key metadata assume the selected cluster categorical basis.'],
                'immutable': True}
    summary.update(metadata=metadata, label_scheme=scheme)
    destination = cluster_path.parent / 'analyses' / ANALYSIS_TYPE / identifier
    destination.parent.mkdir(parents=True, exist_ok=True)
    progress('Writing immutable benchmark artifacts', len(state_ids) + len(contrasts) + 2)
    with tempfile.TemporaryDirectory(prefix='.experimental-benchmark-', dir=cluster_path.parent) as temporary:
        staging = Path(temporary) / identifier
        staging.mkdir()
        np.savez_compressed(staging / 'analysis.npz', **arrays)
        _write_json(staging / 'analysis_metadata.json', clean_json(metadata))
        _write_json(staging / 'summary.json', clean_json(summary))
        _write_json(staging / 'statistics.json', clean_json({'statistics': summary['statistics'], 'axis_redundancy': summary['axis_redundancy'], 'axis_algebra': summary['axis_algebra']}))
        _write_json(staging / 'labels.json', scheme)
        _write_json(staging / 'supports.json', {'residue_keys': residue_keys, 'edges': arrays['edges'].tolist(), 'common_support': arrays['common_support'].tolist(), 'axis_support': arrays['axis_support'].tolist()})
        _write_json(staging / 'provenance.json', clean_json(metadata))
        def write_csv(filename, rows):
            with (staging / filename).open('w', newline='') as file:
                if rows:
                    writer = csv.DictWriter(file, fieldnames=list(rows[0])); writer.writeheader()
                    writer.writerows(clean_json(rows))
        write_csv('scores.csv', score_rows)
        write_csv('structure_pairs.csv', summary['discordance'])
        write_csv('hotspots.csv', [{'axis': h['axis'], 'kind': kind, **row} for h in hotspot_items for kind in ('nodes', 'edges') for row in h[kind]])
        for key in ('euclidean_distance', 'correlation_distance', 'rank_spearman_matrix', 'common_score_pearson_matrix'):
            values = arrays[key]
            identities = state_ids if key.endswith('distance') else [a['name'] for a in axes]
            write_csv(f'{key}.csv', [{'identity': identity, **{f'{j}:{identities[j]}': value for j, value in enumerate(row)}} for identity, row in zip(identities, values)])
        staging.replace(destination)
    progress('Complete', len(state_ids) + len(contrasts) + 3)
    return metadata


def benchmark_context(store, project, system_id, cluster):
    system = store.get_system(project, system_id)
    store.get_cluster_entry(project, system_id, cluster)
    samples = store.list_samples(project, system_id, cluster)
    models = [m for m in store.list_potts_models(project, system_id, cluster) if is_full_model(m)]
    return {'states': [{'state_id': sid, 'name': state.name, 'n_frames': state.n_frames,
                        'descriptor_ready': bool(state.descriptor_file),
                        'assigned_samples': [{'sample_id': s['sample_id'], 'name': s.get('name') or s['sample_id']} for s in samples if s.get('type') == 'md_eval' and s.get('state_id') == sid]}
                       for sid, state in sorted(system.states.items(), key=lambda item: item[1].name.casefold())],
            'models': [{'model_id': m['model_id'], 'name': m.get('name') or m['model_id'], 'suggested_references': suggested_references(m, samples)} for m in models],
            'md_samples': [{'sample_id': s['sample_id'], 'name': s.get('name') or s['sample_id']} for s in samples if s.get('type') == 'md_eval'],
            'label_schemes': list_schemes(store, project, system_id)}


