"""Missing coordinates are not an unassigned/halo cluster (-1)."""
import numpy as np

MISSING_LABEL = -2


def observed_labels(labels, K):
    """Validate labels and return safe indices plus an explicit observation mask."""
    X = np.asarray(labels, dtype=np.int64)
    sizes = np.asarray(K, dtype=np.int64)
    if X.ndim != 2 or X.shape[1] != len(sizes):
        raise ValueError("Labels must be (frames, model residues).")
    missing = X == MISSING_LABEL
    if np.any((X < 0) & ~missing) or np.any(X >= sizes[None, :]):
        raise ValueError("Invalid cluster label; only -2 denotes missing coordinates (-1 is unassigned/halo).")
    return np.where(missing, 0, X), ~missing


def energy_coverage(labels, edges):
    observed = np.asarray(labels) != MISSING_LABEL
    edge_array = np.asarray(edges, dtype=int).reshape((-1, 2))
    edge_observed = observed[:, edge_array[:, 0]] & observed[:, edge_array[:, 1]]
    return {"missing_mask": ~observed, "node_observed_mask": observed,
            "edge_observed_mask": edge_observed,
            "observed_node_counts": observed.sum(axis=1),
            "observed_edge_counts": edge_observed.sum(axis=1)}
