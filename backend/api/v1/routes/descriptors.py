from typing import Any, Dict, List, Optional

import numpy as np
from fastapi import APIRouter, HTTPException, Query

from backend.api.v1.common import get_state_or_404, project_store
from phase.features.extraction import DIHEDRAL_KEYS
from phase.io.descriptors import load_descriptor_npz
from phase.io.residue_names import residue_names_for_keys


router = APIRouter()


@router.get(
    "/projects/{project_id}/systems/{system_id}/states/{state_id}/descriptors",
    summary="Preview descriptor angles for a state (for visualization)",
)
async def get_state_descriptors(
    project_id: str,
    system_id: str,
    state_id: str,
    residue_keys: Optional[str] = Query(
        None,
        description="Comma-separated residue keys to include; defaults to all keys for the state.",
    ),
    metastable_ids: Optional[str] = Query(
        None,
        description="Comma-separated metastable IDs to filter frames; defaults to all frames.",
    ),
    cluster_id: Optional[str] = Query(
        None,
        description="ID of a saved cluster to use for sample-based cluster coloring (optional).",
    ),
    cluster_label_mode: str = Query(
        "halo",
        description="Cluster label mode for coloring: 'halo' (default) or 'assigned'.",
    ),
    cluster_variant_id: Optional[str] = Query(
        None,
        description="Cluster variant to use: 'original' (default) or a preview patch id.",
    ),
    max_points: int = Query(
        2000,
        ge=10,
        le=50000,
        description="Maximum number of points returned per residue (down-sampled evenly).",
    ),
):
    """
    Returns a down-sampled set of descriptor dihedrals (in degrees) for the requested state.
    Intended for client-side scatter plotting; not for bulk export.
    """
    try:
        system = project_store.get_system(project_id, system_id)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"System '{system_id}' not found.")

    state_meta = get_state_or_404(system, state_id)
    if not state_meta.descriptor_file:
        raise HTTPException(status_code=404, detail="No descriptors stored for this state.")

    descriptor_path = project_store.resolve_path(project_id, system_id, state_meta.descriptor_file)
    if not descriptor_path.exists():
        raise HTTPException(status_code=404, detail="Descriptor file missing on disk.")

    try:
        feature_dict = load_descriptor_npz(descriptor_path)
    except Exception as exc:  # pragma: no cover
        raise HTTPException(status_code=500, detail=f"Failed to load descriptor file: {exc}") from exc

    def _residue_sort_key(key: str) -> int:
        tok = str(key).split("_")[-1]
        try:
            return int(tok)
        except Exception:
            return 0

    residue_key_order = sorted(
        [k for k in feature_dict.keys() if str(k).startswith("res_")],
        key=_residue_sort_key,
    )
    if not residue_key_order:
        raise HTTPException(status_code=500, detail="Descriptor file contained no residue keys.")

    key_to_col = {k: i for i, k in enumerate(residue_key_order)}
    keys_to_use = list(residue_key_order)
    if residue_keys:
        requested = [key.strip() for key in residue_keys.split(",") if key.strip()]
        keys_to_use = [k for k in keys_to_use if k in requested]
        if not keys_to_use:
            raise HTTPException(status_code=400, detail="No matching residue keys found in descriptor file.")

    angles_payload: Dict[str, Any] = {}
    residue_labels: Dict[str, str] = {}
    sample_stride = 1

    resname_map = residue_names_for_keys(
        keys_to_use, reference_records=system.reference_residue_records or [],
        pdb_path=(project_store.resolve_path(project_id, system_id, state_meta.pdb_file)
                  if state_meta.pdb_file else None),
        residue_mapping=state_meta.residue_mapping,
    )

    # --- Metastable filtering ---
    metastable_filter_ids = []
    if metastable_ids:
        metastable_filter_ids = [mid.strip() for mid in metastable_ids.split(",") if mid.strip()]
    meta_id_to_index = {}
    index_to_meta_id = {}
    state_metastables = [
        m for m in (system.metastable_states or []) if m.get("macro_state_id") == state_id
    ]
    if state_metastables:
        for m in state_metastables:
            mid = m.get("metastable_id")
            if mid is None:
                continue
            meta_id_to_index[mid] = m.get("metastable_index")
            if m.get("metastable_index") is not None:
                index_to_meta_id[m.get("metastable_index")] = mid

    # --- Cluster / Sample labels (directory-driven; no cluster.npz reads) ---
    cluster_legend: List[Dict[str, Any]] = []
    cluster_variants: List[Dict[str, Any]] = []
    selected_cluster_variant = "original"
    state_labels_arr: Optional[np.ndarray] = None
    state_frame_lookup: Optional[Dict[int, int]] = None
    cluster_entry: Optional[Dict[str, Any]] = None
    label_mode = str(cluster_label_mode or "halo").lower()
    if label_mode not in {"halo", "assigned"}:
        raise HTTPException(status_code=400, detail="cluster_label_mode must be 'halo' or 'assigned'.")

    if cluster_id:
        try:
            cluster_entry = project_store.get_cluster_entry(project_id, system_id, cluster_id)
        except FileNotFoundError:
            raise HTTPException(status_code=404, detail="Cluster not found.")

        cluster_variants = [{"id": "original", "label": "Original cluster", "status": "confirmed"}]
        selected_cluster_variant = "original"

        samples = project_store.list_samples(project_id, system_id, cluster_id)
        md_samples = [
            s
            for s in samples
            if isinstance(s, dict)
            and str(s.get("type") or "") == "md_eval"
            and str(s.get("state_id") or "") == str(state_meta.state_id)
        ]
        if md_samples:
            md_samples.sort(key=lambda s: str(s.get("created_at") or ""))
            sample_entry = md_samples[-1]
            paths = sample_entry.get("paths") if isinstance(sample_entry, dict) else None
            sample_rel = paths.get("summary_npz") if isinstance(paths, dict) else None
            sample_rel = sample_rel or sample_entry.get("path")
            if sample_rel:
                sample_path = project_store.resolve_path(project_id, system_id, str(sample_rel))
                if not sample_path.exists():
                    try:
                        cluster_dirs = project_store.ensure_cluster_directories(project_id, system_id, cluster_id)
                        alt = cluster_dirs["cluster_dir"] / str(sample_rel)
                        if alt.exists():
                            sample_path = alt
                    except Exception:
                        pass
                if sample_path.exists():
                    with np.load(sample_path, allow_pickle=True) as sample_npz:
                        if label_mode == "halo":
                            if "labels_halo" in sample_npz:
                                state_labels_arr = np.asarray(sample_npz["labels_halo"], dtype=int)
                            elif "assigned__labels" in sample_npz:
                                state_labels_arr = np.asarray(sample_npz["assigned__labels"], dtype=int)
                        else:
                            if "labels" in sample_npz:
                                state_labels_arr = np.asarray(sample_npz["labels"], dtype=int)
                            elif "assigned__labels_assigned" in sample_npz:
                                state_labels_arr = np.asarray(sample_npz["assigned__labels_assigned"], dtype=int)
                            elif "assigned__labels" in sample_npz:
                                state_labels_arr = np.asarray(sample_npz["assigned__labels"], dtype=int)
                        frame_indices = None
                        if "frame_indices" in sample_npz:
                            frame_indices = np.asarray(sample_npz["frame_indices"], dtype=int)
                        elif "assigned__frame_indices" in sample_npz:
                            frame_indices = np.asarray(sample_npz["assigned__frame_indices"], dtype=int)
                        if frame_indices is not None:
                            state_frame_lookup = {int(fidx): idx for idx, fidx in enumerate(frame_indices)}

    # Shared frame selection (metastable filter + sampling) computed once
    labels_meta = None
    needs_meta_labels = bool(metastable_filter_ids) or bool(state_metastables)
    if needs_meta_labels:
        labels_meta = feature_dict.get("metastable_labels")
        if labels_meta is None and state_meta.metastable_labels_file:
            label_path = project_store.resolve_path(project_id, system_id, state_meta.metastable_labels_file)
            if label_path.exists():
                labels_meta = np.load(label_path)
        if labels_meta is None and metastable_filter_ids:
            raise HTTPException(status_code=400, detail="Metastable labels missing for this state.")

    first_arr = feature_dict[keys_to_use[0]]
    total_frames = first_arr.shape[0] if hasattr(first_arr, "shape") else 0
    n_dihedral_dims = (
        min(int(first_arr.shape[2]), len(DIHEDRAL_KEYS))
        if hasattr(first_arr, "shape") and len(first_arr.shape) >= 3
        else 0
    )
    dihedral_keys = list(DIHEDRAL_KEYS[:n_dihedral_dims]) if n_dihedral_dims > 0 else []
    indices = np.arange(total_frames)
    if metastable_filter_ids:
        selected_idx = {meta_id_to_index.get(mid) for mid in metastable_filter_ids if mid in meta_id_to_index}
        if not selected_idx:
            raise HTTPException(status_code=400, detail="Selected metastable IDs not found on this system.")
        mask = np.isin(labels_meta, list(selected_idx))
        indices = np.where(mask)[0]
        if indices.size == 0:
            raise HTTPException(status_code=400, detail="No frames match selected metastable states for this state.")

    n_frames_filtered = indices.size
    sample_stride = max(1, n_frames_filtered // max_points) if n_frames_filtered > max_points else 1
    sample_indices = indices[::sample_stride]
    n_frames_out = n_frames_filtered

    for key in keys_to_use:
        arr = feature_dict[key]
        if arr.ndim != 3 or arr.shape[2] < 1:
            continue

        sampled = arr[sample_indices, 0, :]
        angles_payload[key] = {}
        n_dims = min(sampled.shape[1], len(dihedral_keys))
        for dim_idx, dim_name in enumerate(dihedral_keys[:n_dims]):
            values = sampled[:, dim_idx] * 180.0 / np.pi
            angles_payload[key][dim_name] = [float(value) if np.isfinite(value) else None for value in values]
        angles_payload[key]["missing_mask"] = (~np.all(np.isfinite(sampled), axis=1)).tolist()

        if state_labels_arr is not None:
            res_idx = key_to_col.get(key)
            if res_idx is not None and res_idx < int(state_labels_arr.shape[1]):
                if state_frame_lookup is not None:
                    rows = np.array([state_frame_lookup.get(int(fidx), -1) for fidx in sample_indices], dtype=int)
                    labels_for_res = np.full(sample_indices.shape[0], -1, dtype=int)
                    valid = rows >= 0
                    if np.any(valid):
                        labels_for_res[valid] = state_labels_arr[rows[valid], int(res_idx)].astype(int)
                else:
                    if sample_indices.size == 0 or sample_indices.max() < state_labels_arr.shape[0]:
                        labels_for_res = state_labels_arr[sample_indices, int(res_idx)].astype(int)
                    else:
                        safe_rows = np.clip(sample_indices, 0, state_labels_arr.shape[0] - 1)
                        labels_for_res = state_labels_arr[safe_rows, int(res_idx)].astype(int)
                angles_payload[key]["cluster_labels"] = labels_for_res.tolist()

        label = key
        if key in resname_map:
            label = f"{key}_{resname_map[key]}"
        residue_labels[key] = label

    if cluster_id:
        residue_cluster_ids = sorted(
            {
                int(v)
                for payload in angles_payload.values()
                for v in (payload.get("cluster_labels") or [])
                if int(v) >= 0
            }
        )
        cluster_legend = [{"id": cid, "label": f"c{cid}"} for cid in residue_cluster_ids]

    if not angles_payload:
        raise HTTPException(status_code=500, detail="Descriptor file contained no usable angle data.")

    halo_payload = {}
    if cluster_id and isinstance(cluster_entry, dict) and str(selected_cluster_variant or "original") == "original":
        cluster_residue_keys = list(residue_key_order)
        n_residues = len(cluster_residue_keys)
        state_label_map = {str(sid): str(state.name or sid) for sid, state in (system.states or {}).items()}
        meta_label_map = {
            str(m.get("metastable_id")): str(m.get("name") or m.get("default_name") or m.get("metastable_id"))
            for m in (system.metastable_states or [])
            if m.get("metastable_id")
        }
        md_samples = [s for s in project_store.list_samples(project_id, system_id, cluster_id) if isinstance(s, dict) and str(s.get("type") or "") == "md_eval"]
        by_condition: Dict[str, Dict[str, Any]] = {}
        for sample in sorted(md_samples, key=lambda s: str(s.get("created_at") or "")):
            sid = sample.get("state_id")
            mid = sample.get("metastable_id")
            if sid:
                cond_id = f"state:{sid}"
                cond_label = state_label_map.get(str(sid), str(sample.get("name") or sid))
                cond_type = "macro"
            elif mid:
                cond_id = f"meta:{mid}"
                cond_label = meta_label_map.get(str(mid), str(sample.get("name") or mid))
                cond_type = "metastable"
            else:
                sample_id = str(sample.get("sample_id") or "")
                if not sample_id:
                    continue
                cond_id = f"sample:{sample_id}"
                cond_label = str(sample.get("name") or sample_id)
                cond_type = "md_eval"
            paths = sample.get("paths") if isinstance(sample, dict) else None
            sample_rel = paths.get("summary_npz") if isinstance(paths, dict) else None
            sample_rel = sample_rel or sample.get("path")
            if not sample_rel:
                continue
            sample_path = project_store.resolve_path(project_id, system_id, str(sample_rel))
            if not sample_path.exists():
                try:
                    cluster_dirs = project_store.ensure_cluster_directories(project_id, system_id, cluster_id)
                    alt = cluster_dirs["cluster_dir"] / str(sample_rel)
                    if alt.exists():
                        sample_path = alt
                except Exception:
                    pass
            if not sample_path.exists():
                continue
            try:
                with np.load(sample_path, allow_pickle=True) as sample_npz:
                    if "labels_halo" in sample_npz:
                        labels_arr = np.asarray(sample_npz["labels_halo"], dtype=np.int32)
                    elif "labels" in sample_npz:
                        labels_arr = np.asarray(sample_npz["labels"], dtype=np.int32)
                    else:
                        continue
            except Exception:
                continue
            if labels_arr.ndim != 2:
                continue
            if n_residues <= 0:
                n_residues = int(labels_arr.shape[1])
                if not cluster_residue_keys:
                    cluster_residue_keys = [f"res_{i}" for i in range(n_residues)]
            if labels_arr.shape[1] != n_residues:
                continue
            halo_rate = (
                np.mean(labels_arr == -1, axis=0).astype(float, copy=False)
                if labels_arr.shape[0] > 0
                else np.full(n_residues, np.nan, dtype=float)
            )
            by_condition[cond_id] = {
                "id": cond_id,
                "label": str(cond_label),
                "type": cond_type,
                "rate": halo_rate,
            }
        if by_condition:
            rows = list(by_condition.values())
            halo_payload["halo_rate_residue_keys"] = cluster_residue_keys
            halo_payload["halo_rate_matrix"] = np.stack([row["rate"] for row in rows], axis=0).tolist()
            halo_payload["halo_rate_condition_ids"] = [row["id"] for row in rows]
            halo_payload["halo_rate_condition_labels"] = [row["label"] for row in rows]
            halo_payload["halo_rate_condition_types"] = [row["type"] for row in rows]

    response = {
        "missing_residue_keys": state_meta.missing_residue_keys,
        "residue_keys": keys_to_use,
        "dihedral_keys": dihedral_keys,
        "residue_mapping": state_meta.residue_mapping or {},
        "residue_labels": residue_labels,
        "n_frames": n_frames_out,
        "sample_stride": sample_stride,
        "angles": angles_payload,
        "cluster_legend": cluster_legend,
        "cluster_variants": cluster_variants,
        "cluster_variant_id": selected_cluster_variant,
        "metastable_labels": labels_meta[sample_indices].astype(int).tolist() if labels_meta is not None else [],
        "metastable_legend": [
            {
                "id": m.get("metastable_id"),
                "index": m.get("metastable_index"),
                "label": m.get("name") or m.get("default_name") or m.get("metastable_id"),
            }
            for m in state_metastables
            if m.get("metastable_index") is not None
        ],
        "metastable_filter_applied": bool(metastable_filter_ids),
        **halo_payload,
    }
    return response
