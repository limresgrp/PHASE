"""Compact distance queries and background cache/JS submission."""

from datetime import datetime, timezone
import json
import uuid

from fastapi import APIRouter, Depends, HTTPException
from fastapi.concurrency import run_in_threadpool
import numpy as np
from pydantic import BaseModel, Field

from backend.api.v1.common import get_queue, project_store
from phase.services.project_store import _write_json
from phase.workflows.distances import validate_distance_preset, _ca_atoms
from phase.workflows.distance_analysis import distance_context, distance_pair_table, pair_detail


router = APIRouter(prefix="/projects/{project_id}/systems/{system_id}/metastable/clusters/{cluster_id}/distance")


class Condition(BaseModel):
    state_ids: list[str] | None = None
    metastable_ids: list[str] | None = None
    label: str | None = None


class PairQuery(BaseModel):
    condition: Condition = Field(default_factory=Condition)
    thresholds: list[float] = Field(default_factory=list)
    histogram_bins: int = Field(default=60, ge=1, le=500)


class BuildRequest(BaseModel):
    max_distance_cutoff: float = Field(default=12, gt=0, allow_inf_nan=False)
    periodic: bool = True
    trajectory_overrides: dict[str, str] = Field(default_factory=dict)


class JsRequest(BaseModel):
    condition_a: Condition
    condition_b: Condition
    thresholds: list[float] = Field(default_factory=list)
    name: str | None = None
    min_bin_frames: int = Field(default=1, ge=1)


class PresetRequest(BaseModel):
    name: str = Field(min_length=1, max_length=150)
    thresholds: list[float] = Field(default_factory=list)
    distance_type: str = "ca_ca"


def _path(project_id, system_id, cluster_id):
    entry = project_store.get_cluster_entry(project_id, system_id, cluster_id)
    return project_store.resolve_path(project_id, system_id, entry["path"])


async def _call(function, *args, **kwargs):
    try:
        return await run_in_threadpool(function, *args, **kwargs)
    except FileNotFoundError as exc:
        raise HTTPException(404, str(exc)) from exc
    except (ValueError, KeyError) as exc:
        raise HTTPException(400, str(exc)) from exc


def _context(project_id, system_id, cluster_id):
    return distance_context(_path(project_id, system_id, cluster_id))


@router.get("/cache")
async def status(project_id: str, system_id: str, cluster_id: str):
    cluster, cache, _ = await _call(_context, project_id, system_id, cluster_id)
    system = await _call(project_store.get_system, project_id, system_id)
    return {"metadata": cache.metadata, "distance_available": cache.available,
            "pair_count": int(cache.edge_index.shape[1]),
            "residue_names": (cluster.metadata.get("descriptor_symmetry") or {}).get("residue_resnames", {}),
            "residue_keys": cluster.residue_keys.tolist(),
            "states": [{"state_id": sid, "name": getattr(system.states.get(sid), "name", sid), "frame_count": int(np.count_nonzero(cluster.frame_state_ids == sid))}
                       for sid in dict.fromkeys(cluster.frame_state_ids.tolist())],
            "metastable_ids": list(dict.fromkeys(cluster.frame_metastable_ids.tolist()))}


def _enqueue(project_id, system_id, cluster_id, queue, kind, params):
    _path(project_id, system_id, cluster_id)
    identifier = str(uuid.uuid4())
    job = queue.enqueue("backend.distance_tasks.run_distance_job", args=(identifier,
        {"project_id": project_id, "system_id": system_id, "cluster_id": cluster_id}, kind, params),
        job_id=f"{kind}-{identifier}", job_timeout="24h", result_ttl=86400)
    return {"status": "queued", "job_id": job.id}


@router.post("/cache")
async def build(project_id: str, system_id: str, cluster_id: str, payload: BuildRequest, queue=Depends(get_queue)):
    return await _call(_enqueue, project_id, system_id, cluster_id, queue, "distance_cache", payload.dict())


@router.post("/analyses")
async def submit_js(project_id: str, system_id: str, cluster_id: str, payload: JsRequest, queue=Depends(get_queue)):
    _, cache, _ = await _call(_context, project_id, system_id, cluster_id)
    if not cache.edge_index.shape[1]:
        raise HTTPException(400, "No residue pairs are available. Use Build / refresh cache first to compute candidate pairs from the stored trajectories. If a cache already exists, rebuild it with a larger candidate cutoff.")
    await _call(validate_distance_preset, cache, payload.thresholds)
    return await _call(_enqueue, project_id, system_id, cluster_id, queue, "distance_js", payload.dict(exclude_none=True))


@router.get("/pairs")
async def pairs(project_id: str, system_id: str, cluster_id: str, state_ids: str | None = None,
                metastable_ids: str | None = None, contact_cutoff: float = 6,
                min_contact_fraction: float | None = None,
                sort: str = "range", descending: bool = True, search: str = "", offset: int = 0, limit: int = 50):
    path = await _call(_path, project_id, system_id, cluster_id)
    def csv(value):
        return [s for s in value.split(",") if s] if value is not None else None
    return await _call(distance_pair_table, path, state_ids=csv(state_ids), metastable_ids=csv(metastable_ids),
                       contact_cutoff=contact_cutoff, min_contact_fraction=min_contact_fraction,
                       sort=sort, descending=descending, search=search, offset=offset, limit=limit)


@router.get("/structure-map/{state_id}")
async def structure_map(project_id: str, system_id: str, cluster_id: str, state_id: str):
    def read():
        import MDAnalysis as mda
        cluster, _, _ = _context(project_id, system_id, cluster_id)
        system = project_store.get_system(project_id, system_id)
        state = system.states.get(state_id)
        if state is None or not state.pdb_file:
            raise ValueError("Select a state with a reference structure.")
        universe = mda.Universe(str(project_store.resolve_path(project_id, system_id, state.pdb_file)))
        try:
            atoms = _ca_atoms(universe, state, cluster.residue_keys, cluster.metadata.get("residue_mapping") or {})
            return {"residues": [{"key": str(key), "resid": int(atom.resid), "resname": str(atom.resname),
                                  "chain": str(getattr(atom, "chainID", "")), "position": atom.position.tolist()}
                                 for key, atom in zip(cluster.residue_keys, atoms)]}
        finally:
            universe.trajectory.close()
    return await _call(read)


@router.post("/pairs/{r}/{s}")
async def pair(project_id: str, system_id: str, cluster_id: str, r: int, s: int, payload: PairQuery):
    path = await _call(_path, project_id, system_id, cluster_id)
    return await _call(pair_detail, path, r, s, **payload.dict())


def _preset_dir(project_id, system_id, cluster_id):
    return _path(project_id, system_id, cluster_id).parent / "ui_setups"


@router.get("/presets")
async def presets(project_id: str, system_id: str, cluster_id: str):
    def get():
        directory = _preset_dir(project_id, system_id, cluster_id)
        _, cache, _ = _context(project_id, system_id, cluster_id)
        entries = []
        for path in directory.glob("*.json"):
            try:
                item = json.loads(path.read_text())
                if item.get("setup_type") != "distance_bins":
                    continue
                payload = item["payload"]
                try:
                    validate_distance_preset(cache, payload["thresholds"], payload.get("distance_type", "ca_ca"))
                    compatible, reason = True, None
                except ValueError as exc:
                    compatible, reason = False, str(exc)
                entries.append({**item, "compatible": compatible, "compatibility_message": reason})
            except (ValueError, KeyError):
                continue
        return {"presets": sorted(entries, key=lambda item: item["name"].casefold())}
    return await _call(get)


@router.post("/presets")
async def save_preset(project_id: str, system_id: str, cluster_id: str, payload: PresetRequest):
    def save():
        _, cache, _ = _context(project_id, system_id, cluster_id)
        validate_distance_preset(cache, payload.thresholds, payload.distance_type)
        identifier = str(uuid.uuid4())
        item = {"setup_id": identifier, "name": payload.name.strip(), "setup_type": "distance_bins",
                "page": "distance_analysis", "created_at": datetime.now(timezone.utc).isoformat(),
                "payload": payload.dict(exclude={"name"})}
        if not item["name"]:
            raise ValueError("Preset name cannot be blank.")
        _write_json(_preset_dir(project_id, system_id, cluster_id) / f"{identifier}.json", item)
        return item
    return await _call(save)


@router.delete("/presets/{preset_id}")
async def delete_preset(project_id: str, system_id: str, cluster_id: str, preset_id: str):
    def delete():
        identifier = str(uuid.UUID(preset_id))
        path = _preset_dir(project_id, system_id, cluster_id) / f"{identifier}.json"
        item = json.loads(path.read_text())
        if item.get("setup_type") != "distance_bins":
            raise ValueError("Not a distance preset.")
        path.unlink()
        return {"status": "deleted"}
    return await _call(delete)
