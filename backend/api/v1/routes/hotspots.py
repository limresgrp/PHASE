"""Hotspot job submission, lightweight choices and coordinate-only mapping."""
from typing import Literal
import uuid

from fastapi import APIRouter, Depends, HTTPException
from fastapi.concurrency import run_in_threadpool
from pydantic import BaseModel, Field
import numpy as np

from backend.api.v1.common import get_queue, project_store
from phase.potts.hotspots import is_full_model, suggested_references

router = APIRouter(prefix="/projects/{project_id}/systems/{system_id}/metastable/clusters/{cluster_id}/hotspots")


class HotspotRequest(BaseModel):
    model_a_id: str
    model_b_id: str
    mode: Literal["parameter", "ensemble"] = "ensemble"
    reference_a: list[str] | None = None
    reference_b: list[str] | None = None
    sample_ids: list[str] = Field(default_factory=list)
    block_length: int = Field(default=50, ge=1, le=10000000)
    bootstrap: int = Field(default=200, ge=0, le=5000)
    seed: int = Field(default=0, ge=0)


async def call(function, *args):
    try:
        return await run_in_threadpool(function, *args)
    except FileNotFoundError as exc:
        raise HTTPException(404, str(exc)) from exc
    except (ValueError, KeyError) as exc:
        raise HTTPException(400, str(exc)) from exc


def context(project, system, cluster):
    project_store.get_cluster_entry(project, system, cluster)
    models = [m for m in project_store.list_potts_models(project, system, cluster) if is_full_model(m)]
    samples = project_store.list_samples(project, system, cluster)
    return {"models": [{"model_id": m["model_id"], "name": m.get("name") or m["model_id"],
                        "suggested_references": suggested_references(m, samples)} for m in models],
            "samples": [{"sample_id": s["sample_id"], "name": s.get("name") or s["sample_id"], "type": s.get("type")} for s in samples]}


@router.get("/context")
async def get_context(project_id: str, system_id: str, cluster_id: str):
    return await call(context, project_id, system_id, cluster_id)


@router.post("/analyses")
async def submit(project_id: str, system_id: str, cluster_id: str, payload: HotspotRequest, queue=Depends(get_queue)):
    choices = await call(context, project_id, system_id, cluster_id)
    ids = {m["model_id"] for m in choices["models"]}
    if payload.model_a_id == payload.model_b_id or not {payload.model_a_id, payload.model_b_id} <= ids:
        raise HTTPException(400, "Choose two different full/combined models in this cluster.")
    if 0 < payload.bootstrap < 50:
        raise HTTPException(400, "Use zero or at least 50 bootstrap replicates.")
    sample_ids = {s["sample_id"] for s in choices["samples"]}
    if not set(payload.sample_ids) <= sample_ids:
        raise HTTPException(400, "Additional preference samples must belong to this cluster.")
    if payload.mode == "ensemble":
        models = {m["model_id"]: m for m in choices["models"]}
        a = payload.reference_a if payload.reference_a is not None else models[payload.model_a_id]["suggested_references"]
        b = payload.reference_b if payload.reference_b is not None else models[payload.model_b_id]["suggested_references"]
        md_ids = {s["sample_id"] for s in choices["samples"] if s.get("type") == "md_eval"}
        if not a or not b or set(a) & set(b) or not set(a + b) <= md_ids:
            raise HTTPException(400, "Select nonempty, disjoint A/B reference groups of assigned MD samples in this cluster.")
    identifier = str(uuid.uuid4())
    job = await run_in_threadpool(queue.enqueue, "backend.hotspot_tasks.run_hotspot_job",
        args=(identifier, {"project_id": project_id, "system_id": system_id, "cluster_id": cluster_id}, payload.dict()),
        job_id=f"hotspots-{identifier}", job_timeout="24h", result_ttl=86400)
    return {"status": "queued", "job_id": job.id}


def structure_mapping(project, system_id, cluster, state_id):
    import MDAnalysis as mda
    from phase.workflows.reference_states import canonical_mapping, selection_mapping
    from phase.io.alignment import ResidueRef
    system = project_store.get_system(project, system_id)
    entry = project_store.get_cluster_entry(project, system_id, cluster)
    path = project_store.resolve_path(project, system_id, entry["path"])
    if path.is_dir():
        path = path / "cluster.npz"
    with np.load(path, allow_pickle=False) as data:
        keys = data["residue_keys"].astype(str)
    if state_id == "_reference":
        if not system.reference_pdb_file:
            raise ValueError("This system has no canonical reference PDB; select a state.")
        topology = project_store.resolve_path(project, system_id, system.reference_pdb_file)
        mapping = selection_mapping(topology, {r["residue_key"]: str(ResidueRef(r["chain_id"], r["resid"], r["insertion_code"]))
                                              for r in system.reference_residue_records})
    else:
        state = system.states.get(state_id)
        if not state or not state.pdb_file:
            raise ValueError("Select a state with a PDB.")
        topology = project_store.resolve_path(project, system_id, state.pdb_file)
        mapping = (canonical_mapping(project_store, system, state) if system.reference_pdb_file else state.residue_mapping) or {}
    universe = mda.Universe(str(topology))
    try:
        residues = []
        for key in keys:
            selection = mapping.get(str(key))
            # Legacy mappings are selections. Refuse ambiguous fallback numbering.
            atoms = universe.select_atoms(f"({selection}) and name CA") if selection else []
            if len(atoms) != 1:
                residues.append(None)
                continue
            atom = atoms[0]
            residues.append({"key": str(key), "resid": int(atom.resid), "resname": str(atom.resname),
                             "insertion_code": str(getattr(atom.residue, "icode", "")).strip(),
                             "chain": str(getattr(atom, "chainID", "")), "position": atom.position.tolist()})
        return {"residues": residues}
    finally:
        universe.trajectory.close()


@router.get("/structure-map/{state_id}")
async def structure_map(project_id: str, system_id: str, cluster_id: str, state_id: str):
    return await call(structure_mapping, project_id, system_id, cluster_id, state_id)
