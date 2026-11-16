"""HTTP orchestration of the shared reference/alignment ingestion services."""
from pathlib import Path
import tempfile

from fastapi import APIRouter, File, Form, HTTPException, UploadFile
from fastapi.concurrency import run_in_threadpool
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field

from backend.api.v1.common import project_store, serialize_system, stream_upload
from phase.workflows.reference_states import set_reference, preview_state_alignment, load_alignment, import_pdb_ids

router = APIRouter()
BASE = "/projects/{project_id}/systems/{system_id}"


async def call(function, *args, **kwargs):
    try:
        return await run_in_threadpool(function, *args, **kwargs)
    except FileNotFoundError as exc:
        raise HTTPException(404, str(exc)) from exc
    except (ValueError, KeyError) as exc:
        raise HTTPException(400, str(exc)) from exc


@router.post(BASE + "/reference")
async def reference(project_id: str, system_id: str, pdb: list[UploadFile] = File(...),
                    chain_id: str | None = Form(None), selection: str = Form("protein")):
    if len(pdb) != 1:
        raise HTTPException(400, "Upload exactly one reference PDB.")
    pdb = pdb[0]
    with tempfile.TemporaryDirectory(prefix="phase-reference-") as temp:
        path = Path(temp) / "reference.pdb"
        await stream_upload(pdb, path)
        system = await call(set_reference, project_store, project_id, system_id, path,
                            chain_id=chain_id, selection=selection,
                            provenance={"source": "upload", "original_filename": pdb.filename})
    return serialize_system(system)


@router.get(BASE + "/reference")
async def reference_file(project_id: str, system_id: str):
    system = await call(project_store.get_system, project_id, system_id)
    if not system.reference_pdb_file:
        raise HTTPException(404, "No reference structure is configured.")
    return FileResponse(project_store.resolve_path(project_id, system_id, system.reference_pdb_file), filename="reference.pdb")


@router.post(BASE + "/reference/preview")
async def preview(project_id: str, system_id: str, pdb: UploadFile = File(...), chain_id: str | None = Form(None)):
    with tempfile.TemporaryDirectory(prefix="phase-alignment-") as temp:
        path = Path(temp) / "structure.pdb"
        await stream_upload(pdb, path)
        return await call(preview_state_alignment, project_store, project_id, system_id, path,
                          chain_id=chain_id, provenance={"source": "upload", "original_filename": pdb.filename})


@router.get(BASE + "/states/{state_id}/alignment")
async def alignment(project_id: str, system_id: str, state_id: str):
    system = await call(project_store.get_system, project_id, system_id)
    state = system.states.get(state_id)
    if state is None:
        raise HTTPException(404, "State not found.")
    return {"state_id": state_id, "alignment_status": state.alignment_status or "legacy",
            "alignment": await call(load_alignment, project_store, system, state)}


@router.get(BASE + "/alignments")
async def alignments(project_id: str, system_id: str):
    system = await call(project_store.get_system, project_id, system_id)
    return {"reference_chain_id": system.reference_chain_id,
            "states": [{"state_id": sid, "name": state.name, "alignment_status": state.alignment_status or "legacy",
                        "alignment": await call(load_alignment, project_store, system, state)} for sid, state in system.states.items()]}


class PdbIds(BaseModel):
    pdb_ids: list[str] = Field(min_length=1, max_length=20)
    preview: bool = True
    chain_ids: dict[str, str] = Field(default_factory=dict)
    name_template: str = "PDB_{pdb_id}"
    expected_hashes: dict[str, str] | None = None
    build_descriptors: bool = True


@router.post(BASE + "/states/from-pdb-ids")
async def from_pdb_ids(project_id: str, system_id: str, payload: PdbIds):
    await call(project_store.get_system, project_id, system_id)
    results = await call(import_pdb_ids, project_store, project_id, system_id, payload.pdb_ids,
                         preview=payload.preview, chain_ids=payload.chain_ids, name_template=payload.name_template,
                         expected_hashes=payload.expected_hashes, build_descriptors=payload.build_descriptors)
    return {"results": results, "system": serialize_system(await call(project_store.get_system, project_id, system_id))}
