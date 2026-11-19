"""Experimental benchmark jobs and system-scoped, user-authored labels."""
from pathlib import Path
from typing import Literal
import json
import uuid
import zipfile

from fastapi import APIRouter, Depends, HTTPException, BackgroundTasks
from fastapi.concurrency import run_in_threadpool
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field
import tempfile

from backend.api.v1.common import project_store, get_queue
from phase.services.analysis_annotations import list_schemes, save_scheme, delete_scheme, get_scheme, validate_scheme
from phase.analysis.experimental_benchmark import ANALYSIS_TYPE, benchmark_context

router = APIRouter(prefix='/projects/{project_id}/systems/{system_id}/experimental-benchmarks')


def _dump(payload):
    return payload.model_dump() if hasattr(payload, 'model_dump') else payload.dict()


class SchemeRequest(BaseModel):
    name: str = Field(min_length=1, max_length=200)
    kind: Literal['categorical', 'ordinal'] = 'categorical'
    categories: list[str]
    labels: dict[str, str] = Field(default_factory=dict)


class Contrast(BaseModel):
    model_a_id: str
    model_b_id: str
    name: str | None = None
    reference_a: list[str] | None = None
    reference_b: list[str] | None = None


class BenchmarkRequest(BaseModel):
    cluster_id: str
    state_ids: list[str] = Field(min_length=2)
    contrasts: list[Contrast] = Field(min_length=1)
    scheme_id: str | None = None
    label_scheme: SchemeRequest | None = None
    state_samples: dict[str, str] = Field(default_factory=dict)
    name: str | None = None
    assign_missing: bool = True
    hotspots_mode: Literal['parameter', 'ensemble'] = 'parameter'
    top_k: int = Field(default=20, ge=1)
    block_length: int = Field(default=50, ge=1)
    bootstrap: int = Field(default=200, ge=0, le=5000)
    seed: int = Field(default=0, ge=0)
    fingerprint_mode: Literal['standardized', 'raw'] = 'standardized'
    distance: Literal['euclidean', 'correlation'] = 'euclidean'
    positive_label: str | None = None
    negative_label: str | None = None
    best_oriented: bool = False
    bh: bool = True
    global_threshold: float = Field(default=.15, ge=0, le=1, allow_inf_nan=False)
    local_threshold: float | None = Field(default=None, ge=0, allow_inf_nan=False)


async def call(function, *args, **kwargs):
    try:
        return await run_in_threadpool(function, *args, **kwargs)
    except FileNotFoundError as exc:
        raise HTTPException(404, str(exc)) from exc
    except (ValueError, KeyError) as exc:
        raise HTTPException(400, str(exc)) from exc


@router.get('/context')
async def context(project_id: str, system_id: str, cluster_id: str):
    return await call(benchmark_context, project_store, project_id, system_id, cluster_id)


@router.get('/label-schemes')
async def schemes(project_id: str, system_id: str):
    return {'schemes': await call(list_schemes, project_store, project_id, system_id)}


@router.post('/label-schemes')
async def create_scheme(project_id: str, system_id: str, payload: SchemeRequest):
    return await call(save_scheme, project_store, project_id, system_id, _dump(payload))


@router.put('/label-schemes/{scheme_id}')
async def update_scheme(project_id: str, system_id: str, scheme_id: str, payload: SchemeRequest):
    return await call(save_scheme, project_store, project_id, system_id, _dump(payload), scheme_id)


@router.delete('/label-schemes/{scheme_id}')
async def remove_scheme(project_id: str, system_id: str, scheme_id: str):
    await call(delete_scheme, project_store, project_id, system_id, scheme_id)
    return {'status': 'deleted'}


@router.post('/analyses')
async def submit(project_id: str, system_id: str, payload: BenchmarkRequest, queue=Depends(get_queue)):
    context = await call(benchmark_context, project_store, project_id, system_id, payload.cluster_id)
    if not set(payload.state_ids) <= {s['state_id'] for s in context['states']} or len(set(payload.state_ids)) != len(payload.state_ids):
        raise HTTPException(400, 'Choose distinct existing states.')
    if bool(payload.scheme_id) == bool(payload.label_scheme):
        raise HTTPException(400, 'Choose one label scheme or provide one explicit scheme.')
    models = {m['model_id'] for m in context['models']}
    if any(c.model_a_id == c.model_b_id or not {c.model_a_id, c.model_b_id} <= models for c in payload.contrasts):
        raise HTTPException(400, 'Contrasts require two different full/combined models in this cluster.')
    pairs = [(c.model_a_id, c.model_b_id) for c in payload.contrasts]
    if len(set(pairs)) != len(pairs):
        raise HTTPException(400, 'Choose distinct ordered model contrasts.')
    if 0 < payload.bootstrap < 50:
        raise HTTPException(400, 'Bootstrap replicates must be zero or >=50.')
    identifier = str(uuid.uuid4())
    params = _dump(payload); cluster_id = params.pop('cluster_id')
    # Resolve annotations before enqueue: editing a scheme while the job waits
    # must not silently change the submitted experiment.
    if payload.scheme_id:
        snapshot = await call(get_scheme, project_store, project_id, system_id, payload.scheme_id)
        params.update(scheme_id=None, label_scheme=snapshot)
    else:
        snapshot = await call(validate_scheme, _dump(payload.label_scheme), {s['state_id'] for s in context['states']})
    if (payload.positive_label is None) != (payload.negative_label is None):
        raise HTTPException(400, 'Select both positive and negative AUC categories, or neither.')
    if payload.positive_label is not None and (payload.positive_label == payload.negative_label or not {payload.positive_label, payload.negative_label} <= set(snapshot['categories'])):
        raise HTTPException(400, 'AUC requires two distinct declared categories.')
    job = await run_in_threadpool(queue.enqueue, 'backend.experimental_benchmark_tasks.run_benchmark_job',
        args=(identifier, {'project_id': project_id, 'system_id': system_id, 'cluster_id': cluster_id}, params),
        job_id=f'experimental-benchmark-{identifier}', job_timeout='24h', result_ttl=86400)
    return {'job_id': job.id, 'status': 'queued'}


def result_directory(project, system, cluster, identifier):
    # Validate UUID before filesystem access; never accept arbitrary artifact paths.
    identifier = str(uuid.UUID(identifier))
    entry = project_store.get_cluster_entry(project, system, cluster)
    path = project_store.resolve_path(project, system, entry['path'])
    path = path if path.is_dir() else path.parent
    directory = path / 'analyses' / ANALYSIS_TYPE / identifier
    if not (directory / 'analysis_metadata.json').is_file():
        raise FileNotFoundError('Benchmark result not found.')
    return directory


@router.get('/analyses/{analysis_id}')
async def result(project_id: str, system_id: str, analysis_id: str, cluster_id: str):
    directory = await call(result_directory, project_id, system_id, cluster_id, analysis_id)
    return await call(lambda: json.loads((directory / 'summary.json').read_text()))


@router.get('/analyses/{analysis_id}/download')
async def download(project_id: str, system_id: str, analysis_id: str, cluster_id: str, background_tasks: BackgroundTasks):
    directory = await call(result_directory, project_id, system_id, cluster_id, analysis_id)
    def archive():
        with tempfile.NamedTemporaryFile(suffix='.zip', delete=False) as handle:
            path = Path(handle.name)
        try:
            with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as file:
                for artifact in directory.iterdir():
                    if artifact.is_file():
                        file.write(artifact, artifact.name)
            return path
        except Exception:
            path.unlink(missing_ok=True)
            raise
    path = await call(archive)
    background_tasks.add_task(path.unlink, missing_ok=True)
    return FileResponse(path, media_type='application/zip', filename=f'benchmark-{analysis_id}.zip', background=background_tasks)
