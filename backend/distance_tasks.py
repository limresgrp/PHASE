"""RQ distance jobs using the existing per-system job/result directories."""

from datetime import datetime, timezone
from pathlib import Path
import sys

# RQ appends its working directory after site-packages. Prefer mounted sources
# over the older PHASE package installed in the worker image, as backend.tasks
# does, before importing any phase modules.
REPO_ROOT = str(Path(__file__).resolve().parents[1])
if not sys.path or sys.path[0] != REPO_ROOT:
    sys.path.insert(0, REPO_ROOT)

from rq import get_current_job

from phase.services.project_store import ProjectStore, _write_json
from phase.workflows.distances import build_distance_cache
from phase.workflows.distance_analysis import run_distance_js_analysis


def run_distance_job(job_uuid, dataset_ref, kind, params):
    store = ProjectStore()
    job = get_current_job()
    project_id, system_id, cluster_id = (dataset_ref[key] for key in ("project_id", "system_id", "cluster_id"))
    result_path = store.ensure_results_directories(project_id, system_id)["jobs_dir"] / f"{job_uuid}.json"
    result = {"job_id": job_uuid, "rq_job_id": job.id if job else f"{kind}-{job_uuid}",
              "analysis_type": kind, "status": "started", "created_at": datetime.now(timezone.utc).isoformat(),
              "system_reference": dataset_ref, "params": params, "results": None, "error": None}
    def progress(message, current, total):
        if kind == "distance_cache":
            offset = 45 if message.startswith("Computing") else 0
            percentage = offset + int(45 * current / max(1, total))
        else:
            percentage = int(95 * current / max(1, total))
        if job:
            # Avoid thousands of Redis writes for long trajectories.
            if job.meta.get("progress") == percentage and job.meta.get("status") == message:
                return
            job.meta.update(progress=percentage, status=message)
            job.save_meta()
    _write_json(result_path, result)
    try:
        progress("Initializing distance analysis", 0, 1)
        if kind == "distance_cache":
            output = build_distance_cache(project_id, system_id, cluster_id, store=store, progress_callback=progress, **params)
        elif kind == "distance_js":
            output = run_distance_js_analysis(project_id, system_id, cluster_id, store=store, progress_callback=progress, **params)
        else:
            raise ValueError(f"Unknown distance job kind: {kind}")
        result.update(status="finished", results=output)
        if job:
            job.meta.update(progress=100, status="Complete")
            job.save_meta()
        return result
    except Exception as exc:
        result.update(status="failed", error=str(exc))
        if job:
            job.meta.update(status=str(exc), error=str(exc))
            job.save_meta()
        raise
    finally:
        result["completed_at"] = datetime.now(timezone.utc).isoformat()
        _write_json(result_path, result)
