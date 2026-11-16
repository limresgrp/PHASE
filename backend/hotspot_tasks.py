"""Independent RQ entry point; works with one queue worker (no nested jobs)."""
from datetime import datetime, timezone
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rq import get_current_job
from phase.potts.hotspots import ANALYSIS_TYPE, run_hotspot_analysis
from phase.services.project_store import ProjectStore, _write_json


def run_hotspot_job(job_uuid, dataset_ref, params):
    store = ProjectStore()
    job = get_current_job()
    project, system, cluster = (dataset_ref[k] for k in ("project_id", "system_id", "cluster_id"))
    path = store.ensure_results_directories(project, system)["jobs_dir"] / f"{job_uuid}.json"
    result = {"job_id": job_uuid, "rq_job_id": job.id if job else f"hotspots-{job_uuid}",
              "analysis_type": ANALYSIS_TYPE, "status": "started", "params": params,
              "system_reference": dataset_ref, "created_at": datetime.now(timezone.utc).isoformat()}

    def progress(message, current, total):
        if job:
            job.meta.update(progress=int(95 * current / max(1, total)), status=message)
            job.save_meta()
    _write_json(path, result)
    try:
        result["results"] = run_hotspot_analysis(project, system, cluster, store=store, progress_callback=progress, **params)
        result["status"] = "finished"
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
        _write_json(path, result)
