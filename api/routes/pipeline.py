"""Pipeline control endpoints — start runs, monitor status, stream logs.

Supports three execution modes:
  Seqera API mode    — SEQERA_API_TOKEN + SEQERA_WORKSPACE_ID + SEQERA_PIPELINE_ID
                       all set: pipelines are launched on AWS Batch via the Seqera
                       Platform API.  Full status/log monitoring via Seqera.
  Tower monitor mode — only SEQERA_API_TOKEN (or TOWER_ACCESS_TOKEN) set, no
                       workspace/pipeline IDs: Nextflow subprocess runs locally with
                       ``-with-tower`` flag.  Run appears in the user's personal
                       Seqera workspace for monitoring.  Ideal for fresh accounts.
  Local mode         — fallback; launches pipeline/orchestrator.py as a subprocess.

The status and log-streaming endpoints work identically in all modes.
"""

import asyncio
import json
import logging
import os
import subprocess
import sys
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from fastapi import APIRouter, BackgroundTasks, HTTPException
from fastapi.responses import StreamingResponse
from pydantic import BaseModel

from api.services.seqera_client import (
    SeqeraClient,
    _get_token,
    is_configured as seqera_is_configured,
    start_log_poller,
    tower_monitor_mode,
)
from db.models import PipelineRun
from db.session import get_session

log = logging.getLogger(__name__)
router = APIRouter()

_REPO_ROOT = Path(__file__).resolve().parents[2]
_STATE_FILE = _REPO_ROOT / "pipeline_state.json"
_LOG_FILE = _REPO_ROOT / "pipeline.log"
_ORCHESTRATOR = _REPO_ROOT / "pipeline" / "orchestrator.py"

# In-memory handle for the running subprocess (local mode only)
_current_process: Optional[subprocess.Popen] = None


# ---------------------------------------------------------------------------
# State helpers
# ---------------------------------------------------------------------------

def _read_state() -> dict:
    try:
        if _STATE_FILE.exists():
            return json.loads(_STATE_FILE.read_text())
    except Exception:
        pass
    return {"status": "idle", "steps": {}}


def _write_state(state: dict) -> None:
    try:
        _STATE_FILE.write_text(json.dumps(state, indent=2))
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Running-state detection
# ---------------------------------------------------------------------------

def _is_running_local() -> bool:
    """Return True if the local orchestrator subprocess is still alive."""
    global _current_process
    if _current_process is not None:
        if _current_process.poll() is not None:
            _current_process = None
            return False
        return True
    try:
        with get_session() as session:
            run = (
                session.query(PipelineRun)
                .filter(PipelineRun.status == "running")
                .order_by(PipelineRun.started_at.desc())
                .first()
            )
            if run is None or run.pid is None:
                return False
            try:
                os.kill(run.pid, 0)
                return True
            except (ProcessLookupError, PermissionError):
                run.status = "failed"
                session.commit()
                return False
    except Exception:
        return False


def _is_running() -> bool:
    """Return True if a pipeline (Seqera or local) is currently running."""
    state = _read_state()
    if state.get("status") != "running":
        return False
    wf_id = state.get("seqera_workflow_id")
    if wf_id and seqera_is_configured():
        client = SeqeraClient.from_env()
        if client:
            try:
                return client.workflow_status(wf_id) == "running"
            except Exception:
                # Assume running if API is unreachable (fail-safe)
                return True
    return _is_running_local()


# ---------------------------------------------------------------------------
# Seqera background helpers
# ---------------------------------------------------------------------------

def _seqera_poll_and_update(workflow_id: str) -> None:
    """Fetch latest Seqera status + task progress and update pipeline_state.json."""
    client = SeqeraClient.from_env()
    if not client:
        return
    try:
        wf_status = client.workflow_status(workflow_id)
        step_statuses = client.tasks_to_step_statuses(workflow_id)
        state = _read_state()
        state["status"] = wf_status
        state["updated_at"] = datetime.now(timezone.utc).isoformat()
        if step_statuses:
            state["steps"] = step_statuses
        if wf_status in ("complete", "failed", "stopped"):
            state.setdefault("completed_at", datetime.now(timezone.utc).isoformat())
        _write_state(state)
    except Exception as exc:
        log.warning("Seqera status poll failed: %s", exc)


# ---------------------------------------------------------------------------
# Request / Response models
# ---------------------------------------------------------------------------

class RunRequest(BaseModel):
    resume_from: str = "step1"
    dry_run: bool = False


class RunResponse(BaseModel):
    run_id: str
    started_at: str
    resume_from: str
    dry_run: bool
    message: str
    seqera_workflow_id: Optional[str] = None
    seqera_run_url: Optional[str] = None


class PipelineStatus(BaseModel):
    status: str
    current_step: Optional[str] = None
    started_at: Optional[str] = None
    completed_at: Optional[str] = None
    failed_at: Optional[str] = None
    elapsed_total_s: Optional[float] = None
    dry_run: bool = False
    steps: dict = {}
    error: Optional[str] = None
    seqera_workflow_id: Optional[str] = None
    seqera_run_url: Optional[str] = None


# ---------------------------------------------------------------------------
# POST /pipeline/run
# ---------------------------------------------------------------------------

@router.post("/run", response_model=RunResponse)
def start_pipeline(req: RunRequest, background_tasks: BackgroundTasks):
    """Start (or resume) the pipeline.

    If Seqera env vars are configured, launches on Seqera Platform / AWS Batch.
    Otherwise falls back to running ``pipeline/orchestrator.py`` as a subprocess.
    """
    global _current_process

    if _is_running():
        raise HTTPException(status_code=409, detail="Pipeline is already running.")

    import uuid
    run_id = str(uuid.uuid4())
    started_at = datetime.now(timezone.utc).isoformat()

    # ------------------------------------------------------------------
    # Seqera mode
    # ------------------------------------------------------------------
    if seqera_is_configured():
        client = SeqeraClient.from_env()
        pipeline_id = os.getenv("SEQERA_PIPELINE_ID", "")
        if not pipeline_id:
            raise HTTPException(
                status_code=503,
                detail=(
                    "SEQERA_PIPELINE_ID is not set. "
                    "Add it to your environment to launch via Seqera Platform."
                ),
            )

        params = {
            "from_step": req.resume_from,
            "dry_run": str(req.dry_run).lower(),
        }
        try:
            workflow_id = client.launch_workflow(pipeline_id, params)
        except Exception as exc:
            raise HTTPException(
                status_code=502,
                detail=f"Seqera launch failed: {exc}",
            ) from exc

        watch_url = client.watch_url(workflow_id)

        # Initialise state file with Seqera workflow reference
        _write_state({
            "status": "running",
            "run_id": run_id,
            "seqera_workflow_id": workflow_id,
            "seqera_run_url": watch_url,
            "started_at": started_at,
            "updated_at": started_at,
            "resume_from": req.resume_from,
            "dry_run": req.dry_run,
            "steps": {},
        })

        # Background thread: fetch logs from Seqera and write to pipeline.log
        stop_ev = threading.Event()
        threading.Thread(
            target=start_log_poller,
            args=(workflow_id, str(_LOG_FILE), stop_ev),
            daemon=True,
        ).start()

        # Record in DB
        try:
            with get_session() as session:
                run = PipelineRun(
                    id=run_id,
                    status="running",
                    started_at=datetime.now(timezone.utc),
                    step_statuses={},
                )
                session.add(run)
                session.commit()
        except Exception:
            pass

        return RunResponse(
            run_id=run_id,
            started_at=started_at,
            resume_from=req.resume_from,
            dry_run=req.dry_run,
            message="Pipeline launched on Seqera Platform (AWS Batch).",
            seqera_workflow_id=workflow_id,
            seqera_run_url=watch_url,
        )

    # ------------------------------------------------------------------
    # Tower monitor mode: run Nextflow subprocess locally with -with-tower
    # Only SEQERA_API_TOKEN / TOWER_ACCESS_TOKEN is set (no workspace/pipeline IDs).
    # The run is visible on cloud.seqera.io → personal workspace.
    # ------------------------------------------------------------------
    if tower_monitor_mode():
        token = _get_token()
        workspace_id = os.getenv("SEQERA_WORKSPACE_ID", "")

        # Build the nextflow command
        nf_script = _REPO_ROOT / "nextflow" / "main.nf"
        profile = os.getenv("NEXTFLOW_PROFILE", "aws")
        db_url = os.getenv("DB_URL", os.getenv("DATABASE_URL", ""))
        s3_bucket = os.getenv("S3_BUCKET", "bioresilient-data")
        work_dir = f"s3://{s3_bucket}/nf-work"

        cmd = [
            "nextflow", "run", str(nf_script),
            "-profile", profile,
            "-with-tower",
            "-resume",
            "--resume_from", req.resume_from,
            "--outdir", f"s3://{s3_bucket}/results",
            "--work_dir", work_dir,
        ]
        if db_url:
            cmd += ["--db_url", db_url]
        if req.dry_run:
            cmd += ["--dry_run", "true"]

        env = os.environ.copy()
        env["TOWER_ACCESS_TOKEN"] = token
        if workspace_id:
            env["TOWER_WORKSPACE_ID"] = workspace_id

        try:
            with get_session() as session:
                run = PipelineRun(
                    id=run_id,
                    status="running",
                    started_at=datetime.now(timezone.utc),
                    step_statuses={},
                )
                session.add(run)
                session.commit()
        except Exception:
            run_id = None  # type: ignore[assignment]

        if run_id:
            env["PIPELINE_RUN_ID"] = run_id

        if req.resume_from == "step1" and not req.dry_run:
            try:
                _LOG_FILE.write_text("")
            except Exception:
                pass

        _current_process = subprocess.Popen(
            cmd,
            cwd=str(_REPO_ROOT),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            env=env,
        )

        if run_id:
            try:
                with get_session() as session:
                    db_run = session.get(PipelineRun, run_id)
                    if db_run:
                        db_run.pid = _current_process.pid
                        session.commit()
            except Exception:
                pass

        seqera_runs_url = (
            "https://cloud.seqera.io/user/runs"
            if not workspace_id
            else f"https://cloud.seqera.io/orgs/{os.getenv('SEQERA_ORG_NAME', '_')}/workspaces/{os.getenv('SEQERA_WORKSPACE_NAME', workspace_id)}/watch"
        )

        _write_state({
            "status": "running",
            "run_id": run_id,
            "seqera_run_url": seqera_runs_url,
            "started_at": started_at,
            "updated_at": started_at,
            "resume_from": req.resume_from,
            "dry_run": req.dry_run,
            "steps": {},
        })

        def _drain_tower():
            with open(_LOG_FILE, "ab") as lf:
                for line in iter(_current_process.stdout.readline, b""):
                    lf.write(line)
                    lf.flush()
            _current_process.stdout.close()

        threading.Thread(target=_drain_tower, daemon=True).start()

        return RunResponse(
            run_id=run_id or started_at,
            started_at=started_at,
            resume_from=req.resume_from,
            dry_run=req.dry_run,
            message="Pipeline started (Nextflow + Tower monitoring). View at cloud.seqera.io.",
            seqera_run_url=seqera_runs_url,
        )

    # ------------------------------------------------------------------
    # Local subprocess mode
    # ------------------------------------------------------------------
    if req.resume_from == "step1" and not req.dry_run:
        try:
            _LOG_FILE.write_text("")
        except Exception:
            pass

    try:
        with get_session() as session:
            run = PipelineRun(
                id=run_id,
                status="running",
                started_at=datetime.now(timezone.utc),
                step_statuses={},
            )
            session.add(run)
            session.commit()
    except Exception:
        run_id = None  # type: ignore[assignment]

    env = os.environ.copy()
    if run_id:
        env["PIPELINE_RUN_ID"] = run_id

    cmd = [sys.executable, str(_ORCHESTRATOR), f"--resume-from={req.resume_from}"]
    if req.dry_run:
        cmd.append("--dry-run")

    _current_process = subprocess.Popen(
        cmd,
        cwd=str(_REPO_ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=env,
    )

    if run_id:
        try:
            with get_session() as session:
                run = session.get(PipelineRun, run_id)
                if run:
                    run.pid = _current_process.pid
                    session.commit()
        except Exception:
            pass

    def _drain():
        with open(_LOG_FILE, "ab") as lf:
            for line in iter(_current_process.stdout.readline, b""):
                lf.write(line)
                lf.flush()
        _current_process.stdout.close()

    threading.Thread(target=_drain, daemon=True).start()

    response_run_id = run_id if run_id else started_at.replace(":", "-").replace(".", "-")
    return RunResponse(
        run_id=response_run_id,
        started_at=started_at,
        resume_from=req.resume_from,
        dry_run=req.dry_run,
        message="Pipeline started." if not req.dry_run else "Dry run started.",
    )


# ---------------------------------------------------------------------------
# POST /pipeline/stop
# ---------------------------------------------------------------------------

@router.post("/stop")
def stop_pipeline():
    """Terminate the running pipeline (Seqera cancel or local process terminate)."""
    global _current_process

    if not _is_running():
        return {"message": "No pipeline is running."}

    state = _read_state()
    wf_id = state.get("seqera_workflow_id")

    if wf_id and seqera_is_configured():
        client = SeqeraClient.from_env()
        if client:
            try:
                client.cancel_workflow(wf_id)
            except Exception as exc:
                log.warning("Seqera cancel failed: %s", exc)
    elif _current_process is not None:
        _current_process.terminate()
        _current_process = None

    run_id = state.get("run_id")
    try:
        with get_session() as session:
            run = session.get(PipelineRun, run_id) if run_id else None
            if run is None:
                run = (
                    session.query(PipelineRun)
                    .filter(PipelineRun.status == "running")
                    .order_by(PipelineRun.started_at.desc())
                    .first()
                )
            if run:
                run.status = "stopped"
                run.finished_at = datetime.now(timezone.utc)
                session.commit()
    except Exception:
        pass

    state["status"] = "stopped"
    _write_state(state)
    return {"message": "Pipeline stopped."}


# ---------------------------------------------------------------------------
# GET /pipeline/status
# ---------------------------------------------------------------------------

@router.get("/status", response_model=PipelineStatus)
def get_status():
    """Return current pipeline state (step progress, timing, Seqera run URL).

    When a Seqera workflow ID is present and the run is active, polls the
    Seqera API synchronously to refresh the cached state before returning.
    """
    state = _read_state()
    wf_id = state.get("seqera_workflow_id")

    if wf_id and seqera_is_configured() and state.get("status") in ("running", "submitted"):
        # Synchronous poll to get fresh data on each status request
        _seqera_poll_and_update(wf_id)
        state = _read_state()
    elif not wf_id:
        # Local mode: sync running flag from live process
        if _is_running_local() and state.get("status") != "running":
            state["status"] = "running"
        elif not _is_running_local() and state.get("status") == "running":
            state["status"] = "failed"

    return PipelineStatus(**{k: state.get(k) for k in PipelineStatus.model_fields})


# ---------------------------------------------------------------------------
# GET /pipeline/logs  (Server-Sent Events)
# ---------------------------------------------------------------------------

@router.get("/logs")
async def stream_logs(tail: int = 200):
    """Stream pipeline log output as Server-Sent Events.

    For local runs: reads directly from ``pipeline.log``.
    For Seqera runs: the background log-poller thread writes Seqera logs into
    ``pipeline.log``, so this endpoint works identically without any changes.
    """

    async def _generator():
        # Seed with existing tail lines
        if _LOG_FILE.exists():
            lines = _LOG_FILE.read_text(errors="replace").splitlines()
            for line in lines[-tail:]:
                yield f"data: {line}\n\n"

        # Stream new lines as they arrive
        last_size = _LOG_FILE.stat().st_size if _LOG_FILE.exists() else 0
        while True:
            await asyncio.sleep(0.5)
            if not _LOG_FILE.exists():
                continue
            current_size = _LOG_FILE.stat().st_size
            if current_size > last_size:
                with open(_LOG_FILE, errors="replace") as f:
                    f.seek(last_size)
                    new_text = f.read()
                last_size = current_size
                for line in new_text.splitlines():
                    if line.strip():
                        yield f"data: {line}\n\n"
            yield ": heartbeat\n\n"

    return StreamingResponse(
        _generator(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "X-Accel-Buffering": "no",
        },
    )
