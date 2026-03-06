"""Pipeline control endpoints — start runs, monitor status, stream logs."""

import asyncio
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from fastapi import APIRouter, BackgroundTasks, HTTPException
from fastapi.responses import StreamingResponse
from pydantic import BaseModel

router = APIRouter()

_REPO_ROOT = Path(__file__).resolve().parents[2]
_STATE_FILE = _REPO_ROOT / "pipeline_state.json"
_LOG_FILE = _REPO_ROOT / "pipeline.log"
_ORCHESTRATOR = _REPO_ROOT / "pipeline" / "orchestrator.py"

# In-memory handle for the running subprocess (one at a time)
_current_process: Optional[subprocess.Popen] = None


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_state() -> dict:
    try:
        if _STATE_FILE.exists():
            return json.loads(_STATE_FILE.read_text())
    except Exception:
        pass
    return {"status": "idle", "steps": {}}


def _is_running() -> bool:
    global _current_process
    if _current_process is None:
        return False
    poll = _current_process.poll()
    if poll is not None:
        _current_process = None
        return False
    return True


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


# ---------------------------------------------------------------------------
# POST /pipeline/run
# ---------------------------------------------------------------------------

@router.post("/run", response_model=RunResponse)
def start_pipeline(req: RunRequest, background_tasks: BackgroundTasks):
    """Start (or resume) the pipeline. Only one run allowed at a time."""
    global _current_process

    if _is_running():
        raise HTTPException(status_code=409, detail="Pipeline is already running.")

    # Clear old log on a fresh full run
    if req.resume_from == "step1" and not req.dry_run:
        try:
            _LOG_FILE.write_text("")
        except Exception:
            pass

    cmd = [sys.executable, str(_ORCHESTRATOR), f"--resume-from={req.resume_from}"]
    if req.dry_run:
        cmd.append("--dry-run")

    _current_process = subprocess.Popen(
        cmd,
        cwd=str(_REPO_ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    # Stream stdout to log file in background thread
    def _drain():
        with open(_LOG_FILE, "ab") as lf:
            for line in iter(_current_process.stdout.readline, b""):
                lf.write(line)
                lf.flush()
        _current_process.stdout.close()

    import threading
    threading.Thread(target=_drain, daemon=True).start()

    started_at = datetime.now(timezone.utc).isoformat()
    run_id = started_at.replace(":", "-").replace(".", "-")

    return RunResponse(
        run_id=run_id,
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
    """Terminate the running pipeline process."""
    global _current_process
    if not _is_running():
        return {"message": "No pipeline is running."}
    _current_process.terminate()
    _current_process = None
    state = _read_state()
    state["status"] = "stopped"
    _STATE_FILE.write_text(json.dumps(state, indent=2))
    return {"message": "Pipeline stopped."}


# ---------------------------------------------------------------------------
# GET /pipeline/status
# ---------------------------------------------------------------------------

@router.get("/status", response_model=PipelineStatus)
def get_status():
    """Return current pipeline state (step progress, timing)."""
    state = _read_state()
    # Sync running flag from live process
    if _is_running() and state.get("status") != "running":
        state["status"] = "running"
    elif not _is_running() and state.get("status") == "running":
        state["status"] = "failed"
    return PipelineStatus(**{k: state.get(k) for k in PipelineStatus.model_fields})


# ---------------------------------------------------------------------------
# GET /pipeline/logs  (Server-Sent Events)
# ---------------------------------------------------------------------------

@router.get("/logs")
async def stream_logs(tail: int = 200):
    """Stream pipeline log output as Server-Sent Events.

    Sends the last `tail` lines immediately, then streams new lines as they arrive.
    """

    async def _generator():
        # Seed with existing tail
        if _LOG_FILE.exists():
            lines = _LOG_FILE.read_text(errors="replace").splitlines()
            for line in lines[-tail:]:
                yield f"data: {line}\n\n"

        # Then stream new lines
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
            # Send heartbeat every cycle to keep connection alive
            yield ": heartbeat\n\n"

    return StreamingResponse(
        _generator(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "X-Accel-Buffering": "no",
        },
    )
