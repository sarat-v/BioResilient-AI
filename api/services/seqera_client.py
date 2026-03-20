"""Seqera Platform REST API client for BioResilient pipeline orchestration.

Two operating modes — pick the one that matches your Seqera setup:

  FULL API MODE  (launch + monitor via Seqera Launchpad)
  -------------------------------------------------------
  Requires SEQERA_API_TOKEN + SEQERA_WORKSPACE_ID + SEQERA_PIPELINE_ID.
  The API server calls Seqera to launch the pipeline; status + logs are
  fetched via the Seqera API.  Best for production with a registered pipeline.

  TOWER MONITOR MODE  (run Nextflow locally, monitor via Seqera)
  --------------------------------------------------------------
  Requires only SEQERA_API_TOKEN (or TOWER_ACCESS_TOKEN — they are the same).
  The API server spawns a Nextflow subprocess with ``-with-tower``.
  Runs appear automatically in your personal Seqera workspace.
  No workspace ID or pre-registered pipeline needed — ideal for a fresh account.

Environment variables:
  SEQERA_API_TOKEN   — Bearer token (alias: TOWER_ACCESS_TOKEN)
  SEQERA_WORKSPACE_ID — Numeric workspace ID (optional; omit → personal workspace)
  SEQERA_PIPELINE_ID  — Pre-registered Launchpad pipeline ID (full API mode only)
  SEQERA_ORG_NAME     — (optional) org slug for constructing watch URLs
  SEQERA_WORKSPACE_NAME — (optional) workspace slug for constructing watch URLs
"""

import logging
import os
import threading
from typing import Optional

import requests

log = logging.getLogger(__name__)

_SEQERA_BASE = "https://api.cloud.seqera.io"

# ---------------------------------------------------------------------------
# Process → Step ID mapping
# ---------------------------------------------------------------------------

# Maps Nextflow process name (last segment of colon-separated path) → pipeline step ID
PROCESS_TO_STEP: dict[str, str] = {
    "validate_environment":       "step1",
    "download_proteomes":         "step2",
    "run_orthofinder":            "step3",
    "load_orthologs":             "step3b",
    "nucleotide_conservation":    "step3c",
    "align_and_divergence":       "step4",
    "domain_and_consequence":     "step4b",
    "esm1v_scoring":              "step4c",
    "variant_direction":          "step4d",
    "build_species_tree":         "step5",
    "phylo_conservation":         "step3d",
    # Scatter steps: individual task instances mark step as "running"
    "run_meme":                   "step6",
    "collect_meme_results":       "step6",       # completion signal for step6
    "run_fel_busted":             "step6b",
    "collect_fel_busted_results": "step6b",      # completion signal for step6b
    "run_relax":                  "step6c",
    "collect_relax_results":      "step6c",      # completion signal for step6c
    "convergence_scoring":        "step7",
    "convergent_aa":              "step7b",
    "expression_analysis":        "step8",
    "bgee_expression":            "step8b",
    "composite_score_phase1":     "step9",
    "alphagenome_regulatory":     "step10b",
    "disease_annotation":         "step11",
    "rare_variants":              "step11b",
    "literature_search":          "step11c",
    "pathway_convergence":        "step11d",
    "druggability":               "step12",
    "p2rank_pockets":             "step12b",
    "gene_therapy":               "step13",
    "safety_screen":              "step14",
    "depmap_gtex":                "step14b",
    "final_rescore":              "step15",
}

# Collect processes are the definitive "step is done" signal for scatter steps
COLLECT_PROCESSES: frozenset[str] = frozenset({
    "collect_meme_results",
    "collect_fel_busted_results",
    "collect_relax_results",
})

# Scatter runner processes completing does NOT mean the step is finished yet
SCATTER_RUNNER_PROCESSES: frozenset[str] = frozenset({
    "run_meme",
    "run_fel_busted",
    "run_relax",
})

# Normalised status priority for merging multiple tasks into one step status
_STATUS_PRIORITY = {"failed": 3, "running": 2, "complete": 1, "pending": 0}


# ---------------------------------------------------------------------------
# Client
# ---------------------------------------------------------------------------

class SeqeraClient:
    """Thin synchronous client for the Seqera Platform REST API."""

    def __init__(
        self,
        token: str,
        workspace_id: str,
        org_name: str = "",
        workspace_name: str = "",
    ) -> None:
        self.workspace_id = workspace_id
        self.org_name = org_name
        self.workspace_name = workspace_name
        self._session = requests.Session()
        self._session.headers.update({
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/json",
        })

    # -------------------------------------------------------------------
    # Core REST methods
    # -------------------------------------------------------------------

    def launch_workflow(self, pipeline_id: str, params: Optional[dict] = None) -> str:
        """Launch a saved Seqera pipeline. Returns the Seqera workflow ID."""
        body = {"launch": {"id": pipeline_id, "params": params or {}}}
        resp = self._session.post(
            f"{_SEQERA_BASE}/workflow/launch",
            params={"workspaceId": self.workspace_id},
            json=body,
            timeout=30,
        )
        resp.raise_for_status()
        return resp.json()["workflowId"]

    def get_workflow(self, workflow_id: str) -> dict:
        """Return the full workflow object from Seqera."""
        resp = self._session.get(
            f"{_SEQERA_BASE}/workflow/{workflow_id}",
            params={"workspaceId": self.workspace_id},
            timeout=15,
        )
        resp.raise_for_status()
        return resp.json().get("workflow", {})

    def get_workflow_tasks(self, workflow_id: str, max_tasks: int = 500) -> list[dict]:
        """Return the list of task records for a workflow run."""
        tasks: list[dict] = []
        offset = 0
        page_size = min(max_tasks, 100)
        while True:
            resp = self._session.get(
                f"{_SEQERA_BASE}/workflow/{workflow_id}/tasks",
                params={
                    "workspaceId": self.workspace_id,
                    "max": page_size,
                    "offset": offset,
                },
                timeout=15,
            )
            resp.raise_for_status()
            data = resp.json()
            page = data.get("tasks", [])
            tasks.extend(page)
            if len(tasks) >= max_tasks or len(page) < page_size:
                break
            offset += page_size
        return tasks

    def cancel_workflow(self, workflow_id: str) -> None:
        """Cancel a running workflow (no-op if already stopped)."""
        resp = self._session.delete(
            f"{_SEQERA_BASE}/workflow/{workflow_id}",
            params={"workspaceId": self.workspace_id},
            timeout=15,
        )
        if resp.status_code not in (200, 204, 409):  # 409 = already stopped/finished
            resp.raise_for_status()

    def get_workflow_log(self, workflow_id: str) -> str:
        """Return the main Nextflow log as plain text. Empty string if not available."""
        resp = self._session.get(
            f"{_SEQERA_BASE}/workflow/{workflow_id}/log",
            params={"workspaceId": self.workspace_id},
            timeout=30,
        )
        if resp.status_code == 404:
            return ""
        resp.raise_for_status()
        return resp.text

    # -------------------------------------------------------------------
    # Higher-level helpers
    # -------------------------------------------------------------------

    def workflow_status(self, workflow_id: str) -> str:
        """Return a normalised status string.

        Returns one of: "running", "complete", "failed", "stopped".
        """
        wf = self.get_workflow(workflow_id)
        raw = (wf.get("status") or "UNKNOWN").upper()
        return {
            "RUNNING":   "running",
            "SUCCEEDED": "complete",
            "FAILED":    "failed",
            "CANCELLED": "stopped",
            "SUBMITTED": "running",
            "UNKNOWN":   "running",
        }.get(raw, "running")

    def tasks_to_step_statuses(self, workflow_id: str) -> dict:
        """Convert Seqera task list to ``{step_id: {status, elapsed_s}}`` dict.

        Suitable for writing into ``pipeline_state.json["steps"]`` so the
        existing Pipeline page UI can display per-step progress without change.
        """
        tasks = self.get_workflow_tasks(workflow_id)

        # Group tasks by process name (strip module path prefix)
        by_process: dict[str, list[dict]] = {}
        for task in tasks:
            raw_name = task.get("process") or ""
            name = raw_name.split(":")[-1]  # e.g. "PHASE1_SEQUENCE:validate_environment" → "validate_environment"
            by_process.setdefault(name, []).append(task)

        step_info: dict[str, dict] = {}

        for process_name, proc_tasks in by_process.items():
            step_id = PROCESS_TO_STEP.get(process_name)
            if not step_id:
                continue

            is_collect = process_name in COLLECT_PROCESSES
            is_scatter = process_name in SCATTER_RUNNER_PROCESSES

            for task in proc_tasks:
                raw_status = (task.get("status") or "").upper()

                if raw_status in ("SUCCEEDED", "CACHED"):
                    if is_scatter:
                        # Scatter runners finishing means work is in progress
                        # but the step isn't done until the collect process succeeds
                        nf_status = "running"
                    else:
                        nf_status = "complete"
                elif raw_status in ("RUNNING", "SUBMITTED"):
                    nf_status = "running"
                elif raw_status == "FAILED":
                    nf_status = "failed"
                else:
                    nf_status = "pending"

                duration_ms = task.get("duration") or 0
                elapsed_s = round(duration_ms / 1000, 1) if duration_ms else None

                existing = step_info.get(step_id)
                if existing is None or _STATUS_PRIORITY[nf_status] > _STATUS_PRIORITY[existing["status"]]:
                    step_info[step_id] = {"status": nf_status, "elapsed_s": elapsed_s}

        return step_info

    def watch_url(self, workflow_id: str) -> str:
        """Return the Seqera Platform web UI URL to watch this workflow run."""
        if self.org_name and self.workspace_name:
            return (
                f"https://cloud.seqera.io/orgs/{self.org_name}"
                f"/workspaces/{self.workspace_name}/watch/{workflow_id}"
            )
        return f"https://cloud.seqera.io/user/workspace/{self.workspace_id}/watch/{workflow_id}"

    # -------------------------------------------------------------------
    # Factory
    # -------------------------------------------------------------------

    @classmethod
    def from_env(cls) -> Optional["SeqeraClient"]:
        """Instantiate from environment variables. Returns ``None`` if not fully configured."""
        token = _get_token()
        workspace_id = os.getenv("SEQERA_WORKSPACE_ID", "")
        if not token or not workspace_id:
            return None
        return cls(
            token=token,
            workspace_id=workspace_id,
            org_name=os.getenv("SEQERA_ORG_NAME", ""),
            workspace_name=os.getenv("SEQERA_WORKSPACE_NAME", ""),
        )


# ---------------------------------------------------------------------------
# Module-level helper
# ---------------------------------------------------------------------------

def is_configured() -> bool:
    """Return ``True`` if the required Seqera env vars are present (full API mode)."""
    return bool(_get_token()) and bool(os.getenv("SEQERA_WORKSPACE_ID"))


def tower_monitor_mode() -> bool:
    """Return ``True`` if we have a token but NOT a full API-mode setup.

    In this mode Nextflow runs as a local subprocess with ``-with-tower`` and
    reports to the user's personal Seqera workspace.  No pre-registered pipeline
    or workspace ID is needed.
    """
    return bool(_get_token()) and not is_configured()


def _get_token() -> str:
    """Return the Seqera / Tower access token from env, accepting both variable names."""
    return os.getenv("SEQERA_API_TOKEN") or os.getenv("TOWER_ACCESS_TOKEN") or ""


def fetch_workspace_id_from_api(token: str) -> Optional[str]:
    """Auto-fetch the numeric ID of the caller's personal Seqera workspace.

    Calls ``GET /user-info`` and returns the first workspace ID found.
    Returns ``None`` on any failure.
    """
    try:
        resp = requests.get(
            f"{_SEQERA_BASE}/user-info",
            headers={"Authorization": f"Bearer {token}"},
            timeout=10,
        )
        resp.raise_for_status()
        data = resp.json()
        wid = (
            data.get("user", {}).get("id")
            or data.get("workspaceId")
        )
        return str(wid) if wid else None
    except Exception as exc:
        log.debug("fetch_workspace_id_from_api failed: %s", exc)
        return None


def start_log_poller(
    workflow_id: str,
    log_path: str,
    stop_event: threading.Event,
    poll_interval_s: float = 10.0,
) -> None:
    """Background thread target: periodically fetch Seqera logs and append to a local file.

    This keeps the existing SSE ``/pipeline/logs`` endpoint working for Seqera runs
    without any frontend changes — the thread writes fetched log content into
    ``pipeline.log``, and the SSE endpoint tails that file as usual.
    """
    client = SeqeraClient.from_env()
    if not client:
        return
    last_size = 0
    while not stop_event.is_set():
        try:
            log_text = client.get_workflow_log(workflow_id)
            if log_text and len(log_text) > last_size:
                new_content = log_text[last_size:]
                with open(log_path, "a", encoding="utf-8") as fh:
                    fh.write(new_content)
                last_size = len(log_text)
        except Exception as exc:
            log.debug("Seqera log fetch error: %s", exc)
        stop_event.wait(poll_interval_s)
