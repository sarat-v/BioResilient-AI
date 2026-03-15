"""Loads app configuration from config/environment.yml.

Priority: environment variables > config file values.
"""

import json
import logging
import os
from functools import lru_cache
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlencode, urlparse, urlunparse

import yaml

log = logging.getLogger(__name__)

_CONFIG_PATH = Path(__file__).resolve().parents[1] / "config" / "environment.yml"
_EXAMPLE_PATH = _CONFIG_PATH.with_name("environment.example.yml")
_SCORING_WEIGHTS_PATH = Path(__file__).resolve().parents[1] / "config" / "scoring_weights.json"


@lru_cache(maxsize=1)
def _load_raw() -> dict[str, Any]:
    if _CONFIG_PATH.exists():
        path = _CONFIG_PATH
    elif _EXAMPLE_PATH.exists():
        log.warning(
            "config/environment.yml not found — using example config. "
            "Copy environment.example.yml to environment.yml and fill in your values."
        )
        path = _EXAMPLE_PATH
    else:
        raise FileNotFoundError(
            f"No pipeline config found at {_CONFIG_PATH} or {_EXAMPLE_PATH}. "
            "Create config/environment.yml from the example template."
        )
    with open(path) as f:
        return yaml.safe_load(f)


@lru_cache(maxsize=1)
def _load_scoring_weights() -> dict[str, Any]:
    """Load and cache scoring_weights.json once for the process lifetime."""
    if not _SCORING_WEIGHTS_PATH.exists():
        raise FileNotFoundError(
            f"Scoring weights config not found at {_SCORING_WEIGHTS_PATH}. "
            "Ensure config/scoring_weights.json exists."
        )
    with open(_SCORING_WEIGHTS_PATH) as f:
        return json.load(f)


def get_config() -> dict[str, Any]:
    return _load_raw()


def get_deployment() -> str:
    return os.environ.get("DEPLOYMENT", get_config().get("deployment", "local"))


def get_db_url() -> str:
    # Allow full override via DATABASE_URL env var (used in tests and CI)
    if os.environ.get("DATABASE_URL"):
        return os.environ["DATABASE_URL"]
    cfg = get_config()
    deployment = get_deployment()
    db_section = cfg.get("database", {})
    if deployment not in db_section:
        raise KeyError(
            f"Database config for deployment '{deployment}' not found in environment.yml. "
            f"Available deployments: {list(db_section.keys())}"
        )
    url_template = db_section[deployment]
    if not url_template:
        raise ValueError(
            f"Database URL for deployment '{deployment}' is empty in environment.yml."
        )
    # Expand ${RDS_HOST} and similar env vars
    return os.path.expandvars(url_template)


def get_psycopg2_conn():
    """Return a raw psycopg2 connection, correctly handling sslmode from the URL.

    psycopg2 supports ?sslmode=require in the DSN string, but using
    psycopg2.connect(url) can silently hang on SSL handshake with some RDS
    configurations. This helper strips the sslmode param from the URL and
    passes it as an explicit keyword argument to guarantee it is applied.

    The connect_timeout can be overridden via config thresholds.db_connect_timeout.
    """
    import psycopg2

    url = get_db_url()
    parsed = urlparse(url)
    params = parse_qs(parsed.query)

    sslmode = params.pop("sslmode", ["require"])[0]
    clean_query = urlencode({k: v[0] for k, v in params.items()})
    clean_url = urlunparse(parsed._replace(query=clean_query))

    thresholds = get_thresholds()
    connect_timeout = int(thresholds.get("db_connect_timeout", 30))

    return psycopg2.connect(clean_url, sslmode=sslmode, connect_timeout=connect_timeout)


def get_storage_root() -> str:
    cfg = get_config()
    deployment = get_deployment()
    storage_section = cfg.get("storage", {})
    if deployment not in storage_section:
        raise KeyError(
            f"Storage config for deployment '{deployment}' not found in environment.yml. "
            f"Available deployments: {list(storage_section.keys())}"
        )
    return storage_section[deployment]


def get_local_storage_root() -> str:
    """Path to local data root. When storage is S3, returns /tmp/bioresilient so
    tools that need local disk (e.g. OrthoFinder) use the same path as download.py."""
    root = get_storage_root()
    if root.startswith("s3://") or root.startswith("s3:"):
        cfg = get_config()
        return cfg.get("local_tmp_root", "/tmp/bioresilient")
    return root


def get_ncbi_api_key() -> str:
    cfg = get_config()
    return os.environ.get("NCBI_API_KEY", (cfg.get("ncbi") or {}).get("api_key", ""))


def get_ncbi_email() -> str:
    cfg = get_config()
    return os.environ.get("NCBI_EMAIL", (cfg.get("ncbi") or {}).get("email", ""))


def get_anthropic_api_key() -> str:
    """API key for Anthropic Claude (research assistant narrative generation)."""
    cfg = get_config()
    return os.environ.get("ANTHROPIC_API_KEY", (cfg.get("anthropic") or {}).get("api_key", ""))


def get_tool_config() -> dict[str, Any]:
    cfg = get_config()
    tools = cfg.get("tools", {})
    # Support both flat format (legacy) and deployment-split format (local/cloud sub-keys)
    deployment = get_deployment()
    if deployment in tools and isinstance(tools[deployment], dict):
        return tools[deployment]
    # Flat format fallback (no sub-keys)
    return {k: v for k, v in tools.items() if not isinstance(v, dict)}


def get_scoring_weights(phase: str = "phase1") -> dict[str, float]:
    data = _load_scoring_weights()
    if phase not in data:
        raise KeyError(
            f"Scoring phase '{phase}' not found in scoring_weights.json. "
            f"Available phases: {[k for k in data.keys() if k != 'tier_thresholds']}"
        )
    return data[phase]


def get_tier_thresholds() -> dict[str, float]:
    data = _load_scoring_weights()
    if "tier_thresholds" not in data:
        raise KeyError(
            "Key 'tier_thresholds' not found in scoring_weights.json. "
            "Add a tier_thresholds section with tier1 and tier2 keys."
        )
    return data["tier_thresholds"]


def get_thresholds() -> dict[str, Any]:
    return get_config().get("thresholds", {
        "divergence_identity_max": 0.85,
        "divergence_min_species": 2,
        "convergence_min_lineages": 3,
        "expression_log2fc_min": 1.0,
        "expression_padj_max": 0.05,
        "tier1_composite_min": 0.70,
        "tier2_composite_min": 0.40,
    })


def get_api_key() -> str:
    """Expected X-API-Key header value. Empty string = auth disabled."""
    return os.environ.get("API_KEY", (get_config().get("api") or {}).get("key", ""))


def cfg_get(dotpath: str, default: Any = None) -> Any:
    """Access nested config keys with dot notation: cfg_get('ncbi.api_key')."""
    cfg = get_config()
    keys = dotpath.split(".")
    node = cfg
    for k in keys:
        if not isinstance(node, dict):
            return default
        node = node.get(k, default)
        if node is None:
            return default
    return node


def get_s3_bucket() -> str:
    """Return S3 bucket name from storage.cloud URI or S3_BUCKET env var."""
    bucket = os.environ.get("S3_BUCKET", "")
    if bucket:
        return bucket
    cloud_uri = (get_config().get("storage") or {}).get("cloud", "")
    if cloud_uri.startswith("s3://"):
        return cloud_uri.replace("s3://", "").split("/")[0]
    return ""


def sync_to_s3(local_path: "Path | str", s3_key: str) -> bool:
    """Upload a local file to S3. No-op when deployment is not 'cloud'.

    Returns True on success, False on error (never raises).
    """
    import boto3

    if get_deployment() != "cloud":
        return True
    bucket = get_s3_bucket()
    if not bucket:
        return False
    try:
        local_path = Path(local_path)
        if not local_path.exists():
            log.warning("S3 upload skipped — local file does not exist: %s", local_path)
            return False
        boto3.client("s3").upload_file(str(local_path), bucket, s3_key)
        log.info("S3 upload: s3://%s/%s", bucket, s3_key)
        return True
    except Exception as exc:
        log.warning("S3 upload failed for %s → %s: %s", local_path, s3_key, exc)
        return False


def sync_from_s3(s3_key: str, local_path: "Path | str") -> bool:
    """Download a file from S3 to a local path. No-op when deployment is not 'cloud'.

    Returns True on success (or file already exists locally), False on error.
    """
    import boto3

    if get_deployment() != "cloud":
        return True
    bucket = get_s3_bucket()
    if not bucket:
        return False
    try:
        local_path = Path(local_path)
        local_path.parent.mkdir(parents=True, exist_ok=True)
        if local_path.exists():
            return True
        boto3.client("s3").download_file(bucket, s3_key, str(local_path))
        log.info("S3 download: s3://%s/%s → %s", bucket, s3_key, local_path)
        return True
    except Exception as exc:
        log.warning("S3 download failed for %s → %s: %s", s3_key, local_path, exc)
        return False
