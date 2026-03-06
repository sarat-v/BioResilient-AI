"""Loads app configuration from config/environment.yml.

Priority: environment variables > config file values.
"""

import os
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

_CONFIG_PATH = Path(__file__).resolve().parents[1] / "config" / "environment.yml"
_EXAMPLE_PATH = _CONFIG_PATH.with_name("environment.example.yml")


@lru_cache(maxsize=1)
def _load_raw() -> dict[str, Any]:
    path = _CONFIG_PATH if _CONFIG_PATH.exists() else _EXAMPLE_PATH
    with open(path) as f:
        return yaml.safe_load(f)


def get_config() -> dict[str, Any]:
    return _load_raw()


def get_deployment() -> str:
    return os.environ.get("DEPLOYMENT", get_config().get("deployment", "local"))


def get_db_url() -> str:
    cfg = get_config()
    deployment = get_deployment()
    url_template = cfg["database"][deployment]
    # Expand ${RDS_HOST} if present
    return os.path.expandvars(url_template)


def get_storage_root() -> str:
    cfg = get_config()
    deployment = get_deployment()
    return cfg["storage"][deployment]


def get_ncbi_api_key() -> str:
    cfg = get_config()
    return os.environ.get("NCBI_API_KEY", (cfg.get("ncbi") or {}).get("api_key", ""))


def get_ncbi_email() -> str:
    cfg = get_config()
    return os.environ.get("NCBI_EMAIL", (cfg.get("ncbi") or {}).get("email", ""))


def get_tool_config() -> dict[str, Any]:
    return get_config().get("tools", {})


def get_scoring_weights(phase: str = "phase1") -> dict[str, float]:
    import json
    weights_path = Path(__file__).resolve().parents[1] / "config" / "scoring_weights.json"
    with open(weights_path) as f:
        data = json.load(f)
    return data[phase]


def get_tier_thresholds() -> dict[str, float]:
    import json
    weights_path = Path(__file__).resolve().parents[1] / "config" / "scoring_weights.json"
    with open(weights_path) as f:
        data = json.load(f)
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
