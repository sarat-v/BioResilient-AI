"""API route tests using FastAPI TestClient.

Uses in-memory SQLite (DATABASE_URL) so no PostgreSQL is required.
Run with: pytest tests/test_api_routes.py -v
"""

import asyncio
import os
import uuid

import pytest
from httpx import ASGITransport, AsyncClient

# Use a file DB so all connections (test + app) share the same data (in-memory SQLite is per-connection)
_test_db_path = os.path.join(os.path.dirname(__file__), "test_api_routes.db")
os.environ["DATABASE_URL"] = f"sqlite:///{_test_db_path}"

from db.models import Base, CandidateScore, Gene, EvolutionScore
from db.session import get_engine


class _SyncASGIClient:
    """Sync wrapper around httpx AsyncClient + ASGITransport for testing."""

    def __init__(self, app):
        self._app = app
        self._base = "http://testserver"

    def _run(self, coro):
        return asyncio.run(coro)

    async def _request(self, method, path, **kwargs):
        async with AsyncClient(
            transport=ASGITransport(app=self._app),
            base_url=self._base,
        ) as client:
            return await client.request(method, path, **kwargs)

    def get(self, path, params=None, **kwargs):
        return self._run(self._request("GET", path, params=params, **kwargs))

    def post(self, path, json=None, **kwargs):
        return self._run(self._request("POST", path, json=json, **kwargs))


@pytest.fixture(scope="module")
def create_tables():
    """Create all tables in the test DB (module scope, once)."""
    # Remove existing file so SQLite gets a fresh schema with all current columns
    if _test_db_path and os.path.isfile(_test_db_path):
        try:
            os.remove(_test_db_path)
        except Exception:
            pass
    engine = get_engine()
    Base.metadata.create_all(engine)


@pytest.fixture(scope="module")
def app_client(create_tables):
    """HTTP client for the FastAPI app with tables created."""
    from api.main import app
    return _SyncASGIClient(app)


@pytest.fixture
def seeded_client(app_client, create_tables):
    """Client with one gene and one candidate score seeded (for detail tests)."""
    from db.session import get_session
    gene_id = str(uuid.uuid4())
    human_gene_id = f"GENE_TEST_{uuid.uuid4().hex[:8]}"
    with get_session() as session:
        session.add(Gene(
            id=gene_id,
            human_gene_id=human_gene_id,
            gene_symbol="TP53",
            human_protein="P04637",
        ))
        session.add(CandidateScore(
            gene_id=gene_id,
            trait_id="",
            composite_score=0.75,
            tier="Tier1",
            convergence_score=0.8,
            selection_score=0.7,
            expression_score=0.6,
            disease_score=0.0,
            druggability_score=0.0,
            safety_score=1.0,
            regulatory_score=0.0,
        ))
        session.add(EvolutionScore(
            gene_id=gene_id,
            dnds_ratio=1.5,
            dnds_pvalue=0.01,
            selection_model="test",
            convergence_count=3,
        ))
        session.commit()
    yield app_client, gene_id


class TestCandidatesAPI:
    """GET /candidates and GET /candidates/{gene_id}."""

    def test_get_candidates_returns_200_and_list(self, app_client):
        r = app_client.get("/candidates", params={"limit": 10})
        assert r.status_code == 200
        data = r.json()
        assert isinstance(data, list)

    def test_get_candidates_accepts_tier_and_species_params(self, app_client):
        r = app_client.get("/candidates", params={"tier": "Tier1", "species_id": "naked_mole_rat", "limit": 5})
        assert r.status_code == 200
        assert isinstance(r.json(), list)

    def test_get_candidate_404_when_not_found(self, app_client):
        r = app_client.get("/candidates/nonexistent-gene-id")
        assert r.status_code == 404

    def test_get_candidate_200_with_detail(self, seeded_client):
        client, gene_id = seeded_client
        r = client.get(f"/candidates/{gene_id}")
        assert r.status_code == 200
        data = r.json()
        assert data["gene_id"] == gene_id
        assert data["gene_symbol"] == "TP53"
        assert "evolution" in data
        assert "orthologs" in data
        assert data.get("disease") is None or isinstance(data["disease"], dict)
        assert data.get("drug_target") is None or isinstance(data["drug_target"], dict)


class TestScoresAPI:
    """GET /scores/{gene_id}."""

    def test_get_scores_404_when_not_found(self, app_client):
        r = app_client.get("/scores/nonexistent-gene-id")
        assert r.status_code == 404

    def test_get_scores_200_with_breakdown(self, seeded_client):
        client, gene_id = seeded_client
        r = client.get(f"/scores/{gene_id}")
        assert r.status_code == 200
        data = r.json()
        assert data["gene_id"] == gene_id
        assert "sub_scores" in data
        assert "regulatory" in data["sub_scores"]
        assert "evolution" in data


class TestPipelineAPI:
    """POST /pipeline/run and /pipeline/status."""

    def test_pipeline_status_returns_200(self, app_client):
        r = app_client.get("/pipeline/status")
        assert r.status_code == 200
        data = r.json()
        assert "status" in data
        assert data["status"] in ("idle", "running", "complete", "failed", "stopped")

    def test_pipeline_run_returns_200_with_run_id(self, app_client):
        r = app_client.post("/pipeline/run", json={"resume_from": "step1", "dry_run": True})
        assert r.status_code == 200
        data = r.json()
        assert "run_id" in data
        assert "started_at" in data
        assert data.get("dry_run") is True


class TestHealthAPI:
    """GET /api/health."""

    def test_health_returns_200(self, app_client):
        r = app_client.get("/api/health")
        assert r.status_code == 200
        data = r.json()
        assert "status" in data
