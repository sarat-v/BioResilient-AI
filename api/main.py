"""BioResilient AI — FastAPI application.

Phase 1 endpoints:
  GET /candidates            Ranked candidate list (filterable by tier)
  GET /candidates/{gene_id}  Full candidate detail with motifs and orthologs
  GET /species               Species registry
  GET /pipeline/status       Pipeline run status
  POST /pipeline/run         Start or resume the pipeline
  GET /pipeline/logs         SSE log stream

Start with:
  uvicorn api.main:app --host 0.0.0.0 --port 8000 --reload
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from api.routes import candidates, scores, species
from api.routes import pipeline as pipeline_routes

_REPO_ROOT = Path(__file__).resolve().parents[1]
_STATIC_DIR = _REPO_ROOT / "api" / "static" / "dist"

app = FastAPI(
    title="BioResilient AI",
    description=(
        "Bioinformatics pipeline API for identifying therapeutic target candidates "
        "from naturally disease-resistant animals. Phase 1: Layers 1 + 2 (sequence "
        "divergence and evolutionary selection)."
    ),
    version="1.0.0-phase1",
    docs_url="/api/docs",
    redoc_url="/api/redoc",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],   # Restrict before external deployment
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(candidates.router, prefix="/candidates", tags=["candidates"])
app.include_router(species.router, prefix="/species", tags=["species"])
app.include_router(scores.router, prefix="/scores", tags=["scores"])
app.include_router(pipeline_routes.router, prefix="/pipeline", tags=["pipeline"])


@app.get("/api/health", tags=["health"])
def health():
    """Database connectivity check."""
    try:
        from db.session import get_engine
        engine = get_engine()
        with engine.connect():
            pass
        return {"status": "ok", "db": "connected"}
    except Exception as exc:
        return {"status": "error", "db": str(exc)}


# ---------------------------------------------------------------------------
# Serve React SPA — must be LAST so it doesn't swallow API routes
# ---------------------------------------------------------------------------

if _STATIC_DIR.exists():
    app.mount("/assets", StaticFiles(directory=str(_STATIC_DIR / "assets")), name="assets")

    @app.get("/{full_path:path}", include_in_schema=False)
    async def spa_fallback(request: Request, full_path: str):
        # API routes are already handled above — this only fires for UI paths
        index = _STATIC_DIR / "index.html"
        if index.exists():
            return FileResponse(str(index))
        return {"service": "BioResilient AI", "ui": "not built yet — run: cd frontend && npm run build"}
else:
    @app.get("/", include_in_schema=False)
    def root():
        return {
            "service": "BioResilient AI",
            "docs": "/api/docs",
            "ui": "not built yet — run: cd frontend && npm run build",
            "status": "running",
        }
