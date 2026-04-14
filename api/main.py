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

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from starlette.middleware.base import BaseHTTPMiddleware
from starlette.requests import Request
from starlette.responses import JSONResponse

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from api.routes import candidates, research, scores, species, auth as auth_routes
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

# JWT auth: protect all data API routes. Auth routes, health, docs, and SPA paths are exempt.
# Routes that require a valid Bearer token:
_JWT_PROTECTED_PREFIXES = ("/candidates", "/species", "/scores", "/research", "/pipeline")


class JWTMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request: Request, call_next):
        path = request.url.path

        # Only protect data API routes — leave auth, health, docs, and SPA paths open
        if not any(path.startswith(p) for p in _JWT_PROTECTED_PREFIXES):
            return await call_next(request)

        secret = os.environ.get("JWT_SECRET_KEY", "")
        if not secret:
            # JWT not configured — open access (dev / local mode)
            return await call_next(request)

        auth_header = request.headers.get("Authorization", "")
        if not auth_header.startswith("Bearer "):
            return JSONResponse(status_code=401, content={"detail": "Authentication required"})

        token = auth_header[len("Bearer "):]
        try:
            from jose import jwt, JWTError
            jwt.decode(token, secret, algorithms=["HS256"])
        except Exception:
            return JSONResponse(status_code=401, content={"detail": "Invalid or expired token"})

        return await call_next(request)


app.add_middleware(JWTMiddleware)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],   # Restrict before external deployment
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(auth_routes.router, prefix="/auth", tags=["auth"])
app.include_router(candidates.router, prefix="/candidates", tags=["candidates"])
app.include_router(species.router, prefix="/species", tags=["species"])
app.include_router(scores.router, prefix="/scores", tags=["scores"])
app.include_router(research.router, prefix="/research", tags=["research"])
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
