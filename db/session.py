"""SQLAlchemy session factory."""

from contextlib import contextmanager
from typing import Generator

from sqlalchemy import create_engine
from sqlalchemy.orm import Session, sessionmaker

from pipeline.config import get_db_url

_engine = None
_SessionLocal = None
_engine_url = None


def _get_engine():
    global _engine, _engine_url, _SessionLocal
    current_url = get_db_url()
    # Recreate engine if URL has changed (e.g. DATABASE_URL env var set after import)
    if _engine is None or _engine_url != current_url:
        _engine = create_engine(
            current_url,
            pool_pre_ping=True,
            connect_args={
                # 15 min — reporter queries do ORDER BY on 1.69M-row tables joined to
                # ortholog+gene; without a covering index these can take 8-12 min.
                # The consequence_score index (ix_divergent_motif_consequence_score)
                # should bring this under 1s once built, but the higher timeout is a
                # safety net for future large-table reporters.
                "options": "-c statement_timeout=900000"
            },
        )
        _engine_url = current_url
        _SessionLocal = None  # reset session factory too
    return _engine


def _get_session_factory():
    global _SessionLocal
    if _SessionLocal is None:
        _SessionLocal = sessionmaker(bind=_get_engine(), expire_on_commit=False)
    return _SessionLocal


@contextmanager
def get_session() -> Generator[Session, None, None]:
    """Context manager that yields a database session and handles commit/rollback."""
    factory = _get_session_factory()
    session: Session = factory()
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()


def get_engine():
    return _get_engine()
