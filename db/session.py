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
                # Prevent reporters / ad-hoc queries from hanging indefinitely.
                # 5 min is generous for GROUP-BY on 500K+ row tables.
                "options": "-c statement_timeout=300000"
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
