FROM python:3.12-slim AS base

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential libpq-dev curl wget git ca-certificates \
        procps \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt /app/requirements.txt
RUN pip install --upgrade pip && pip install -r requirements.txt

RUN pip install \
    pandas>=2.1.0 \
    scipy>=1.11.0 \
    boto3>=1.34.0 \
    google-cloud-storage>=2.14.0 \
    GEOparse>=2.0.3 \
    ete3>=3.1.3 \
    alphagenome

COPY pipeline/ /app/pipeline/
COPY db/ /app/db/
COPY alembic.ini /app/alembic.ini
COPY scripts/nf_wrappers/ /app/scripts/nf_wrappers/
COPY scripts/pipeline_health_check.py /app/scripts/pipeline_health_check.py

# Copy config last so changes to environment.yml (which is gitignored and not
# tracked by Docker's content-hash cache) always land in a fresh layer.
# This prevents stale esm_chunk_size / threshold values from being baked in.
COPY config/ /app/config/

# Force rebuild: 2026-04-10 step14-final-fixes-phewas-selectivity-depmap
# Force rebuild: 2026-04-10 step14-phewas-otp-page-arg-fix-depmap-unboundlocal-fix
# Force rebuild: 2026-04-10 fix-scoring-depmap-probability-thresholds
# Force rebuild: 2026-04-10 fix-safety-score-always-computed-fix-protac-modality-pr-code
ENV PYTHONPATH=/app
# Force rebuild: 2026-04-16 hyphy-self-sufficient-architecture
