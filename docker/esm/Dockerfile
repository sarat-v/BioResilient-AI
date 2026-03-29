FROM nvidia/cuda:12.1.1-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

RUN apt-get update && apt-get install -y --no-install-recommends \
        python3.11 python3.11-venv python3.11-dev python3-pip \
        build-essential libpq-dev curl wget ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && ln -sf /usr/bin/python3.11 /usr/bin/python3 \
    && ln -sf /usr/bin/python3.11 /usr/bin/python

WORKDIR /app

COPY requirements.txt /app/requirements.txt
RUN pip install --upgrade pip && pip install -r requirements.txt

RUN pip install \
    torch>=2.1.0 \
    --extra-index-url https://download.pytorch.org/whl/cu121

RUN pip install \
    fair-esm>=2.0.0 \
    pandas>=2.1.0 \
    scipy>=1.11.0 \
    boto3>=1.34.0 \
    google-cloud-storage>=2.14.0

COPY pipeline/ /app/pipeline/
COPY db/ /app/db/
COPY config/ /app/config/
COPY alembic.ini /app/alembic.ini
COPY scripts/nf_wrappers/ /app/scripts/nf_wrappers/

ENV PYTHONPATH=/app
