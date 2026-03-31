ARG BASE_IMAGE=bioresilient/base:latest

FROM --platform=linux/amd64 mambaorg/micromamba:1.5.8 AS hyphy-builder

USER root
RUN micromamba install -y -n base \
        -c bioconda -c conda-forge \
        hyphy \
    && micromamba clean -afy

FROM ${BASE_IMAGE}

COPY --from=hyphy-builder /opt/conda /opt/conda

ENV PATH="/opt/conda/bin:$PATH"

RUN pip install ete3

RUN hyphy --version || echo "hyphy version check skipped (cross-compile)"
