ARG BASE_IMAGE=bioresilient/base:latest

# PAML codeml image for step6 (branch-site model A positive selection).
# HyPhy is retained as fallback for aBSREL proxy path only; FEL/BUSTED/RELAX
# are no longer part of the pipeline.
# 2026-04-17: rebuilt with df=1 LRT, NULL HyPhy fields, fix load_selection_scores.
FROM --platform=linux/amd64 mambaorg/micromamba:1.5.8 AS tools

USER root
RUN micromamba install -y -n base \
        -c bioconda -c conda-forge \
        paml \
        hyphy \
    && micromamba clean -afy

FROM ${BASE_IMAGE}

COPY --from=tools /opt/conda /opt/conda

ENV PATH="/opt/conda/bin:$PATH"

RUN command -v codeml && command -v hyphy
RUN hyphy --version 2>&1 | head -1
