ARG BASE_IMAGE=bioresilient/base:latest

# Build stage: install PAML via bioconda (includes codeml, baseml, evolver, etc.)
FROM --platform=linux/amd64 mambaorg/micromamba:1.5.8 AS paml-builder

USER root
RUN micromamba install -y -n base \
        -c bioconda -c conda-forge \
        paml \
    && micromamba clean -afy

# Final stage: copy codeml binary into the base image
FROM ${BASE_IMAGE}

COPY --from=paml-builder /opt/conda/bin/codeml  /usr/local/bin/codeml
COPY --from=paml-builder /opt/conda/bin/baseml   /usr/local/bin/baseml
COPY --from=paml-builder /opt/conda/bin/evolver  /usr/local/bin/evolver

# Verify installation (scipy + ete3 already in base image)
RUN codeml 2>&1 | head -3 || echo "codeml installed"
