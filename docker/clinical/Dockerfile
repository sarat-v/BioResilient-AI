ARG BASE_IMAGE=bioresilient/base:latest

FROM --platform=linux/amd64 mambaorg/micromamba:1.5.8 AS fpocket-builder

USER root
RUN micromamba install -y -n base \
        -c bioconda -c conda-forge \
        fpocket \
    && micromamba clean -afy

FROM ${BASE_IMAGE}

COPY --from=fpocket-builder /opt/conda /opt/conda

ENV PATH="/opt/conda/bin:$PATH"

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget unzip default-jre-headless r-base r-base-dev r-bioc-deseq2 \
    && rm -rf /var/lib/apt/lists/*

ARG P2RANK_VERSION=2.4.2
RUN wget -qO /tmp/p2rank.tar.gz \
        "https://github.com/rdk/p2rank/releases/download/${P2RANK_VERSION}/p2rank_${P2RANK_VERSION}.tar.gz" \
    && tar xzf /tmp/p2rank.tar.gz -C /opt/ \
    && ln -s /opt/p2rank_${P2RANK_VERSION} /opt/p2rank \
    && rm /tmp/p2rank.tar.gz

# Use base image Python (not conda pip on PATH) + quote >= for shell
RUN /usr/local/bin/python3 -m pip install 'rpy2>=3.5.0'

ENV P2RANK_JAR=/opt/p2rank/p2rank.jar

RUN fpocket --help > /dev/null 2>&1 || true
RUN java -jar /opt/p2rank/p2rank.jar --help > /dev/null 2>&1 || true
