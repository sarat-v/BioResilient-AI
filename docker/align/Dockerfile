ARG BASE_IMAGE=bioresilient/base:latest
FROM ${BASE_IMAGE}

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget ncbi-blast+ \
    && rm -rf /var/lib/apt/lists/*

ARG MAFFT_VERSION=7.520
ARG MINIMAP2_VERSION=2.28
ARG TRIMAL_VERSION=1.5.0

RUN wget -qO /tmp/mafft.deb \
        "https://mafft.cbrc.jp/alignment/software/mafft_${MAFFT_VERSION}-1_amd64.deb" \
    && dpkg -i /tmp/mafft.deb || apt-get install -f -y \
    && rm /tmp/mafft.deb

RUN wget -qO /tmp/minimap2.tar.bz2 \
        "https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2" \
    && tar xjf /tmp/minimap2.tar.bz2 -C /opt/ \
    && ln -s /opt/minimap2-${MINIMAP2_VERSION}_x64-linux/minimap2 /usr/local/bin/minimap2 \
    && rm /tmp/minimap2.tar.bz2

RUN wget -qO /tmp/trimal.tar.gz \
        "https://github.com/inab/trimal/archive/refs/tags/v${TRIMAL_VERSION}.tar.gz" \
    && tar xzf /tmp/trimal.tar.gz -C /tmp/ \
    && cd /tmp/trimal-${TRIMAL_VERSION}/source \
    && make -j$(nproc) \
    && cp trimal /usr/local/bin/ \
    && rm -rf /tmp/trimal*

RUN apt-get update && apt-get install -y --no-install-recommends \
        lastz phast \
    && rm -rf /var/lib/apt/lists/* \
    || echo "LASTZ/PHAST not in apt — will try conda fallback if needed"

RUN mafft --version && minimap2 --version
