ARG BASE_IMAGE=bioresilient/base:latest
FROM ${BASE_IMAGE}

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget unzip mafft fasttree mcl \
    && rm -rf /var/lib/apt/lists/*

ARG ORTHOFINDER_VERSION=2.5.5
ARG DIAMOND_VERSION=2.1.9

RUN wget -qO /tmp/diamond.tar.gz \
        "https://github.com/bbuchfink/diamond/releases/download/v${DIAMOND_VERSION}/diamond-linux64.tar.gz" \
    && tar xzf /tmp/diamond.tar.gz -C /usr/local/bin/ \
    && chmod +x /usr/local/bin/diamond \
    && rm /tmp/diamond.tar.gz

RUN wget -qO /tmp/orthofinder.tar.gz \
        "https://github.com/davidemms/OrthoFinder/releases/download/${ORTHOFINDER_VERSION}/OrthoFinder.tar.gz" \
    && tar xzf /tmp/orthofinder.tar.gz -C /opt/ \
    && ln -s /opt/OrthoFinder/orthofinder /usr/local/bin/orthofinder \
    && chmod +x /opt/OrthoFinder/orthofinder \
    && rm /tmp/orthofinder.tar.gz

RUN orthofinder --help > /dev/null 2>&1 || true \
    && diamond version
