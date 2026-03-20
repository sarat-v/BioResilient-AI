FROM bioresilient/base:latest

ARG IQTREE_VERSION=2.3.6

RUN apt-get update && apt-get install -y --no-install-recommends wget \
    && rm -rf /var/lib/apt/lists/*

RUN wget -qO /tmp/iqtree.tar.gz \
        "https://github.com/iqtree/iqtree2/releases/download/v${IQTREE_VERSION}/iqtree-${IQTREE_VERSION}-Linux-intel.tar.gz" \
    && tar xzf /tmp/iqtree.tar.gz -C /opt/ \
    && ln -s /opt/iqtree-${IQTREE_VERSION}-Linux-intel/bin/iqtree2 /usr/local/bin/iqtree2 \
    && ln -s /usr/local/bin/iqtree2 /usr/local/bin/iqtree \
    && rm /tmp/iqtree.tar.gz

RUN iqtree2 --version
