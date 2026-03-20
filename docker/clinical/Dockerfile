FROM bioresilient/base:latest

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget unzip default-jre-headless r-base \
    && rm -rf /var/lib/apt/lists/*

ARG P2RANK_VERSION=2.4.2

RUN wget -qO /tmp/p2rank.tar.gz \
        "https://github.com/rdk/p2rank/releases/download/${P2RANK_VERSION}/p2rank_${P2RANK_VERSION}.tar.gz" \
    && tar xzf /tmp/p2rank.tar.gz -C /opt/ \
    && ln -s /opt/p2rank_${P2RANK_VERSION} /opt/p2rank \
    && rm /tmp/p2rank.tar.gz

RUN apt-get update && apt-get install -y --no-install-recommends \
        fpocket \
    && rm -rf /var/lib/apt/lists/* \
    || echo "fpocket not in apt — installing from source" \
    && apt-get update && apt-get install -y --no-install-recommends \
        libnetcdf-dev \
    && wget -qO /tmp/fpocket.tar.gz \
        "https://github.com/Discngine/fpocket/archive/refs/tags/4.2.tar.gz" \
    && tar xzf /tmp/fpocket.tar.gz -C /tmp/ \
    && cd /tmp/fpocket-4.2 \
    && make -j$(nproc) && make install \
    && rm -rf /tmp/fpocket*

RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install("DESeq2", ask=FALSE, update=FALSE)'

RUN pip install rpy2>=3.5.0

ENV P2RANK_JAR=/opt/p2rank/p2rank.jar

RUN java -jar /opt/p2rank/p2rank.jar --help > /dev/null 2>&1 || true
