FROM bioresilient/base:latest

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget cmake g++ libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

ARG HYPHY_VERSION=2.5.62

RUN wget -qO /tmp/hyphy.tar.gz \
        "https://github.com/veg/hyphy/archive/refs/tags/${HYPHY_VERSION}.tar.gz" \
    && tar xzf /tmp/hyphy.tar.gz -C /tmp/ \
    && cd /tmp/hyphy-${HYPHY_VERSION} \
    && cmake -DCMAKE_INSTALL_PREFIX=/usr/local . \
    && make -j$(nproc) HYPHYMP \
    && make install \
    && rm -rf /tmp/hyphy*

RUN ln -sf /usr/local/bin/HYPHYMP /usr/local/bin/hyphy || true

RUN pip install ete3

RUN hyphy --version || HYPHYMP --version
