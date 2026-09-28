FROM rust:1-bookworm AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

RUN git clone https://github.com/nf-osi/mafsmith.git

WORKDIR /opt/mafsmith

RUN cargo build --release


FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    ca-certificates \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/mafsmith/target/release/mafsmith /usr/local/bin/mafsmith

RUN chmod +x /usr/local/bin/mafsmith \
    && mkdir -p /data/mafsmith

ENV PATH="/usr/local/bin:${PATH}"
ENV MAFSMITH_HOME="/data/mafsmith"
ENV HOME="/data"

CMD ["/bin/bash"]
