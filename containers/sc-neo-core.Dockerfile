FROM mambaorg/micromamba:1.5.10

USER root

ENV DEBIAN_FRONTEND=noninteractive
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:$PATH

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        gzip \
        less \
        procps \
        tcsh \
    && rm -rf /var/lib/apt/lists/*

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 \
        samtools=1.19 \
        bcftools=1.19 \
        htslib=1.19 \
        blast=2.15.0 \
        subread \
        pandas=2.2.3 \
        anndata=0.10.8 \
        numpy=1.26.4 \
        h5py=3.12.1 \
        pysam \
        biopython \
    && micromamba clean -a -y

RUN python3 -c "import anndata, h5py, numpy, pandas; print('anndata', anndata.__version__); print('pandas', pandas.__version__); print('numpy', numpy.__version__); print('h5py', h5py.__version__)"

# Licensed/model tools are intentionally not installed here. Bind mount local
# installations and set absolute *_cmd paths in conf/local.config.

WORKDIR /work

CMD ["/bin/bash"]
