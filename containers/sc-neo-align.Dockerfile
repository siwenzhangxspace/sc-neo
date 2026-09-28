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
    && rm -rf /var/lib/apt/lists/*

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 \
        openjdk=17 \
        fastqc=0.12.1 \
        fastp=0.23.4 \
        star=2.7.11b \
        samtools=1.19 \
    && micromamba clean -a -y

WORKDIR /work

CMD ["/bin/bash"]
