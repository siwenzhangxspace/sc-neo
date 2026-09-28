# syntax=docker/dockerfile:1
FROM frankligy123/altanalyze:0.7.0.1

USER root

# Keep the dependency layer independent from repository scripts and reuse the
# pip download cache across image builds.
RUN --mount=type=cache,target=/root/.cache/pip \
    python -m pip install --disable-pip-version-check lxml

COPY containers/altanalyze-nextflow.sh /usr/local/bin/altanalyze-nextflow
COPY containers/build_altanalyze_groups.py /usr/local/bin/build-altanalyze-groups.py
COPY containers/prepare_altanalyze_runtime.py /usr/local/bin/prepare-altanalyze-runtime.py
COPY containers/prune-altanalyze.py /usr/local/bin/prune-altanalyze.py
COPY containers/run_altanalyze.py /usr/local/bin/run-altanalyze.py
RUN chmod +x \
    /usr/local/bin/altanalyze-nextflow \
    /usr/local/bin/build-altanalyze-groups.py \
    /usr/local/bin/prepare-altanalyze-runtime.py \
    /usr/local/bin/prune-altanalyze.py \
    /usr/local/bin/run-altanalyze.py

# Nextflow supplies the process command; retain the upstream runtime and assets.
ENTRYPOINT []
CMD ["/bin/bash"]
