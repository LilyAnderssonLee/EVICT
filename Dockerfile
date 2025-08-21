# syntax=docker/dockerfile:1.6
FROM mambaorg/micromamba:1.5.10

COPY --chown=$MAMBA_USER:$MAMBA_USER conda_env.yml /tmp/env.yml
RUN --mount=type=cache,target=/opt/conda/pkgs \
    --mount=type=cache,target=/root/.cache/mamba \
    micromamba create -y -n enterovirus_env -f /tmp/env.yml && micromamba clean -a -y

ENV PATH="/opt/conda/envs/enterovirus_env/bin:${PATH}"
ENV NXF_HOME="/home/${MAMBA_USER}/.nextflow"

# Pre-pull the pipeline (non-fatal if network blocked)
RUN mkdir -p "$NXF_HOME" && \
    micromamba run -n enterovirus_env nextflow -version && \
    micromamba run -n enterovirus_env nextflow pull nf-core/taxprofiler -r 1.2.3 || true

WORKDIR /work

# docker buildx build --platform linux/amd64 -t evtyping:latest . --load