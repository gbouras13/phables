#
# phables -- ONE monolithic image for Setonix: ROCm base + conda + every
# per-rule conda env pre-built INSIDE the image. Build it once, convert it to a
# single .sif, and every `phables run` inside it finds its conda envs already
# there. Same strategy as hybracter's own container (miniforge installed into
# the image, envs materialised at build time), NOT the per-rule `container:`
# directive approach -- phables has no --container/--prostt5-container flags
# any more, and its Snakemake rules carry plain `conda:` directives only.
#
# Why this exists: a real, hit-in-production race condition. Every sample's
# `phables run` is its own separate process (one SLURM array task each), and
# Snakemake only serialises conda-env creation WITHIN a single process's DAG.
# Many array tasks all needing the same not-yet-built shared env at once
# corrupt each other's `mamba env create` ("Fatal Python error:
# init_fs_encoding ... no codec search functions registered"). With every env
# already baked into the image there is nothing left to create at runtime, and
# so nothing left to race.
#
# Base is Pawsey's own verified-working ROCm image, as requested: it supplies
# the ROCm userspace matching Setonix's MI250X (gfx90a) kernel driver. NOTE
# that the base's own system-python torch (2.7.1) is NOT what the workflow
# uses -- predict_3di runs inside the prostt5-rocm conda env built below, which
# brings its own torch (2.9.1+rocm6.3, see workflow/envs/prostt5-rocm.yaml).
# The base image's torch is left untouched and simply goes unused.
#
# See ../docs/container.md for how to build the .sif and run it on Setonix.
#
FROM quay.io/pawsey/pytorch:2.7.1-rocm6.3.3

ARG DEBIAN_FRONTEND="noninteractive"

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget \
        bzip2 \
        ca-certificates \
        git \
        procps \
    && rm -rf /var/lib/apt/lists/*

#
# Miniforge inside the image -- hybracter's containers/HPC/Dockerfile does
# exactly this (wget the installer, `bash ... -b -p <prefix>`) rather than
# relying on a conda that may or may not exist on the host. conda-forge is
# miniforge's default channel, so no extra channel config is needed.
#
ENV CONDA_DIR=/opt/miniforge3
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p "${CONDA_DIR}" && \
    rm /tmp/miniforge.sh
ENV PATH="${CONDA_DIR}/bin:${PATH}"

# Snakemake lives in the base conda env alongside phables itself -- it is the
# thing that READS the pre-built per-rule envs, so it must not be one of them.
RUN conda install -y -c conda-forge -c bioconda \
        "snakemake>=8" \
    && conda clean -a -y

#
# phables from the exact source in this build context. CI builds this image per
# commit and tags it by commit SHA, so this is the checked-out tree -- not a
# separately-pinned git ref that could drift from the rules and env files the
# image bakes in.
#
COPY . /opt/phables_src
RUN pip install /opt/phables_src

# Pre-build every per-rule conda env. See the script for the full rationale
# (default --conda-prefix, placeholder databases, one pass per flag combo).
RUN bash /opt/phables_src/container/prebuild_envs.sh

# Fail the build loudly rather than pushing a broken image to quay.io.
RUN bash /opt/phables_src/container/test_image.sh
