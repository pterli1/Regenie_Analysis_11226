# Consolidated Docker image for all regenie WDL tasks + the UKB meta-analysis applet.
# Base: Ubuntu 22.04 (same major as current runtime) with Google Cloud SDK
# Provides: gcloud + gsutil, Python 3.10, fedcat + fedcat_rs, regenie + plink2,
#           the `fedcat collinearity` numerics (scipy + pgenlib), zstd, utilities
#
# Each supported regenie version is copied out of its own hwright0/ubuntu_regenie_gcloud
# image, built from regenie_workflows/aou_regenie_workflows/docker/Dockerfile, whose tag
# is fully determined by its pinned REGENIE_VERSION/CLOUDSDK_VERSION. Each binary is
# installed as /usr/local/bin/regenie-<version>; tasks select one with the config's
# regenie_version, and plain `regenie` points at the default.
# To add a regenie version, follow these steps:
# 1. Update the REGENIE_VERSION (and CLOUDSDK_VERSION if wanted) in
#    regenie_workflows/aou_regenie_workflows/docker/Dockerfile.
# 2. Let the regenie-mkl workflow publish the new tag (older tags stay published).
# 3. Add a REGENIE_IMAGE_<version> ARG, FROM stage and COPY for it below, keeping the others.
# 4. Add the version to REGENIE_VERSIONS in fedcat/core/structures.py.
# Ensure both dockers keep same OS version
# Changed back to this split setup since the cache was too large and building this docker was taking 40 min
#
# Built automatically by .github/workflows/fedcat-meta-docker.yml
# Pushed to Docker Hub as hwright0/fedcat_meta:<fedcat-version>
# AoU pulls via GAR remote-repo: us-central1-docker.pkg.dev/.../hwright0/fedcat_meta:<version>

# ── regenie source images (pinned; see header) ──
ARG REGENIE_IMAGE_4_1_3=hwright0/ubuntu_regenie_gcloud:22.04_4.1.3.gz.mkl_584.0.0
ARG REGENIE_IMAGE_4_1=hwright0/ubuntu_regenie_gcloud:22.04_4.1.gz.mkl_510.0.0

FROM ${REGENIE_IMAGE_4_1_3} AS regenie_4_1_3
FROM ${REGENIE_IMAGE_4_1} AS regenie_4_1

FROM ubuntu:22.04

ARG FEDCAT_WHL=fedcat-*.whl
ARG FEDCAT_RS_WHL=fedcat_rs-*.whl
ARG CLOUDSDK_VERSION=584.0.0

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC \
    PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    ca-certificates \
    gnupg \
    lsb-release \
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    python3-crcmod \
    build-essential \
    pkg-config \
    unzip \
    zip \
    tar \
    gzip \
    bzip2 \
    xz-utils \
    zstd \
    wget \
    jq \
    git \
    coreutils \
    locales \
    && rm -rf /var/lib/apt/lists/*

# Runtime shared libraries required by the regenie + plink2 binaries copied below.
# Mirrors the runtime deps installed in the final stage of the original regenie
# image (Dockerfile_regenie_mkl): libgomp1, gfortran (→ libgfortran5), and the
# htslib compression/crypto libs (z/bz2/lzma/curl/ssl/deflate). Runtime packages
# are used here rather than the -dev variants / full gfortran compiler.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
    libgfortran5 \
    libquadmath0 \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    libcurl4 \
    libssl3 \
    libdeflate0 \
    && rm -rf /var/lib/apt/lists/*

# Install Google Cloud SDK (includes gcloud & gsutil)
RUN curl -s https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor > /usr/share/keyrings/cloud.google.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee /etc/apt/sources.list.d/google-cloud-sdk.list && \
    apt-get update && apt-get install -y --no-install-recommends google-cloud-cli=${CLOUDSDK_VERSION}-0 && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set storage/parallel_composite_upload_enabled False && \
    rm -rf /var/lib/apt/lists/* && \
    gcloud --version

# regenie is statically linked (STATIC=1), so the binary alone is enough.
# `regenie` is the default (DEFAULT_REGENIE_VERSION in fedcat/core/structures.py).
COPY --from=regenie_4_1_3 /usr/local/bin/regenie /usr/local/bin/regenie-4.1.3
COPY --from=regenie_4_1 /usr/local/bin/regenie /usr/local/bin/regenie-4.1
RUN ln -s regenie-4.1.3 /usr/local/bin/regenie

# Download plink2 (linux AVX2 stable build)
ARG PLINK2_URL=https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_latest.zip
RUN curl -fsSL "${PLINK2_URL}" -o /tmp/plink2.zip \
    && unzip /tmp/plink2.zip -d /tmp/plink_extracted/\
    && mv /tmp/plink_extracted/plink2 /usr/local/bin/ \
    && rm -r /tmp/plink2.zip /tmp/plink_extracted/\
    && chmod +x /usr/local/bin/plink2

RUN set -e; for r in /usr/local/bin/regenie*; do echo "$r"; "$r" --version; done; plink2 --version

# pgenlib for the collinearity checker (built from the plink-ng git subdirectory
# — not published to PyPI). Build separately from fedcat for cache.
RUN python3 -m pip install --upgrade pip && \
    pip install --no-cache-dir 'git+https://github.com/chrchang/plink-ng.git#egg=Pgenlib&subdirectory=2.0/Python'

# Install fedcat and its native Rust extension from pre-built wheels.
COPY ${FEDCAT_RS_WHL} /tmp/
COPY ${FEDCAT_WHL} /tmp/
RUN pip install --no-cache-dir /tmp/*.whl && rm -f /tmp/*.whl

RUN echo "fedcat check" && fedcat --help > /dev/null && \
    python3 -c "import fedcat; print('fedcat OK')" && \
    fedcat collinearity --help > /dev/null && \
    python3 -c "import scipy, pgenlib; print('collinearity deps OK')"
