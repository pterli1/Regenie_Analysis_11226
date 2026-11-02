# regenie v4.1.3 (MKL) + plink2 + gcloud on Ubuntu 22.04
# Maintained by Harry Wright h.wright4@exeter.ac.uk
# Built automatically by .github/workflows/regenie-mkl-docker.yml
# Pushed to Docker Hub as hwright0/ubuntu_regenie_gcloud:<tag>
#
# The tag is 22.04_<REGENIE_VERSION>.gz.mkl_<CLOUDSDK_VERSION>, i.e. fully
# determined by the two ARGs below. Both are pinned so the tag is knowable before
# the build and downstream images can pin it in the same commit — bump them here
# and the tag follows.
#
# Manual build:
#   docker build -t hwright0/ubuntu_regenie_gcloud:22.04_4.1.3.gz.mkl_584.0.0 .

ARG REGENIE_VERSION=4.1.3
ARG CLOUDSDK_VERSION=584.0.0

FROM ubuntu:22.04 AS builder1

ARG BOOST_IO=1
ARG STATIC=1
ARG CMAKE_VERSION_MAJOR=3.13
ARG CMAKE_VERSION_MINOR=0
ARG HTSLIB_VERSION=1.18
ARG BGEN_VERSION=1.1.7
ARG REGENIE_VERSION

WORKDIR /src

ADD http://cmake.org/files/v${CMAKE_VERSION_MAJOR}/cmake-${CMAKE_VERSION_MAJOR}.${CMAKE_VERSION_MINOR}-Linux-x86_64.sh cmake_install.sh
ADD https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2 htslib-$HTSLIB_VERSION.tar.bz2
ADD https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB intel_key.PUB
ADD https://github.com/rgcgithub/regenie/archive/refs/tags/v$REGENIE_VERSION.tar.gz v$REGENIE_VERSION.tar.gz

# install BGEN and HTSlib libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc-9 \
        g++-9 \
        gfortran-9 \
        libz-dev \
        bzip2 \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        curl \
        make \
        gnupg \
        python3 \
        zlib1g-dev \
        apt-transport-https \
        ca-certificates \
        libboost-iostreams-dev \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 60 --slave /usr/bin/g++ g++ /usr/bin/g++-9 \
    && update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-9 70 \
    && update-alternatives --install /usr/bin/f77 f77 /usr/bin/gfortran-9 70 \
    && update-alternatives --install /usr/bin/cc cc /usr/bin/gcc-9 70 \
    && tar -xf htslib-$HTSLIB_VERSION.tar.bz2 \
    && cd htslib-$HTSLIB_VERSION/ \
    && ./configure \
    && make \
    && make install \
    && cd .. \
    && sh cmake_install.sh --prefix=/usr/local --skip-license --exclude-subdir \
    && rm cmake_install.sh \
    && curl -fsSL --retry 5 --retry-delay 10 --retry-all-errors \
            -o bgen.tgz "https://enkre.net/cgi-bin/code/bgen/tarball/release/v$BGEN_VERSION" \
    && gzip -t bgen.tgz \
    && tar -xzf bgen.tgz \
    && rm bgen.tgz \
    && cd v$BGEN_VERSION \
    && python3 waf configure \
    && python3 waf \
    && cd .. \
    && tar -xzf v$REGENIE_VERSION.tar.gz \
    && sh -c 'cat intel_key.PUB | gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null' \
    && sh -c 'echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list' \
    && apt-get update \
    && apt-get install intel-oneapi-mkl-devel -y --no-install-recommends \
    && . /opt/intel/oneapi/setvars.sh \
    && echo "MKL_THREADING_LAYER=GNU" >> /etc/environment \
    && cd /src/regenie-$REGENIE_VERSION \
    && BGEN_PATH=/src/v$BGEN_VERSION HAS_BOOST_IOSTREAM=$BOOST_IO MKLROOT=$MKLROOT HTSLIB_PATH=/usr/local/lib/ STATIC=$STATIC cmake . \
    && make


# no need to install Boost IO and MKL here (static linking)
FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC
ARG CLOUDSDK_VERSION

RUN apt-get update -qqy && apt-get -qqy upgrade && apt-get install -qqy \
        curl \
        gcc \
        python3-dev \
        python3-crcmod \
        python3-pip \
        apt-transport-https \
        lsb-release \
        openssh-client \
        gnupg \
        unzip && \
    curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg  | gpg --dearmor -o /etc/apt/keyrings/google-cloud-cli.gpg && \
    echo "deb [signed-by=/etc/apt/keyrings/google-cloud-cli.gpg] https://packages.cloud.google.com/apt "cloud-sdk" main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    apt-get update && apt-get install -y google-cloud-cli=${CLOUDSDK_VERSION}-0 && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set storage/parallel_composite_upload_enabled False && \
    gcloud --version

RUN apt-get install -qqy \
        libgomp1 \
        gfortran-9 \
        libz-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libdeflate-dev

ARG REGENIE_VERSION
COPY --from=builder1 /src/regenie-$REGENIE_VERSION/regenie /usr/local/bin

# Download plink2 (linux AVX2 stable build)
ARG PLINK2_URL=https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_latest.zip
RUN curl -fsSL "${PLINK2_URL}" -o /tmp/plink2.zip \
    && unzip /tmp/plink2.zip -d /tmp/plink_extracted/\
    && mv /tmp/plink_extracted/plink2 /usr/local/bin/ \
    && rm -r /tmp/plink2.zip /tmp/plink_extracted/\
    && chmod +x /usr/local/bin/plink2
