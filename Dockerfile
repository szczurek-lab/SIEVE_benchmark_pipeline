FROM condaforge/mambaforge:latest as conda_bin

FROM mambaorg/micromamba:jammy

LABEL maintainer="Senbai Kang <kang@mimuw.edu.pl>"

ENV ENV_NAME snake
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV DEBIAN_FRONTEND noninteractive

ARG JAVA_VERSION=8
ARG CONAN_VERSION=1.54.0
ARG SIEVE_SIMULATOR_VERSION=v1.3.0
ARG DATAFILTER_VERSION=v0.1.0
ARG SCIPHI_VERSION=v0.1.7
ARG BEAST_2_VERSION_MAJOR=2.6
ARG BEAST_2_VERSION_MINOR=7
ARG SIEVE_VERSION=v0.15.6
ARG CELLPHY_VERSION=v0.9.2

USER root

# Copy conda binaries from condaforge/mambaforge:latest as it is required by snakemake --use-conda
COPY --from=conda_bin /opt/conda /opt/conda

WORKDIR /tmp

# Install system dependencies
RUN apt-get update && \
    apt-get install --no-install-recommends -y \
    software-properties-common \
    dirmngr \
    apt-utils \
    locales \
    wget \
    curl \
    vim \
    htop \
    ca-certificates \
    apt-transport-https \
    build-essential \
    autotools-dev \
    libicu-dev \
    libbz2-dev \
    libboost-all-dev \
    libomp-dev \
    gsfonts \
    gnupg2 && \
    apt clean && \
	rm -rf /var/lib/apt/lists/*

# Add R repository
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
	locale-gen en_US.utf8 && \
	/usr/sbin/update-locale LANG=en_US.UTF-8 && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    add-apt-repository ppa:c2d4u.team/c2d4u4.0+ && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

# Add zulu Java repository
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 0xB1998361219BD9C9 && \
    curl -O https://cdn.azul.com/zulu/bin/zulu-repo_1.0.0-3_all.deb && \
    dpkg -i zulu-repo_1.0.0-3_all.deb && \
    rm zulu-repo_1.0.0-3_all.deb

# Install system dependencies
RUN apt-get update && \
    apt-get install --no-install-recommends -y \
    cmake \
    git \
    zulu${JAVA_VERSION}-ca-jre \
    r-base \
    r-base-core \
    r-base-dev \
    r-recommended \
    r-cran-dplyr \
    r-cran-optparse \
    r-cran-ape \
    r-cran-phangorn && \
    apt clean

# Create a micromamba environment containing snakemake
COPY environment.yml environment.yml
RUN micromamba create -y -f environment.yml && \
    micromamba clean -afy && \
    rm -f environment.yml

# Install and configure conan
RUN wget -O conan-ubuntu-64.deb https://github.com/conan-io/conan/releases/download/${CONAN_VERSION}/conan-ubuntu-64.deb && \
    dpkg -i conan-ubuntu-64.deb && \
    rm -f conan-ubuntu-64.deb && \
    conan profile new default --detect && \
    conan profile update settings.compiler.libcxx=libstdc++11 default && \
    conan profile update env.CC=/usr/bin/gcc default && \
    conan profile update env.CXX=/usr/bin/g++ default

WORKDIR /root/pkgs

# Install SIEVE_simulator
RUN git clone --depth 1 --branch ${SIEVE_SIMULATOR_VERSION} https://github.com/szczurek-lab/SIEVE_simulator.git SIEVE_simulator && \
    cd SIEVE_simulator && \
    mkdir build && \
    cd build && \
    cmake -D CMAKE_C_COMPILER=/usr/bin/gcc -D CMAKE_CXX_COMPILER=/usr/bin/g++ -D CMAKE_BUILD_TYPE=Release .. && \
    cmake --build . && \
    ln bin/SIEVE_simulator /usr/local/bin/SIEVE_simulator && \
    cd ../.. && \
    rm -rf SIEVE_simulator

# Install DataFilter
RUN git clone --depth 1 --branch ${DATAFILTER_VERSION} https://github.com/szczurek-lab/DataFilter.git && \
    cd DataFilter && \
    mkdir build && \
    cd build && \
    conan install .. -s build_type=Release --build=missing && \
    cmake -D CMAKE_C_COMPILER=/usr/bin/gcc -D CMAKE_CXX_COMPILER=/usr/bin/g++ -D CMAKE_BUILD_TYPE=Release .. && \
    cmake --build . && \
    ln bin/datafilter /usr/local/bin/datafilter && \
    cd ../.. && \
    rm -rf DataFilter

# Install SCIPhI
RUN git clone --depth 1 --branch ${SCIPHI_VERSION} --recurse-submodules https://github.com/cbg-ethz/SCIPhI.git && \
    cd SCIPhI && \
    mkdir build && \
    cd build && \
    cmake .. && \
    cmake -D CMAKE_C_COMPILER=/usr/bin/gcc -D CMAKE_CXX_COMPILER=/usr/bin/g++ -D CMAKE_BUILD_TYPE=Release .. && \
    cmake --build . && \
    ln sciphi /usr/local/bin/sciphi && \
    cd ../.. && \
    rm -rf SCIPhI

# Install BEAST 2 and SIEVE
RUN wget -O BEAST.v${BEAST_2_VERSION_MAJOR}.${BEAST_2_VERSION_MINOR}.Linux.tgz https://github.com/CompEvol/beast2/releases/download/v${BEAST_2_VERSION_MAJOR}.${BEAST_2_VERSION_MINOR}/BEAST.v${BEAST_2_VERSION_MAJOR}.${BEAST_2_VERSION_MINOR}.Linux.tgz && \
    tar -xzf BEAST.v${BEAST_2_VERSION_MAJOR}.${BEAST_2_VERSION_MINOR}.Linux.tgz && \
    rm -f BEAST.v${BEAST_2_VERSION_MAJOR}.${BEAST_2_VERSION_MINOR}.Linux.tgz && \
    ln -s "${PWD}/beast/bin/beast" /usr/local/bin/beast && \
    ln -s "${PWD}/beast/bin/applauncher" /usr/local/bin/applauncher && \
    ln -s "${PWD}/beast/bin/packagemanager" /usr/local/bin/packagemanager && \
    packagemanager -add ORC && \
    mkdir -p /root/.beast/${BEAST_2_VERSION_MAJOR}/SIEVE && \
    cd /root/.beast/${BEAST_2_VERSION_MAJOR}/SIEVE && \
    wget -O SIEVE-${SIEVE_VERSION}.tar.gz https://codeload.github.com/szczurek-lab/SIEVE/tar.gz/refs/tags/${SIEVE_VERSION} && \
    tar -xzf SIEVE-${SIEVE_VERSION}.tar.gz && \
    rm -f SIEVE-${SIEVE_VERSION}.tar.gz && \
    cp SIEVE*/dist/SIEVE.${SIEVE_VERSION}.zip ${SIEVE_VERSION}.zip && \
    rm -rf SIEVE* && \
    unzip ${SIEVE_VERSION}.zip && \
    rm -f ${SIEVE_VERSION}.zip

# Install CellPhy
RUN git clone --depth 1 --branch ${CELLPHY_VERSION} https://github.com/amkozlov/cellphy.git && \
    ln cellphy/bin/raxml-ng-cellphy-linux /usr/local/bin/raxml-ng-cellphy-linux && \
    rm -rf cellphy

# Install SiFit
RUN git clone --depth 1 https://github.com/KChen-lab/SiFit.git && \
    mv SiFit/SiFit.jar /usr/local/lib/SiFit.jar && \
    rm -rf SiFit

VOLUME /root/data

WORKDIR /root/data

# ENTRYPOINT is provided in the base image of micromamba to activate the conda environment and is inherited here.
# CMD is used to provide the default parameters to the ENTRYPOINT.
CMD ["snakemake", "--cores", "all", "--use-conda", "--rerun-incomplete", "--keep-going", "--printshellcmds", "--rerun-triggers", "mtime"]