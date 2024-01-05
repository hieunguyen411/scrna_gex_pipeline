FROM rocker/rstudio:latest
RUN apt-get update && apt-get upgrade -y
RUN apt-get install build-essential -y
RUN apt install libudunits2-dev libgdal-dev libgeos-dev libproj-dev libfontconfig1-dev -y
RUN apt-get clean all && \
  apt-get update && \
  apt-get upgrade -y && \
  apt-get install -y \
    libhdf5-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    libxt-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libglpk40 \
    libgit2-28 \
  && apt-get clean all && \
  apt-get purge && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt install libfontconfig1-dev 
RUN apt-get update && apt-get install libmagick++-6.q16-dev -y
RUN apt-get install -y libgsl23

RUN mkdir -p installation_src
COPY install_R_packages.R ./installation_src/
RUN update-alternatives --set libblas.so.3-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/blas/libblas.so.3
RUN update-alternatives --set liblapack.so.3-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3

RUN Rscript ./installation_src/install_R_packages.R

RUN mkdir -p /home/storage
RUN mkdir -p /home/src
RUN mkdir -p /home/src_pipeline
RUN mkdir -p /home/outdir

