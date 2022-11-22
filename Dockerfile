FROM ubuntu:22.04

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update --yes \
  && apt-get install --yes --no-install-recommends \
  wget \
  locales \
  tar \
  gpg

# R Package dependencies
RUN apt-get install --yes --no-install-recommends \
  git \
  libgmp3-dev \
  default-jdk \
  make \
  libcurl4-openssl-dev \
  libicu-dev \
  libssl-dev \
  libpng-dev \
  libjpeg-dev \
  libxml2-dev \
  libglpk-dev \
  zlib1g-dev \
  libcairo2-dev \
  sudo \
  pandoc

ENV TZ UTC

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

RUN export _RIG_VERSION=0.5.2 \
  && wget https://github.com/r-lib/rig/releases/download/v${_RIG_VERSION}/rig-linux-${_RIG_VERSION}.tar.gz \
  && mkdir -p /data \
  && tar -C /data/ -xf rig-linux-${_RIG_VERSION}.tar.gz \
  && ln -sf /data/bin/rig /usr/local/bin/rig \
  && rm -rf /rig-linux-${_RIG_VERSION}.tar.gz

RUN rig add release

RUN addgroup --system app \
    && adduser --disabled-password \
    --gecos "" \
    --system \
    --ingroup app \
    app \
    && adduser app sudo

RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

RUN chown app:app -R /home/app

USER app

ENV TZ UTC
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

WORKDIR /home/app/pdxplorer

RUN R -q -s -e \
  'options(repos="https://packagemanager.rstudio.com/all/__linux__/jammy/latest"); \
  install.packages("pak");'

RUN R -q -s -e \
  'options(repos="https://packagemanager.rstudio.com/all/__linux__/jammy/latest"); \
  pak::pkg_install(pkg=c( \
  "shiny", \
  "shinydashboard", \
  "shinyjs", \
  "shinythemes", \
  "bslib", \
  "rhino", \
  "cachem", \
  "data.table", \
  "reshape2", \
  "fs", \
  "dplyr", \
  "tidyr", \
  "readr", \
  "ggdendro", \
  "ggplot2", \
  "grid", \
  "gridExtra", \
  "DT", \
  "scales", \
  "pheatmap", \
  "RColorBrewer", \
  "bioc::maftools", \
  "bioc::limma", \
  "bioc::DESeq2", \
  "bioc::scran", \
  "bioc::systemPipeR", \
  "bioc::EnhancedVolcano", \
  "bioc::chimeraviz", \
  "github::statistikat/codeModules", \
  "github::luciorq/ABCutilities" \
  ),upgrade=TRUE,ask=FALSE);'

USER root

RUN rm -rf /tmp/downloaded_packages/ \
  && rm -rf /tmp/*.rds \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /data

USER app

ENV TZ UTC

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

COPY --chown=app:app app /home/app/pdxplorer/app

COPY --chown=app:app data /home/app/pdxplorer/data

EXPOSE 80

WORKDIR /home/app/pdxplorer

CMD ["R", "-q", "-s", "-e", "shiny::runApp(appDir='/home/app/pdxplorer/app',host='0.0.0.0',port=80,launch.browser=FALSE)"]