FROM r-base:3.5.2

RUN apt-get update --fix-missing -qq && apt-get install -y -q \
        libssl-dev \
        libcurl4-openssl-dev \
        libnlopt-dev \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## Install R packages
RUN R -e 'install.packages(c("car", "CompQuadForm", "testthat"), repos="http://cloud.r-project.org/")'
