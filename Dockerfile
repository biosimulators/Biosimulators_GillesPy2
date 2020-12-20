# Base OS
FROM python:3.7.9-slim-buster

ARG VERSION=0.1.10
ARG SIMULATOR_VERSION="1.5.7"

# metadata
LABEL \
    org.opencontainers.image.title="GillesPy2" \
    org.opencontainers.image.version="${SIMULATOR_VERSION}" \
    org.opencontainers.image.description="Python 3 package for stochastic simulation of biochemical systems" \
    org.opencontainers.image.url="https://github.com/StochSS/GillesPy2" \
    org.opencontainers.image.documentation="https://stochss.github.io/GillesPy2/" \
    org.opencontainers.image.source="https://github.com/biosimulators/Biosimulators_GillesPy2" \
    org.opencontainers.image.authors="BioSimulators Team <info@biosimulators.org>" \
    org.opencontainers.image.vendor="BioSimulators Team" \
    org.opencontainers.image.licenses="GPL-3.0-only" \
    \
    base_image="python:3.7.9-slim-buster" \
    version="${VERSION}" \
    software="gillespy2" \
    software.version="${SIMULATOR_VERSION}" \
    about.summary="Python 3 package for stochastic simulation of biochemical systems" \
    about.home="https://github.com/StochSS/GillesPy2" \
    about.documentation="https://stochss.github.io/GillesPy2/" \
    about.license_file="https://raw.githubusercontent.com/stochss/GillesPy2/main/LICENSE" \
    about.license="SPDX:GPL-3.0-only" \
    about.tags="systems biology,biochemical networks,dynamical modeling,stochastic simulation,SBML,SED-ML,COMBINE,OMEX,BioSimulators" \
    maintainer="BioSimulators Team <info@biosimulators.org>"

# Install requirements
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        build-essential \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# Copy code for command-line interface into image and install it
COPY . /root/Biosimulators_gillespy2
RUN pip install /root/Biosimulators_gillespy2
RUN pip install "gillespy2==${SIMULATOR_VERSION}"

# Entrypoint
ENTRYPOINT ["gillespy2"]
CMD []
