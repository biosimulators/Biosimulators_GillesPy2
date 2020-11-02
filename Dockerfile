# Base OS
FROM ubuntu:18.04

# metadata
LABEL base_image="ubuntu:18.04"
LABEL version="1.0.0"
LABEL software="gillespy2"
LABEL software.version="1.5.3"
LABEL about.summary="GillesPy2 is a Python 3 package for stochastic simulation of biochemical systems"
LABEL about.home="https://github.com/GillesPy2/GillesPy2"
LABEL about.documentation="https://gillespy2.github.io/GillesPy2/https://gillespy2.github.io/GillesPy2/"
LABEL about.license_file="https://raw.githubusercontent.com/GillesPy2/GillesPy2/main/LICENSE"
LABEL about.license="GPL-3.0-only"
LABEL about.tags="systems biology,biochemical networks,dynamical modeling,stochastic simulation,SBML,SED-ML,COMBINE,OMEX,BioSimulators"
LABEL maintainer="BioSimulators Team <info@biosimulators.org>"

# Install requirements
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
    && pip3 install -U pip \
    && pip3 install -U setuptools \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# Copy code for command-line interface into image and install it
COPY . /root/Biosimulators_gillespy2
RUN pip3 install /root/Biosimulators_gillespy2

# Entrypoint
ENTRYPOINT ["gillespy2"]
CMD []
