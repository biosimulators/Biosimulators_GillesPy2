# Base OS
FROM python:3.7.9-slim-buster

# metadata
LABEL base_image="python:3.7.9-slim-buster"
LABEL version="0.0.1"
LABEL software="gillespy2"
LABEL software.version="1.5.6"
LABEL about.summary="GillesPy2 is a Python 3 package for stochastic simulation of biochemical systems"
LABEL about.home="https://github.com/GillesPy2/GillesPy2"
LABEL about.documentation="https://gillespy2.github.io/GillesPy2/https://gillespy2.github.io/GillesPy2/"
LABEL about.license_file="https://raw.githubusercontent.com/GillesPy2/GillesPy2/main/LICENSE"
LABEL about.license="GPL-3.0-only"
LABEL about.tags="systems biology,biochemical networks,dynamical modeling,stochastic simulation,SBML,SED-ML,COMBINE,OMEX,BioSimulators"
LABEL maintainer="BioSimulators Team <info@biosimulators.org>"

# Copy code for command-line interface into image and install it
COPY . /root/Biosimulators_gillespy2
RUN pip install /root/Biosimulators_gillespy2

# Entrypoint
ENTRYPOINT ["gillespy2"]
CMD []
