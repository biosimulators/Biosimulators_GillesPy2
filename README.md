[![Latest release](https://img.shields.io/github/v/tag/biosimulators/Biosimulators_GillesPy2)](https://github.com/biosimulations/Biosimulators_GillesPy2/releases)
[![PyPI](https://img.shields.io/pypi/v/biosimulators_gillespy2)](https://pypi.org/project/biosimulators_gillespy2/)
[![CI status](https://github.com/biosimulators/Biosimulators_GillesPy2/workflows/Continuous%20integration/badge.svg)](https://github.com/biosimulators/Biosimulators_GillesPy2/actions?query=workflow%3A%22Continuous+integration%22)
[![Test coverage](https://codecov.io/gh/biosimulators/Biosimulators_GillesPy2/branch/dev/graph/badge.svg)](https://codecov.io/gh/biosimulators/Biosimulators_GillesPy2)
[![All Contributors](https://img.shields.io/github/all-contributors/biosimulators/Biosimulators_GillesP2/HEAD)](#contributors-)

# BioSimulators-GillesPy2
BioSimulators-compliant command-line interface and Docker image for the [GillesPy2](https://stochss.github.io/GillesPy2) simulation program.

This command-line interface and Docker image enable users to use GillesPy2 to execute [COMBINE/OMEX archives](https://combinearchive.org/) that describe one or more simulation experiments (in [SED-ML format](https://sed-ml.org)) of one or more models (in [SBML format](http://sbml.org])).

A list of the algorithms and algorithm parameters supported by GillesPy2 is available at [BioSimulators](https://biosimulators.org/simulators/gillespy2).

A simple web application and web service for using GillesPy2 to execute COMBINE/OMEX archives is also available at [runBioSimulations](https://run.biosimulations.org).

## Installation

### Dependencies

* Python >= 3.7
* pip
* build-essential

### Install Python package
```
pip install biosimulators-gillespy2
```

### Install Docker image
```
docker pull ghcr.io/biosimulators/gillespy2
```

## Usage

### Local usage
```
usage: biosimulators-gillespy2 [-h] [-d] [-q] -i ARCHIVE [-o OUT_DIR] [-v]

BioSimulators-compliant command-line interface to the GillesPy2 <https://stochss.github.io/GillesPy2> simulation program.

optional arguments:
  -h, --help            show this help message and exit
  -d, --debug           full application debug mode
  -q, --quiet           suppress all console output
  -i ARCHIVE, --archive ARCHIVE
                        Path to OMEX file which contains one or more SED-ML-
                        encoded simulation experiments
  -o OUT_DIR, --out-dir OUT_DIR
                        Directory to save outputs
  -v, --version         show program's version number and exit
```

### Usage through Docker container
The entrypoint to the Docker image supports the same command-line interface described above.

For example, the following command could be used to use the Docker image to execute the COMBINE/OMEX archive `./modeling-study.omex` and save its outputs to `./`.

```
docker run \
  --tty \
  --rm \
  --mount type=bind,source="$(pwd)",target=/root/in,readonly \
  --mount type=bind,source="$(pwd)",target=/root/out \
  ghcr.io/biosimulators/gillespy2:latest \
    -i /root/in/modeling-study.omex \
    -o /root/out
```

## Documentation
Documentation is available at https://docs.biosimulators.org/Biosimulators_GillesPy2/.

## License
This package is released under the [MIT](LICENSE).

## Development team
This package was developed by the [Karr Lab](https://www.karrlab.org) at the Icahn School of Medicine at Mount Sinai and the [Center for Reproducible Biomedical Modeling](https://reproduciblebiomodels.org/) with assistance from the contributors listed [here](CONTRIBUTORS.md).

## Questions and comments
Please contact the [BioSimulators Team](mailto:info@biosimulators.org) with any questions or comments.
