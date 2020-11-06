# Biosimulators_GillesPy2
BioSimulators-compliant command-line interface and Docker image for the [Gillespy2](https://stochss.github.io/GillesPy2) simulation program.

## Contents
* [Installation](#installation)
* [Usage](#usage)
* [License](#license)
* [Development team](#development-team)
* [Questions and comments](#questions-and-comments)

## Installation

### Install Python package
```
pip install git+https://github.com/biosimulators/Biosimulators_GillesPy2
```

### Install Docker image
```
docker pull ghcr.io/biosimulators/gillespy2
```

## Local usage
```
usage: gillespy2 [-h] [-d] [-q] -i ARCHIVE [-o OUT_DIR] [-v]

BioSimulators-compliant command-line interface to the Gillespy2 <https://stochss.github.io/GillesPy2> simulation program.

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

## Usage through Docker container
```
docker run \
  --tty \
  --rm \
  --mount type=bind,source="$(pwd)"/tests/fixtures,target=/root/in,readonly \
  --mount type=bind,source="$(pwd)"/tests/results,target=/root/out \
  ghcr.io/biosimulators/gillespy2:latest \
    -i /root/in/BIOMD0000000297.edited.omex \
    -o /root/out
```

## License
This package is released under the [MIT](LICENSE).

## Development team
This package was developed by the [Karr Lab](https://www.karrlab.org) at the Icahn School of Medicine at Mount Sinai and the [Center for Reproducible Biomedical Modeling](http://reproduciblebiomodels.org/).

## Questions and comments
Please contact [Bilal Shaikh](mailto:BilalShaikh42@gmail.com) or the [BioSimulators Team](mailto:info@biosimulators.org) with any questions or comments.
