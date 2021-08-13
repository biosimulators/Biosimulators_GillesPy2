""" BioSimulators-compliant command-line interface to the `GillesPy2 <https://github.com/StochSS/GillesPy2>`_ simulation program.

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Bilal Shaikh <bilalshaikh42@gmail.com>
:Date: 2020-12-13
:Copyright: 2020, BioSimulators Team
:License: MIT
"""

from . import get_simulator_version
from ._version import __version__
from .core import exec_sedml_docs_in_combine_archive
from biosimulators_utils.simulator.cli import build_cli

App = build_cli('biosimulators-gillespy2', __version__,
                'GillesPy2', get_simulator_version(), 'https://github.com/StochSS/GillesPy2',
                exec_sedml_docs_in_combine_archive)


def main():
    with App() as app:
        app.run()
