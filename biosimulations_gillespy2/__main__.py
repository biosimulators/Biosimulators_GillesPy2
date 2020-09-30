""" BioSimulators-compliant command-line interface to the `Gillespy2 <https://github.com/GillesPy2/GillesPy2>`_ simulation program.

:Author: Bilal Shaikh <bilalshaikh42@gmail.com>
:Date: 2020-08-18
:Copyright: 2020, Biosimulations Developers
:License: MIT
"""

from .core import exec_combine_archive
import Biosimulators_gillespy2
import cement


class BaseController(cement.Controller):
    """ Base controller for command line application """

    class Meta:
        label = 'base'
        description = ("BioSimulators-compliant command-line interface to the "
                       "<Gillespy2> simulation program <https://github.com/GillesPy2/GillesPy2>.")
        help = "gillespy2"
        arguments = [
            (['-i', '--archive'], dict(type=str,
                                       required=True,
                                       help='Path to OMEX file which contains one or more SED-ML-encoded simulation experiments')),
            (['-o', '--out-dir'], dict(type=str,
                                       default='.',
                                       help='Directory to save outputs')),
            (['-v', '--version'], dict(action='version',
                                       version=Biosimulators_gillespy2.__version__)),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        exec_combine_archive(args.archive, args.out_dir)


class App(cement.App):
    """ Command line application """
    class Meta:
        label = 'gillespy2'
        base_controller = 'base'
        handlers = [
            BaseController,
        ]


def main():
    with App() as app:
        app.run()
