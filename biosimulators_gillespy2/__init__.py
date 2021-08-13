from ._version import __version__  # noqa: F401
# :obj:`str`: version

from .core import exec_sedml_docs_in_combine_archive  # noqa: F401

import gillespy2


def get_simulator_version():
    """ Get the version of GillesPy2

    Returns:
        :obj:`str`: version
    """
    return gillespy2.__version__
