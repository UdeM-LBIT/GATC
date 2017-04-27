"""
:mod:`evolve` -- the main evolve namespace
================================================================

This is the main module of the evolve, every other module
is above this namespace, for example, to import :mod:`Mutators`:

This code was forked from Pyevolve : https://github.com/perone/Pyevolve
written by Christian S. Perone. Modifications were made on almost all files
of the package to accomodate it for GaPol.

"""
__all__ = ["Consts", "FunctionSlot",
                     "GenomeBase", "GPopulation",
                     "GSimpleGA", "Scaling", "Selectors",
                     "Statistics", "Util"]

from . import Consts
import sys

if sys.version_info[:2] < Consts.CDefPythonRequire:
    raise Exception("Python 2.5+ required, the version %s was found on your system !" % (sys.version_info[:2],))

del sys

def logEnable(filename=Consts.CDefLogFile, level=Consts.CDefLogLevel):
    """ Enable the log system for evolve
    :param filename: the log filename
    :param level: the debugging level

    Example:
           >>> evolve.logEnable()

    """
    import logging
    logging.basicConfig(level=level,
                        format='%(asctime)s [%(module)s:%(funcName)s:%(lineno)d] %(levelname)s %(message)s',
                        filename=filename,
                        filemode='w')
    logging.info("Log was enabled by user.")
