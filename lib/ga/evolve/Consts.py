"""

:mod:`Consts` -- constants module
============================================================================

evolve have defaults in all genetic operators, settings and etc, this is an issue to helps the user in the API use and minimize the source code needed to make simple things. In the module :mod:`Consts`, you will find those defaults settings. You are encouraged to see the constants, but not to change directly on the module, there are methods for this.

General constants
----------------------------------------------------------------------------

.. attribute:: CDefPythonRequire

    The mininum version required to run evolve.

.. attribute:: CDefLogFile

    The default log filename.

.. attribute:: CDefLogLevel

    Default log level.

.. attribute:: CDefESCKey

    The ESC key ASCII code. Used to start Interactive Mode.

.. attribute:: CDefRangeMin

    Minimum range. This constant is used as integer and real max/min.

.. attribute:: CDefRangeMax

    Maximum range. This constant is used as integer and real max/min.

.. attribute:: CDefBroadcastAddress

    The broadcast address for UDP, 255.255.255.255

.. attribute:: CDefImportList

    The import list and messages

.. attribute:: nodeType

    The genetic programming node types, can be "TERMINAL":0 or "NONTERMINAL":1

.. attribute:: CDefGPGenomes

    The classes which are used in Genetic Programming, used to detected the
    correct mode when starting the evolution

Selection methods constants (:mod:`Selectors`)
----------------------------------------------------------------------------

.. attribute:: CDefTournamentPoolSize

    The default pool size for the Tournament Selector (:func:`Selectors.GTournamentSelector`).

Scaling scheme constants (:mod:`Scaling`)
----------------------------------------------------------------------------

.. attribute:: CDefScaleLinearMultiplier

    The multiplier of the Linear (:func:`Scaling.LinearScaling`) scaling scheme.

.. attribute:: CDefScaleSigmaTruncMultiplier

    The default Sigma Truncation (:func:`Scaling.SigmaTruncScaling`) scaling scheme.

.. attribute:: CDefScalePowerLawFactor

    The default Power Law (:func:`Scaling.PowerLawScaling`) scaling scheme factor.

.. attribute:: CDefScaleBoltzMinTemp

    The default mininum temperature of the (:func:`Scaling.BoltzmannScaling`) scaling scheme factor.

.. attribute:: CDefScaleBoltzFactor

    The default Boltzmann Factor of (:func:`Scaling.BoltzmannScaling`) scaling scheme factor.
    This is the factor that the temperature will be subtracted.

.. attribute:: CDefScaleBoltzStart

    The default Boltzmann start temperature (:func:`Scaling.BoltzmannScaling`).
    If you don't set the start temperature parameter, this will be the default initial
    temperature for the Boltzmann scaling scheme.

Population constants (:class:`GPopulation.GPopulation`)
----------------------------------------------------------------------------

.. attribute:: CDefPopScale

    Default scaling scheme.


GA Engine constants (:class:`GSimpleGA.GSimpleGA`)
----------------------------------------------------------------------------

.. attribute:: CDefGAGenerations

    Default number of generations.

.. attribute:: CDefGAMutationRate

    Default mutation rate.

.. attribute:: CDefGACrossoverRate

    Default crossover rate.

.. attribute:: CDefGAPopulationSize

    Default population size.

.. attribute:: CDefGASelector

    Default selector method.

"""
import Scaling
import Selectors
import logging
import Util

# Required python version 2.5+
CDefPythonRequire = (2, 5)

# Logging system
CDefLogFile = "evolve.log"
# Optimization type

CDefESCKey = 27

CDefImportList = {"visual.graph": "you must install VPython !",
                        "csv": "csv module not found !",
                        "sqlite3": "sqlite3 module not found, are you using Jython or IronPython ?",
                        "xmlrpclib": "xmlrpclib module not found !",
                        "pydot": "Pydot module not found, you must install Pydot to plot graphs !"}

####################
# Defaults section #
####################

# - Tournament selector
CDefTournamentPoolSize = 2

# - Scale methods defaults
CDefScaleLinearMultiplier = 1.2
CDefScaleSigmaTruncMultiplier = 2.0
CDefScalePowerLawFactor = 1.0005
CDefScaleBoltzMinTemp = 1.0
CDefScaleBoltzFactor = 0.05
# 40 temp. = 500 generations
CDefScaleBoltzStart = 40.0

# - Population Defaults
CDefPopScale = Scaling.NoScaling
CDefRawPopSort = Util.RawSorting

# - GA Engine defaults
CDefGAGenerations = 100
CDefGAMutationRate = 0.3
CDefGACrossoverRate = 0.7
CDefGASelector = Selectors.GRankSelector
CDefGAElitismReplacement = 1

# - This is general used by integer/real ranges defaults
CDefRangeMin = 0
CDefRangeMax = 100

# Gaussian Gradient
CDefGaussianGradientMU = 1.0
CDefGaussianGradientSIGMA = (1.0 / 3.0)  # approx. +/- 3-sigma is +/- 10%

# - DB Adapters URL Post defaults
CDefURLPostStatsGenFreq = 100

# - DB Adapters CSV File defaults
CDefCSVFileName = "evolve.csv"
CDefCSVFileStatsGenFreq = 1