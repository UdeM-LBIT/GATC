"""

:mod:`Scaling` -- scaling schemes module
===========================================================

This module have the *scaling schemes* like Linear scaling, etc.

"""
import Consts
import Util
import math
import logging
import numpy as np

def LinearScaling(pop, sp=0):
    """ Linear Scaling scheme

    .. warning :: Linear Scaling is only for positive raw scores

    """
    logging.debug("Running linear scaling.")
    pop.statistics(ftfunc=sp)
    c = Consts.CDefScaleLinearMultiplier
    a = b = delta = 0.0

    pop_rawAve = pop.stats["rawAve"]
    pop_rawMax = pop.stats["rawMax"]
    pop_rawMin = pop.stats["rawMin"]

    if pop_rawAve == pop_rawMax:
        a = 1.0
        b = 0.0
    elif pop_rawMin > (c * pop_rawAve - pop_rawMax / c - 1.0):
        delta = pop_rawMax - pop_rawAve
        a = (c - 1.0) * pop_rawAve / delta
        b = pop_rawAve * (pop_rawMax - (c * pop_rawAve)) / delta
    else:
        delta = pop_rawAve - pop_rawMin
        a = pop_rawAve / delta
        b = -pop_rawMin * pop_rawAve / delta

    for i in xrange(len(pop)):
        f = pop[i].score[sp]
        if f < 0.0:
            Util.raiseException("Score %r is negative, linear scaling not supported !" % (f,), ValueError)
        f = f * a + b
        if f < 0:
            f = 0.0
        pop[i].fitness = f


def NoScaling(pop):
    """Fitness is raw score"""
    for i, p in enumerate(pop):
        pop[i].fitness = sum(p.score)


def SigmaTruncScaling(pop, sp=0):
    """ Sigma Truncation scaling scheme, allows negative scores """
    logging.debug("Running sigma truncation scaling.")
    pop.statistics(ftfunc=sp)
    c = Consts.CDefScaleSigmaTruncMultiplier
    pop_rawAve = pop.stats["rawAve"]
    pop_rawDev = pop.stats["rawDev"]
    for i in xrange(len(pop)):
        f = pop[i].score[sp] - pop_rawAve
        f += c * pop_rawDev
        if f < 0:
            f = 0.0
        pop[i].fitness = f

def PowerLawScaling(pop, sp=0):
    """ Power Law scaling scheme

    .. warning :: Power Law Scaling is only for positive raw scores

    """
    logging.debug("Running power law scaling.")
    k = Consts.CDefScalePowerLawFactor
    for i in xrange(len(pop)):
        f = pop[i].score[sp]
        if f < 0.0:
            Util.raiseException("Score %r is negative, power law scaling not supported !" % (f,), ValueError)
        f = math.pow(f, k)
        pop[i].fitness = f


def BoltzmannScaling(pop, sp=0):
    """ Boltzmann scaling scheme. You can specify the **boltz_temperature** to the
    population parameters, this parameter will set the start temperature. You
    can specify the **boltz_factor** and the **boltz_min** parameters, the **boltz_factor**
    is the value that the temperature will be subtracted and the **boltz_min** is the
    mininum temperature of the scaling scheme.

    .. versionadded: 0.6
        The `BoltzmannScaling` function.

    """
    boltz_temperature = pop.getParam("boltz_temperature", Consts.CDefScaleBoltzStart)
    boltz_factor = pop.getParam("boltz_factor", Consts.CDefScaleBoltzFactor)
    boltz_min = pop.getParam("boltz_min", Consts.CDefScaleBoltzMinTemp)

    boltz_temperature -= boltz_factor
    boltz_temperature = max(boltz_temperature, boltz_min)
    pop.setParams(boltzTemperature=boltz_temperature)

    boltz_e = []
    avg = 0.0

    for i in xrange(len(pop)):
        val = math.exp(pop[i].score[sp] / boltz_temperature)
        boltz_e.append(val)
        avg += val

    avg /= len(pop)

    for i in xrange(len(pop)):
        pop[i].fitness = boltz_e[i] / avg

def ExponentialScaling(pop, sp=0):
    """ Exponential Scaling Scheme. The fitness will be the same as (e^score).

    .. versionadded: 0.6
        The `ExponentialScaling` function.
    """
    for i in xrange(len(pop)):
        pop[i].fitness = math.exp(pop[i].score[sp])

def SaturatedScaling(pop, sp=0):
    """ Saturated Scaling Scheme. The fitness will be the same as 1.0-(e^score)

    .. versionadded: 0.6
        The `SaturatedScaling` function.
    """
    for i in xrange(len(pop)):
        pop[i].fitness = 1.0 - math.exp(pop[i].score[sp])


def sigmoid(x, sfunc='tanh'):
    if sfunc == 'exp':
        return 1 / (1 + np.exp(-x))
    else:
        return (math.tanh(x) + 1) /2


def LklCostScaling(pop):
    """This should not be used"""
    stat = pop.getStatistics()
    lkldict = stat.getRawScore(0)
    costdict = stat.getRawScore(1) 
    a = (costdict["rawMax"] - costdict["rawMin"])/ (lkldict["rawMax"] - lkldict["rawMin"])
    b = costdict["rawMax"] - a*lkldict["rawMax"]
    for p in pop:
        normalized_lklcost = (p.score[0]*a + b)
        p.fitness =  normalized_lklcost + p.score[1]


def __scaleLKLToZero(pop, keepraw=True):
    stat = pop.getStatistics()
    lkldict = stat.getRawScore(0)
    if keepraw:
        min_lkl = 0
    else:
        min_lkl =  lkldict["rawMin"]
    lklscaled_vec = []
    for p in pop:
        x =  p.score[0] - min_lkl
        lklscaled_vec.append(x)
        p.setParams(lklscaled=x) # should be the same if keepraw
    return lklscaled_vec


def WeightSigmoidScaling(pop, weight=[], keepraw=True):
    
    lklscaled = __scaleLKLToZero(pop, keepraw)    
    if not weight:
        weight =  pop.getOrSetWeight(lklscaled)
    for p in pop:
        p.fitness = sigmoid(weight[0]*p.getParam('lklscaled'), 'exp') + sigmoid(weight[1]*p.score[1], 'exp')


def WeightScaling(pop, weight=[], keepraw=True):
    # fix the weight for all generation
    lklscaled = __scaleLKLToZero(pop, keepraw)
    if not weight:
        weight =  pop.getOrSetWeight(lklscaled)
    for p in pop:
        p.fitness = np.dot(weight, [p.getParam('lklscaled'), p.score[1]])