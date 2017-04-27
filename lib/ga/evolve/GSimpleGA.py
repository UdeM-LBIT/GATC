"""

:mod:`GSimpleGA` -- the genetic algorithm by itself
=====================================================================

This module contains the GA Engine, the GA Engine class is responsible
for all the evolutionary process. It contains the GA Engine related
funtions, like the Termination Criteria functions for convergence analysis, etc.

Default Parameters
-------------------------------------------------------------

*Number of Generations*

   Default is 100 generations

*Mutation Rate*

   Default is 0.02, which represents 2%

*Crossover Rate*

   Default is 0.9, which represents 90%

*Elitism Replacement*

   Default is 1 individual

*Population Size*

   Default is 80 individuals

*Minimax*

   >>> Consts.minimaxType["maximize"]

   Maximize the evaluation function

*Selector (Selection Method)*

   :func:`Selectors.GRankSelector`

   The Rank Selection method

Class
-------------------------------------------------------------

"""
import random
import math
import logging
from time import time
from types import BooleanType
from sys import stdout as sys_stdout

from .GPopulation import GPopulation
from .FunctionSlot import FunctionSlot
from .GenomeBase import GenomeBase
from . import Consts
from . import Util


def ConvergenceCriteria(ga_engine):
    """ Terminate the evolution when the population have converged

    Example:
       >>> ga_engine.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)

    """
    pop = ga_engine.getPopulation()
    return pop[0] == pop[len(pop) - 1]



def FitnessStatsCriteria(ga_engine):
    """ Terminate the evoltion based on the fitness stats

    Example:
       >>> ga_engine.terminationCriteria.set(GSimpleGA.FitnessStatsCriteria)


    """
    stats = ga_engine.getStatistics()
    if stats["fitMax"] == stats["fitMin"]:
        if stats["fitAve"] == stats["fitMax"]:
            return True
    return False

class GSimpleGA(object):
    """ GA Engine Class - The Genetic Algorithm Core

    Example:
       >>> ga = GSimpleGA.GSimpleGA(genome)
       >>> ga.selector.set(Selectors.GRouletteWheel)
       >>> ga.setGenerations(120)

    :param genome: the :term:`Sample Genome`
    :param seed: the random seed value

    .. note:: if you use the same random seed, all the runs of algorithm will be the same

    """

    selector = None
    """ This is the function slot for the selection method
    if you want to change the default selector, you must do this: ::

       ga_engine.selector.set(Selectors.GRouletteWheel) """

    stepCallback = None
    """ This is the :term:`step callback function` slot,
    if you want to set the function, you must do this: ::

       def your_func(ga_engine):
          # Here you have access to the GA Engine
          return False

       ga_engine.stepCallback.set(your_func)

    now *"your_func"* will be called every generation.
    When this function returns True, the GA Engine will stop the evolution and show
    a warning, if False, the evolution continues.
    """

    terminationCriteria = None
    """ This is the termination criteria slot, if you want to set one
    termination criteria, you must do this: ::

       ga_engine.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)

    Now, when you run your GA, it will stop when the termination criteria is met

    When this function returns True, the GA Engine will stop the evolution and show
    a warning, if False, the evolution continues, this function is called every
    generation.
    """

    def __init__(self, genomes, seed=None):
        """ Initializator of GSimpleGA """
        if seed:
            random.seed(seed)

        self.internalPop = None
        if isinstance(genomes, GenomeBase):
            self.internalPop = GPopulation(genomes, single=True)
        elif isinstance(genomes, GPopulation):
            self.internalPop = genomes
        else:
            Util.raiseException("genomes should be a single Genome from GenomeBase or a GPopulation", TypeError)

        
        self.nGenerations = Consts.CDefGAGenerations
        self.pMutation = Consts.CDefGAMutationRate
        self.pHGT = 0
        self.pCrossover = Consts.CDefGACrossoverRate
        self.nElitismReplacement = Consts.CDefGAElitismReplacement
        self.minimax = Consts.CDefPopMinimax
        self.elitism = True

        self.time_init = None
        self.max_time = None

        self.selector = FunctionSlot("Selector")
        self.stepCallback = FunctionSlot("Generation Step Callback")
        self.terminationCriteria = FunctionSlot("Termination Criteria")
        self.selector.set(Consts.CDefGASelector)
        self.allSlots = (self.selector, self.stepCallback, self.terminationCriteria)

        self.internalParams = {}

        self.currentGeneration = 0

        logging.debug("A GA Engine was created, nGenerations=%d", self.nGenerations)

    def __call__(self, *args, **kwargs):
        """ A method to implement a callable object

        Example:
           >>> ga_engine(freq_stats=10)

        .. versionadded:: 0.6
           The callable method.
        """
        if kwargs.get("freq_stats", None):
            return self.evolve(kwargs.get("freq_stats"))
        else:
            return self.evolve()

    def setParams(self, **args):
        """ Set the internal params

        Example:
           >>> ga.setParams(gp_terminals=['x', 'y'])


        :param args: params to save

        ..versionaddd:: 0.6
           Added the *setParams* method.
        """
        self.internalParams.update(args)

    def getParam(self, key, nvl=None):
        """ Gets an internal parameter

        Example:
           >>> ga.getParam("gp_terminals")
           ['x', 'y']

        :param key: the key of param
        :param nvl: if the key doesn't exist, the nvl will be returned

        ..versionaddd:: 0.6
           Added the *getParam* method.
        """
        return self.internalParams.get(key, nvl)

    def setElitismReplacement(self, numreplace):
        """ Set the number of best individuals to copy to the next generation on the elitism

        :param numreplace: the number of individuals

        .. versionadded:: 0.6
           The *setElitismReplacement* method.

        """
        if numreplace < 1:
            if numreplace > 0 and self.internalPop.popSize:
                numreplace = int(math.ceil(numreplace*self.internalPop.popSize))
            else:
                Util.raiseException("Replacement number must be >= 1", ValueError)
        self.nElitismReplacement = numreplace

    def __repr__(self):
        """ The string representation of the GA Engine """
        minimax_type = Consts.minimaxType.keys()[Consts.minimaxType.values().index(self.minimax)]
        ret = "- GSimpleGA\n"
        ret += "\tPopulation Size:\t %d\n" % self.internalPop.popSize
        ret += "\tGenerations:\t\t %d\n" % self.nGenerations
        ret += "\tCurrent Generation:\t %d\n" % self.currentGeneration
        ret += "\tMutation Rate:\t\t %.2f\n" % self.pMutation
        ret += "\tCrossover Rate:\t\t %.2f\n" % self.pCrossover
        ret += "\tMinimax Type:\t\t %s\n" % minimax_type.capitalize()
        ret += "\tElitism:\t\t %s\n" % self.elitism
        ret += "\tElitism Replacement:\t %d\n" % self.nElitismReplacement
        for slot in self.allSlots:
            ret += "\t" + slot.__repr__()
        ret += "\n"
        return ret

    def setMultiProcessing(self, flag=True, full_copy=False, max_processes=None):
        """ Sets the flag to enable/disable the use of python multiprocessing module.
        Use this option when you have more than one core on your CPU and when your
        evaluation function is very slow.

        evolve will automaticly check if your Python version has **multiprocessing**
        support and if you have more than one single CPU core. If you don't have support
        or have just only one core, evolve will not use the **multiprocessing**
        feature.

        evolve uses the **multiprocessing** to execute the evaluation function over
        the individuals, so the use of this feature will make sense if you have a
        truly slow evaluation function (which is commom in GAs).

        The parameter "full_copy" defines where the individual data should be copied back
        after the evaluation or not. This parameter is useful when you change the
        individual in the evaluation function.

        :param flag: True (default) or False
        :param full_copy: True or False (default)
        :param max_processes: None (default) or an integer value

        .. warning:: Use this option only when your evaluation function is slow, so you'll
                     get a good tradeoff between the process communication speed and the
                     parallel evaluation. The use of the **multiprocessing** doesn't means
                     always a better performance.

        .. note:: To enable the multiprocessing option, you **MUST** add the *__main__* check
                  on your application, otherwise, it will result in errors. See more on the
                  `Python Docs <http://docs.python.org/library/multiprocessing.html#multiprocessing-programming>`__
                  site.

        .. versionadded:: 0.6
           The `setMultiProcessing` method.

        """
        if type(flag) != BooleanType:
            Util.raiseException("Multiprocessing option must be True or False", TypeError)

        if type(full_copy) != BooleanType:
            Util.raiseException("Multiprocessing 'full_copy' option must be True or False", TypeError)

        self.internalPop.setMultiProcessing(flag, full_copy, max_processes)


    def setPopulationSize(self, size):
        """ Sets the population size, calls setPopulationSize() of GPopulation

        :param size: the population size

        .. note:: the population size must be >= 2

        """
        if size < 2:
            Util.raiseException("population size must be >= 2", ValueError)
        self.internalPop.setPopulationSize(size)


    def setScaleMethod(self, fn):
        """Setting scale method for the GA"""
        self.internalPop.setScaleMethod(fn)

    
    def setRawSortMethod(self, fn):
        """Setting scale method for the GA"""
        self.internalPop.setRawSortMethod(fn)


    def setMutationRate(self, rate):
        """ Sets the mutation rate, between 0.0 and 1.0

        :param rate: the rate, between 0.0 and 1.0

        """
        if (rate > 1.0) or (rate < 0.0):
            Util.raiseException("Mutation rate must be >= 0.0 and <= 1.0", ValueError)
        self.pMutation = rate


    def setTransferParam(self, param):
        """ Sets the mutation rate, between 0.0 and 1.0

        :param rate: the rate, between 0.0 and 1.0

        """
        self.pHGT = param.rate        
        self.internalPop.setScaleMethod(param.getScaling())
        self.internalPop.setParams(**param.__dict__)


    def setCrossoverRate(self, rate):
        """ Sets the crossover rate, between 0.0 and 1.0

        :param rate: the rate, between 0.0 and 1.0

        """
        if (rate > 1.0) or (rate < 0.0):
            Util.raiseException("Crossover rate must be >= 0.0 and <= 1.0", ValueError)
        self.pCrossover = rate

    def setGenerations(self, num_gens):
        """ Sets the number of generations to evolve

        :param num_gens: the number of generations

        """
        if num_gens < 1:
            Util.raiseException("Number of generations must be >= 1", ValueError)
        self.nGenerations = num_gens

    def getGenerations(self):
        """ Return the number of generations to evolve

        :rtype: the number of generations

        .. versionadded:: 0.6
           Added the *getGenerations* method
        """
        return self.nGenerations

    def getMinimax(self):
        """ Gets the minimize/maximize mode

        :rtype: the Consts.minimaxType type

        """
        return self.minimax

    def setMinimax(self, mtype):
        """ Sets the minimize/maximize mode, use Consts.minimaxType

        :param mtype: the minimax mode, from Consts.minimaxType

        """
        if mtype not in Consts.minimaxType.values():
            Util.raiseException("Minimax must be maximize or minimize", TypeError)
        self.minimax = mtype

    def getCurrentGeneration(self):
        """ Gets the current generation

        :rtype: the current generation

        """
        return self.currentGeneration

    def setElitism(self, flag):
        """ Sets the elitism option, True or False

        :param flag: True or False

        """
        if type(flag) != BooleanType:
            Util.raiseException("Elitism option must be True or False", TypeError)
        self.elitism = flag


    def setMaxTime(self, seconds):
        """ Sets the maximun evolve time of the GA Engine

        :param seconds: maximum time in seconds
        """
        self.max_time = seconds

    def getMaxTime(self):
        """ Get the maximun evolve time of the GA Engine

        :rtype: True or False
        """
        return self.max_time

    def bestIndividual(self):
        """ Returns the population best individual

        :rtype: the best individual

        """
        return self.internalPop.bestFitness()

    def worstIndividual(self):
        """ Returns the population worst individual

        :rtype: the best individual

        """
        return self.internalPop.worstFitness()


    def initialize(self):
        """ Initializes the GA Engine. Create and initialize population """
        self.internalPop.create(minimax=self.minimax)
        self.internalPop.initialize(ga_engine=self)
        logging.debug("The GA Engine was initialized !")

    def getPopulation(self):
        """ Return the internal population of GA Engine

        :rtype: the population (:class:`GPopulation.GPopulation`)

        """
        return self.internalPop

    def getStatistics(self):
        """ Gets the Statistics class instance of current generation

        :rtype: the statistics instance (:class:`Statistics.Statistics`)

        """
        return self.internalPop.getStatistics()

    def step(self):
        """ Just do one step in evolution, one generation """
        newPop = GPopulation(self.internalPop)
        logging.debug("Population was cloned.")

        size_iterate = len(self.internalPop)

        # Odd population size
        if size_iterate % 2 != 0:
            size_iterate -= 1

        crossover_empty = self.select(popID=self.currentGeneration).crossover.isEmpty()

        for i in xrange(0, size_iterate, 2):
            genomeMom = self.select(popID=self.currentGeneration)
            genomeDad = self.select(popID=self.currentGeneration)
            t = time()
            if not crossover_empty and self.pCrossover >= 1.0:
                for it in genomeMom.crossover.applyFunctions(mom=genomeMom, dad=genomeDad, count=2):
                    (sister, brother) = it
                    break
            else:
                if not crossover_empty and Util.randomFlipCoin(self.pCrossover):
                    for it in genomeMom.crossover.applyFunctions(mom=genomeMom, dad=genomeDad, count=2):
                        (sister, brother) = it
                        break
                else:
                    sister = genomeMom.clone()
                    brother = genomeDad.clone()
            
            #self.printTimeElapsed("Crossover")
            sister.mutate(pmut=self.pMutation, ga_engine=self)
            brother.mutate(pmut=self.pMutation, ga_engine=self)
            #self.printTimeElapsed("Mutation")
            newPop.internalPop.append(sister)
            newPop.internalPop.append(brother)

        if len(self.internalPop) % 2 != 0:
            genomeMom = self.select(popID=self.currentGeneration)
            genomeDad = self.select(popID=self.currentGeneration)

            if Util.randomFlipCoin(self.pCrossover):
                for it in genomeMom.crossover.applyFunctions(mom=genomeMom, dad=genomeDad, count=1):
                    (sister, brother) = it
                    break
            else:
                sister = random.choice([genomeMom, genomeDad])
                sister = sister.clone()
                sister.mutate(pmut=self.pMutation, ga_engine=self)

            newPop.internalPop.append(sister)

        logging.debug("Evaluating the new created population.")
        newPop.evaluate()
        # self.printTimeElapsed("Step")


        if self.getParam('fastconv', False):
            # lazyness overload
            t = time()
            logging.debug("Running in fastconv mode.")  
            # we are going to remove all transfer that are worst than the previous gen
            # in other words replace all transfered that increase cost then likelihood
            # transfer case
            lower_worst = max(1, int(math.log(len(newPop))))

            if self.pHGT > 0:
                # remove genome with worst cost and lkl compared to previous generation
                # we look at fitness only, in order to not sort multiple times
                if self.getMinimax() == Consts.minimaxType["maximize"]:
                    for i in xrange(lower_worst):
                        if self.internalPop.worstFitness(i).score[1] > newPop.worstFitness().score[1] and  self.internalPop.worstFitness(i).score[0] > newPop.worstFitness().score[0]:
                                newPop[-1] = self.internalPop.worstFitness(i)

                elif self.getMinimax() == Consts.minimaxType["minimize"]:
                    for i in xrange(lower_worst):
                        if self.internalPop.worstFitness(i).score[1] < newPop.worstFitness().score[1] and  self.internalPop.worstFitness(i).score[0] < newPop.worstFitness().score[0]:
                            newPop[-1] = self.internalPop.worstFitness(i)

            # check the half worst population of current generation
            # and bring them to replace the worst of next generation
            else:
                if self.getMinimax() == Consts.minimaxType["maximize"]:
                    for i in xrange(lower_worst):
                        if self.internalPop.worstRaw(i).score[0] > newPop.worstRaw().score[0]:
                            newPop[-1] = self.internalPop.worstRaw(i)
                elif self.getMinimax() == Consts.minimaxType["minimize"]:
                    for i in xrange(lower_worst):
                        if self.internalPop.worstRaw(i).score[0] < newPop.worstRaw().score[0]:
                            newPop[-1] = self.internalPop.worstRaw(i)

            logging.debug("FastConv time %f" % (time() - t))

        if self.elitism:
            t = time()
            logging.debug("Doing elitism.")
            if self.getMinimax() == Consts.minimaxType["maximize"]:
                for i in xrange(self.nElitismReplacement):
                    #re-evaluate before being sure this is the best
                    #self.internalPop.bestRaw(i).evaluate()
                    if self.internalPop.bestFitness(i).fitness > newPop.bestFitness(i).fitness:
                        newPop[len(newPop) - 1 - i] = self.internalPop.bestFitness(i)
                        # this force sorting again
            elif self.getMinimax() == Consts.minimaxType["minimize"]:
                for i in xrange(self.nElitismReplacement):
                    # re-evaluate before being sure this is the best
                    # self.internalPop.bestRaw(i).evaluate()
                    # not need since score did not changed
                    if self.internalPop.bestFitness(i).fitness < newPop.bestFitness(i).fitness:
                        newPop[len(newPop) - 1 - i] = self.internalPop.bestFitness(i)
            logging.debug("Elitism time %f" % (time() - t))
        self.internalPop = newPop
        self.internalPop.sort()

        logging.debug("The generation %d was finished.", self.currentGeneration)

        self.currentGeneration += 1

        if self.max_time:
           total_time = time() - self.time_init
           if total_time > self.max_time:
              return True
        return self.currentGeneration == self.nGenerations

    def printStats(self):
        """ Print generation statistics

        :rtype: the printed statistics as string

        .. versionchanged:: 0.6
           The return of *printStats* method.
        """
        percent = self.currentGeneration * 100 / float(self.nGenerations)
        message = "Gen. %d (%.2f%%):" % (self.currentGeneration, percent)
        logging.info(message)
        print message,
        sys_stdout.flush()
        self.internalPop.statistics()
        stat_ret = self.internalPop.printStats()
        return message + stat_ret

    def printTimeElapsed(self, msg=""):
        """ Shows the time elapsed since the begin of evolution """
        total_time = time() - self.time_init
        if msg:
            msg = "for " + msg
        print "Total time elapsed %s : %.3f seconds." % (msg, total_time)
        return total_time


    def evolve(self, freq_stats=0):
        """ Do all the generations until the termination criteria, accepts
        the freq_stats (default is 0) to dump statistics at n-generation

        Example:
           >>> ga_engine.evolve(freq_stats=10)
           (...)

        :param freq_stats: if greater than 0, the statistics will be
                           printed every freq_stats generation.
        :rtype: returns the best individual of the evolution

        .. versionadded:: 0.6
           the return of the best individual

        """

        stopFlagCallback = False
        stopFlagTerminationCriteria = False

        self.time_init = time()
        
        self.initialize()
        # self.printTimeElapsed("Initialize")
        self.internalPop.evaluate()
        # self.printTimeElapsed("Evaluate")
        self.internalPop.sort()
        # self.printTimeElapsed("Sort")
        logging.debug("Starting loop over evolutionary algorithm.")
        try:
            while True:
                if not self.stepCallback.isEmpty():
                    for it in self.stepCallback.applyFunctions(self):
                        stopFlagCallback = it

                if not self.terminationCriteria.isEmpty():
                    for it in self.terminationCriteria.applyFunctions(self):
                        stopFlagTerminationCriteria = it

                if freq_stats:
                    if (self.currentGeneration % freq_stats == 0) or (self.getCurrentGeneration() == 0):
                        self.printStats()
                        self.printTimeElapsed("Generation")

                if stopFlagTerminationCriteria:
                    logging.debug("Evolution stopped by the Termination Criteria !")
                    if freq_stats:
                        print "\n\tEvolution stopped by Termination Criteria function !\n"
                        print self.terminationCriteria
                    break

                if stopFlagCallback:
                    logging.debug("Evolution stopped by Step Callback function !")
                    if freq_stats:
                        print "\n\tEvolution stopped by Step Callback function !\n"
                    break

                if self.step():
                    break

        except KeyboardInterrupt:
            logging.debug("CTRL-C detected, finishing evolution.")
            if freq_stats:
                print "\n\tA break was detected, you have interrupted the evolution !\n"

        if freq_stats != 0:
            self.printStats()
            self.printTimeElapsed("Generation")

        return self.bestIndividual()

    def select(self, **args):
        """ Select one individual from population

        :param args: this parameters will be sent to the selector

        """
        for it in self.selector.applyFunctions(self.internalPop, **args):
            return it
