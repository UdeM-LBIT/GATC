"""
:mod:`GPopulation` -- the population module
================================================================

This module contains the :class:`GPopulation.GPopulation` class, which is reponsible
to keep the population and the statistics.

Default Parameters
-------------------------------------------------------------

*Minimax*

    >>> Consts.minimaxType["maximize"]

    Maximize the evaluation function

*Scale Method*

    :func:`Scaling.LinearScaling`

    The Linear Scaling scheme

Class
-------------------------------------------------------------


"""

import Consts
import Util
import numpy as np
from .FunctionSlot import FunctionSlot
from .Statistics import Statistics
import logging
from functools import partial

try:
    import pathos.multiprocessing as mp
    CPU_COUNT = mp.cpu_count()
    MULTI_PROCESSING = True if CPU_COUNT > 1 else False
    logging.debug("You have %d CPU cores, so the multiprocessing state is %s", CPU_COUNT, MULTI_PROCESSING)
except ImportError:
    MULTI_PROCESSING = False
    logging.debug("You don't have multiprocessing support for your Python version !")


def key_raw_score(individual, k=0):
    """ A key function to return raw score

    :param individual: the individual instance
    :rtype: the individual raw score

    .. note:: this function is used by the max()/min() python functions

    """
    return individual.score[k]

def key_fitness_score(individual):
    """ A key function to return fitness score, used by max()/min()

    :param individual: the individual instance
    :rtype: the individual fitness score

    .. note:: this function is used by the max()/min() python functions

    """
    return individual.fitness


def multiprocessing_eval(ind, args):
    """ Internal used by the multiprocessing """
    pos, ind = ind
    args['ext']=str(pos)
    ind.evaluate(**args)
    return ind.score

def multiprocessing_eval_full(ind, args):
    """ Internal used by the multiprocessing (full copy)"""
    pos, ind = ind
    args['ext']=str(pos)
    ind.evaluate(**args)
    return ind


class GPopulation(object):
    """ GPopulation Class - The container for the population

    **Examples**
     Get the population from the :class:`GSimpleGA.GSimpleGA` (GA Engine) instance
         >>> pop = ga_engine.getPopulation()

     Get the best fitness individual
         >>> bestIndividual = pop.bestFitness()

     Get the best raw individual
         >>> bestIndividual = pop.bestRaw()

     Get the statistics from the :class:`Statistics.Statistics` instance
         >>> stats = pop.getStatistics()
         >>> print stats["rawMax"]
         10.4

     Iterate, get/set individuals
         >>> for ind in pop:
         >>>   print ind
         (...)

         >>> for i in xrange(len(pop)):
         >>>    print pop[i]
         (...)

         >>> pop[10] = newGenome
         >>> pop[10].fitness
         12.5

    :param genome: the :term:`Sample genome`, or a GPopulation object, when cloning.

    """

    def __init__(self, genome, single=True, bulkEval=False):
        """ The GPopulation Class creator """

        if isinstance(genome, GPopulation):
            self.oneSelfGenome = genome.oneSelfGenome
            self.blkevaluator = genome.blkevaluator
            self.bulkEval = genome.bulkEval
            self.internalPop = []
            self.internalPopRaw = []
            self.popSize = genome.popSize
            self.sorted = False
            self.minimax = genome.minimax
            self.scaleMethod = genome.scaleMethod
            self.rawSortMethod =  genome.rawSortMethod
            self.allSlots = [self.scaleMethod, self.rawSortMethod]
            self.bulkEval = genome.bulkEval
            self.internalParams = genome.internalParams
            self.multiProcessing = genome.multiProcessing
            self.statted = False
            self.stats = Statistics()
            return

        logging.debug("New population instance, %s class genomes.", genome.__class__.__name__)
        self.blkevaluator = FunctionSlot("Whole pop evaluator")
        self.internalPop = []
        self.internalPopRaw = []
        self.popSize = 0
        self.sorted = False
        self.minimax = Consts.CDefPopMinimax
        self.scaleMethod = FunctionSlot("Scale Method")
        self.scaleMethod.set(Consts.CDefPopScale)
        self.rawSortMethod = FunctionSlot("Raw Sort Method")
        self.rawSortMethod.set(Consts.CDefRawPopSort)
        self.allSlots = [self.scaleMethod, self.rawSortMethod]

        self.internalParams = {}
        self.multiProcessing = (False, False, None)

        # Bulk evaluation
        self.bulkEval = bulkEval

        # Statistics
        self.statted = False
        self.stats = Statistics()

        self.oneSelfGenome = genome if not isinstance(genome, list) else np.random.choice(genome)

        if isinstance(genome, list) and not single:
            self.internalPop = genome

    def setBulkEval(self, value):
        """ Set the bulk evaluation value """
        self.bulkEval = value

    def setPopulationEvaluator(self, fn):
        """Use input function to set whole population evaluator"""
        self.blkevaluator.set(fn)


    def setScaleMethod(self, fn):
        """Setting scale method for the GA"""
        self.scaleMethod.set(fn)


    def getOrSetWeight(self, lklscaled=None):
        # current value
        weight = self.getParam('weight', None)
        # if this value is not found
        if not weight:
            try:
                costdict = self.getStatistics(1) 
                ave_cost = costdict["rawAve"]
                if not lklscaled:
                    lkldict = self.getStatistics(0)
                    ave_lkl = lkldict["rawAve"]
                else:
                    ave_lkl = np.median(lklscaled) # should I take the mean instead ?
                # get the closest power of 10 and take its inverse
                weight = [1/(10**(round(np.log10(ave_lkl)) -1)), 1/(10**round(np.log10(ave_cost)))]
            except Exception as e:
                # fixed default value
                weight = [0.001, 0.1]
            self.setParams(weight=weight)
        return weight
    
    def setRawSortMethod(self, fn):
        """Setting scale method for the GA"""
        self.rawSortMethod.set(fn)


    def setMultiProcessing(self, flag=True, full_copy=False, max_processes=None):
        """ Sets the flag to enable/disable the use of python multiprocessing module.
         Use this option when you have more than one core on your CPU and when your
         evaluation function is very slow.
         The parameter "full_copy" defines where the individual data should be copied back
         after the evaluation or not. This parameter is useful when you change the
         individual in the evaluation function.

         :param flag: True (default) or False
         :param full_copy: True or False (default)
         :param max_processes: None (default) or an integer value
 
         .. warning:: Use this option only when your evaluation function is slow, se you
                 will get a good tradeoff between the process communication speed and the
                 parallel evaluation.

         The `setMultiProcessing` method.

        """
        self.multiProcessing = (flag, full_copy, max_processes)

    def setMinimax(self, minimax):
        """ Sets the population minimax

        Example:
         >>> pop.setMinimax(Consts.minimaxType["maximize"])

        :param minimax: the minimax type

        """
        self.minimax = minimax

    def __repr__(self):
        """ Returns the string representation of the population """
        ret = "- GPopulation\n"
        ret += "\tPopulation Size:\t %d\n" % (self.popSize,)
        ret += "\tMinimax Type:\t\t %s\n" % (Consts.minimaxType.keys()[Consts.minimaxType.values().index(self.minimax)].capitalize(),)
        for slot in self.allSlots:
            ret += "\t" + slot.__repr__()
        ret += "\n"
        ret += self.stats.__repr__()
        return ret

    def __len__(self):
        """ Return the length of population """
        return len(self.internalPop)

    def __getitem__(self, key):
        """ Returns the specified individual from population """
        return self.internalPop[key]

    def __iter__(self):
        """ Returns the iterator of the population """
        return iter(self.internalPop)

    def __setitem__(self, key, value):
        """ Set an individual of population """
        self.internalPop[key] = value
        self.clearFlags()

    def clearFlags(self):
        """ Clear the sorted and statted internal flags """
        self.sorted = False
        self.statted = False

    def getStatistics(self, ftfunc=None):
        """ Return a Statistics class for statistics

        :rtype: the :class:`Statistics.Statistics` instance

        """
        self.statistics()
        if ftfunc is None:
            return self.stats
        if not isinstance(ftfunc, int):
            Util.raiseException("ftfunc of getStatistics should be an integer")
        return self.stats.getRawScore(ftfunc)

    def statistics(self):
        """ Do statistical analysis of population and set 'statted' to True"""
        if self.statted:
            return
        logging.debug("Running statistical calculations")
        nscore = len(self[0].score)
        len_pop = len(self)

        for ftfunc in xrange(nscore):
            raw_tot = []

            for ind in xrange(len_pop):
                raw_tot.append(self[ind].score[ftfunc])

            keyfn = partial(key_raw_score, k=ftfunc)
            self.stats["rawMax"].append(max(self, key=keyfn).score[ftfunc])
            self.stats["rawMin"].append(min(self, key=keyfn).score[ftfunc])
            self.stats["rawAve"].append(np.mean(raw_tot))
            self.stats["rawMed"].append(np.median(raw_tot))

            tmpvar = 0.0
            for ind in xrange(len_pop):
                s = self[ind].score[ftfunc] - self.stats["rawAve"][ftfunc]
                s *= s
                tmpvar += s

            tmpvar /= float((len(self) - 1))
            try:
                self.stats["rawDev"].append(np.sqrt(tmpvar))
            except:
                self.stats["rawDev"].append(0.0)

            self.stats["rawVar"].append(tmpvar)

        self.statted = True

    def bestFitness(self, index=0):
        """ Return the best fitness individual of population

        :param index: the *index* best individual
        :rtype: the individual

        """
        if not self.sorted:
            self.sort()
        return self.internalPop[index]

    def worstFitness(self, index=0):
        """ Return the worst fitness individual of the population

        :rtype: the individual

        """
        if not self.sorted:
            self.sort()
        return self.internalPop[-index-1]

    def bestRaw(self, index=0, score=0):
        """ Return the best raw score individual of population

        :param index: the *index* best raw individual
        :rtype: the individual

        """
        if not self.sorted:
            self.sort(raw_score=score)
        return self.internalPopRaw[index]

    
    def worstRaw(self, index=0, score=0):
        """ Return the worst raw score individual of population

        :rtype: the individual

        """
        if not self.sorted:
            self.sort(raw_score=score)
        return self.internalPopRaw[-index-1]


    def sort(self, raw_score=0):
        """ Sort the population """
        if self.sorted:
            return
        rev = (self.minimax == Consts.minimaxType["maximize"])

        for it in self.rawSortMethod.applyFunctions(self, reverse=rev, score=raw_score):
            pass
        
        # copy current order into rawSort
        self.internalPopRaw = self.internalPop[:]
        # now scale and sort by fitness
        self.scale()
        self.internalPop.sort(cmp=Util.cmp_individual_fitness, reverse=rev)

        self.sorted = True


    def setPopulationSize(self, size):
        """ Set the population size

        :param size: the population size

        """
        self.popSize = size


    def create(self, **args):
        """ Clone the example genome to fill the population """
        self.minimax = args["minimax"]
        # not defined
        if not self.internalPop:
            self.internalPop = [self.oneSelfGenome.clone() for i in xrange(self.popSize)]
        # if larger than wanted popsize
        elif len(self.internalPop) > self.popSize:
            self.internalPop = np.random.choice(self.internalPop, self.popSize, replace=False)
        # fill with a clone of some genome
        # here we enable the repetition of some of them
        else:
            missing = self.popSize - len(self.internalPop)
            self.internalPop.extend([x.clone() for x in np.random.choice(self.internalPop, missing)])
        self.internalPop = list(self.internalPop)
        self.clearFlags()

    def __findIndividual(self, individual, end):
        for i in xrange(end):
            if individual.compare(self.internalPop[i]) == 0:
                return True

    def initialize(self, **args):
        """ Initialize all individuals of population,
        this calls the initialize() of individuals """
        logging.debug("Initializing the population")
        limit =  args.get('limit', 10)
        if self.oneSelfGenome.getParam("full_diversity", True) and hasattr(self.oneSelfGenome, "compare"):
            for i in xrange(len(self.internalPop)):
                curr = self.internalPop[i]
                curr.initialize(**args)
            while self.__findIndividual(curr, i):
                curr.initialize(**args)
                # infinite loop breaker
                limit -= 1
                if limit < 1:
                    break
        else:
            for gen in self.internalPop:
                gen.initialize(**args)
        self.clearFlags()

    def evaluate(self, **args):
        """ Evaluate all individuals in population, calls the evaluate() method of individuals

        :param args: this params are passed to the evaluation function

        """
        # We have multiprocessing
        args.update(self.internalParams)
        m_e_f = partial(multiprocessing_eval_full, args=args)
        m_e = partial(multiprocessing_eval, args=args)
        if self.multiProcessing[0] and MULTI_PROCESSING:
            # print("Multiprocessing evaluation chosen")
            logging.debug("Evaluating the population using the multiprocessing method")
            proc_pool = mp.Pool(processes=self.multiProcessing[2])

        # Multiprocessing full_copy parameter
            if self.multiProcessing[1]:
                results = proc_pool.map(m_e_f, enumerate(self.internalPop))
                proc_pool.close()
                proc_pool.join()
                for i in xrange(len(self.internalPop)):
                    self.internalPop[i] = results[i]
            else:
                results = proc_pool.map(m_e, enumerate(self.internalPop))
                proc_pool.close()
                proc_pool.join()
                for individual, score in zip(self.internalPop, results):
                    individual.score = score
        
        elif self.bulkEval and not self.blkevaluator.isEmpty():
            # print("*** Bulk evaluate chosen")
            logging.debug("Evaluating the population using bulk evaluator")

            scores = np.zeros(len(self.internalPop))
            for it in self.blkevaluator.applyFunctions(self, **args):
                scores += np.asarray(it)
            for pos, ind in enumerate(self.internalPop):
                ind.resetStats()
                ind.score.append(scores[pos])

            # look at recon cost function now
            if len(ind.evaluator) > 1:
                for ind in self.internalPop:
                    reccost = ind.evaluator.apply(-1, ind, **args)
                    ind.score.append(reccost)
        
        else:
            print("*** Single evaluate chosen")
            for ind in self.internalPop:
                ind.evaluate(**args)
        self.clearFlags()



    def scale(self, **args):
        """ Scale the population using the scaling method

        :param args: this parameter is passed to the scale method

        """
        for it in self.scaleMethod.applyFunctions(self, **args):
            pass

        fit_sum = 0
        for ind in xrange(len(self)):
            fit_sum += self[ind].fitness

        self.stats["fitMax"] = max(self, key=key_fitness_score).fitness
        self.stats["fitMin"] = min(self, key=key_fitness_score).fitness
        self.stats["fitAve"] = fit_sum / float(len(self))

        self.sorted = False

    def printStats(self):
        """ Print statistics of the current population """
        message = "Max/Min/Avg Fitness [%(fitMax).2f/%(fitMin).2f/%(fitAve).2f]" % self.stats
        print self.stats
        logging.info(message)
        return message

    def copy(self, pop):
        """ Copy current population to 'pop'
        :param pop: the destination population
         .. warning:: this method do not copy the individuals, only the population logic

        """
        pop.popSize = self.popSize
        pop.minimax = self.minimax
        pop.scaleMethod = self.scaleMethod
        pop.internalParams = self.internalParams
        pop.multiProcessing = self.multiProcessing
        pop.bulkEval = self.bulkEval


    def getParam(self, key, nvl=None):
        """ Gets an internal parameter
        Example:
        >>> population.getParam("tournamentPool")
        5
        :param key: the key of param
        :param nvl: if the key doesn't exist, the nvl will be returned
        """
        return self.internalParams.get(key, nvl)

    def setParams(self, **args):
        """ Gets an internal parameter

        Example:
            >>> population.setParams(tournamentPool=5)

        :param args: parameters to set
        """
        self.internalParams.update(args)

    def clear(self):
        """ Remove all individuals from population """
        del self.internalPop[:]
        del self.internalPopRaw[:]
        self.clearFlags()

    def clone(self):
        """ Return a brand-new cloned population """
        newpop = GPopulation(self.oneSelfGenome)
        self.copy(newpop)
        return newpop
