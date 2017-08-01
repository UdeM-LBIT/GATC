"""

:mod:`GenomeBase` -- the genomes base module
================================================================

This module have the class which every representation extends,
if you are planning to create a new representation, you must
take a inside look into this module.

"""
from FunctionSlot import FunctionSlot

class GenomeBase(object):
    """ GenomeBase Class - The base of all chromosome representation """
    __slots__ = ["evaluator", "initializator", "mutator", "crossover", "internalParams", "score", "fitness"]

    def __init__(self):
        """Genome Constructor"""
        self.evaluator = FunctionSlot("Evaluator")
        self.initializator = FunctionSlot("Initializator")
        self.mutator = FunctionSlot("Mutator")
        self.crossover = FunctionSlot("Crossover")

        self.internalParams = {}
        self.score = []
        self.fitness = 0.0

    def getRawScore(self):
        """ Get the Raw Score of the genome

        :rtype: genome raw score, which is a list

        """
        return self.score

    def getFitnessScore(self):
        """ Get the Fitness Score of the genome

        :rtype: genome fitness score

        """
        return self.fitness

    def __repr__(self):
        """String representation of Genome"""
        allSlots = [self.evaluator, self.initializator, self.mutator,
                        self.crossover]

        ret = "- GenomeBase\n"
        ret += "\tScore:\t\t\t %s\n" % (", ".join(["%.5f"%x for x in  self.score]),)
        ret += "\tFitness:\t\t %.5f\n\n" % (self.fitness,)
        ret += "\tParams:\t\t %s\n\n" % (self.internalParams,)

        for slot in allSlots:
            ret += "\t" + slot.__repr__()
        ret += "\n"

        return ret

    def setParams(self, **args):
        """ Set the internal params

        Example:
            >>> genome.setParams(rangemin=0, rangemax=100, gauss_mu=0, gauss_sigma=1)

        .. note:: All the individuals of the population shares this parameters and uses
                     the same instance of this dict.

        :param args: this params will saved in every chromosome for genetic op. use

        """
        self.internalParams.update(args)


    def getParam(self, key, nvl=None):
        """ Gets an internal parameter

        Example:
            >>> genome.getParam("rangemax")
            100

        .. note:: All the individuals of the population shares this parameters and uses
                     the same instance of this dict.

        :param key: the key of param
        :param nvl: if the key doesn't exist, the nvl will be returned

        """
        return self.internalParams.get(key, nvl)

    def resetStats(self):
        """ Clear score and fitness of genome """
        self.score = []
        self.fitness = 0

    def evaluate(self, **args):
        """ Called to evaluate genome

        :param args: this parameters will be passes to the evaluator

        """
        self.resetStats()
        for it in self.evaluator.applyFunctions(self, **args):
            self.score.append(it)

    def initialize(self, **args):
        """ Called to initialize genome

        :param args: this parameters will be passed to the initializator

        """
        for it in self.initializator.applyFunctions(self, **args):
            pass

    def mutate(self, **args):
        """ Called to mutate the genome

        :param args: this parameters will be passed to the mutator
        :rtype: the number of mutations returned by mutation operator

        """
        nmuts = 0
        for it in self.mutator.applyFunctions(self, **args):
            nmuts += it
        return nmuts

    def copy(self, g):
        """ Copy the current GenomeBase to 'g'

        :param g: the destination genome

        .. note:: If you are planning to create a new chromosome representation, you
                     **must** implement this method on your class.

        """
        g.score = self.score
        g.fitness = self.fitness
        g.evaluator = self.evaluator
        g.initializator = self.initializator
        g.mutator = self.mutator
        g.crossover = self.crossover
        g.internalParams = self.internalParams

    def clone(self):
        """ Clone this GenomeBase

        :rtype: the clone genome

        .. note:: If you are planning to create a new chromosome representation, you
                     **must** implement this method on your class.
        """
        newcopy = GenomeBase()
        self.copy(newcopy)
        return newcopy
