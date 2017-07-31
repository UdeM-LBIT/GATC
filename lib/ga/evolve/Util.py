"""

:mod:`Util` -- utility module
============================================================================

This is the utility module, with some utility functions of general
use, like list item swap, random utilities and etc.

"""

from random import random as rand_random
from math import sqrt as math_sqrt
import logging
import Consts

def compNextGen(newPop, oldPop, requiredSize):
    Pstar = list(newPop + oldPop)
    Waves = []
    prevlen = len(Pstar)
    wave_ind = 1
    dranklist = __computeDominanceRank(Pstar)

    while(True):
        wave_i = []
        keeplist = []
        for i, p in enumerate(Pstar):
            if dranklist[i]==0:
                p.fitness = wave_ind
                wave_i.append(p)
            else:
                keeplist.append(p)
        if len(wave_i) == 0:
            break
        Waves.extend(wave_i)
        Pstar = keeplist
        
        newlen = len(Pstar)
        assert(prevlen > newlen) , "Wrong"
        newlen = prevlen
        
        dranklist = __computeDominanceRank(Pstar)
        wave_ind += 1
    
    for i, p in enumerate(Pstar):
        p.fitness = wave_ind + dranklist[i]
        Waves.append(p)

    Waves.sort(key=lambda x: x.score[-1])
    Waves.sort(cmp=cmp_individual_fitness)

    nextPop = Waves[:requiredSize]
    return nextPop

def __computeDominanceRank(pop):
    dranklist = []
    for i, pi in enumerate(pop):
        drank = 0
        for pj in pop:
            if (pj.score[0] < pi.score[0]) and (pj.score[1] < pi.score[1]):
                drank += 1
        dranklist.append(drank)
    return dranklist


def randomFlipCoin(p):
    """Returns True with the *p* probability. If *p* is 1, the
    function will always return True. If *p* is 0, the function will
    return always False.

    Example:
       >>> Util.randomFlipCoin(1.0)
       True

    :param p: probability, between 0.0 and 1.0
    :rtype: True or False

    """
    if p == 1.0:
        return True
    if p == 0.0:
        return False

    return rand_random() <= p



def raiseException(message, expt=None):
    """ Raise an exception and logs the message.

    Example:
       >>> Util.raiseException('The value is not an integer', ValueError)

    :param message: the message of exception
    :param expt: the exception class
    :rtype: None

    """
    logging.critical(message)
    if expt is None:
        raise Exception(message)
    else:
        raise (expt, message)


def cmp_individual_fitness(a, b):
    """ Compares two individual fitness scores, used for sorting population

    Example:
       >>> GPopulation.cmp_individual_scaled(a, b)

    :param a: the A individual instance
    :param b: the B individual instance
    :rtype: 0 if the two individuals fitness score are the same,
            -1 if the B individual fitness score is greater than A and
            1 if the A individual fitness score is greater than B.

    .. note:: this function is used to sorte the population individuals

    """
    if a.fitness < b.fitness:
        return -1
    if a.fitness > b.fitness:
        return 1
    return 0


def RawSorting(pop, **args):
    """Will sort using the score as key"""
    pop.internalPop.sort(key=lambda x: x.score[args.get('score', 0)], reverse=args.get('reverse', False))


class ErrorAccumulator(object):
    """ An accumulator for the Root Mean Square Error (RMSE) and the
    Mean Square Error (MSE)
    """
    def __init__(self):
        self.acc = 0.0
        self.acc_square = 0.0
        self.acc_len = 0

    def reset(self):
        """ Reset the accumulator """
        self.acc_square = 0.0
        self.acc = 0.0
        self.acc_len = 0

    def append(self, target, evaluated):
        """ Add value to the accumulator

        :param target: the target value
        :param evaluated: the evaluated value
        """
        self.acc_square += (target - evaluated) ** 2
        self.acc += abs(target - evaluated)
        self.acc_len += 1

    def __iadd__(self, value):
        """ The same as append, but you must pass a tuple """
        self.append(*value)
        return self

    def getMean(self):
        """ Return the mean of the non-squared accumulator """
        return self.acc / self.acc_len

    def getSquared(self):
        """ Returns the squared accumulator """
        return self.acc_square

    def getNonSquared(self):
        """ Returns the non-squared accumulator """
        return self.acc

    def getAdjusted(self):
        """ Returns the adjusted fitness
        This fitness is calculated as 1 / (1 + standardized fitness)
        """
        return 1.0 / (1.0 + self.acc)

    def getRMSE(self):
        """ Return the root mean square error

        :rtype: float RMSE
        """
        return math_sqrt(self.acc_square / float(self.acc_len))

    def getMSE(self):
        """ Return the mean square error

        :rtype: float MSE
        """
        return self.acc_square / float(self.acc_len)


class Graph(object):
    """ The Graph class

    Example:
       >>> g = Graph()
       >>> g.addEdge("a", "b")
       >>> g.addEdge("b", "c")
       >>> for node in g:
       ...    print node
       a
       b
       c

    .. versionadded:: 0.6
       The *Graph* class.
    """

    def __init__(self):
        """ The constructor """
        self.adjacent = {}

    def __iter__(self):
        """ Returns an iterator to the all graph elements """
        return iter(self.adjacent)

    def addNode(self, node):
        """ Add the node

        :param node: the node to add
        """
        if node not in self.adjacent:
            self.adjacent[node] = {}

    def __iadd__(self, node):
        """ Add a node using the += operator """
        self.addNode(node)
        return self

    def addEdge(self, a, b):
        """ Add an edge between two nodes, if the nodes
        doesn't exists, they will be created

        :param a: the first node
        :param b: the second node
        """
        if a not in self.adjacent:
            self.adjacent[a] = {}

        if b not in self.adjacent:
            self.adjacent[b] = {}

        self.adjacent[a][b] = True
        self.adjacent[b][a] = True

    def getNodes(self):
        """ Returns all the current nodes on the graph

        :rtype: the list of nodes
        """
        return self.adjacent.keys()

    def reset(self):
        """ Deletes all nodes of the graph """
        self.adjacent.clear()

    def getNeighbors(self, node):
        """ Returns the neighbors of the node

        :param node: the node
        """
        return self.adjacent[node].keys()

    def __getitem__(self, node):
        """ Returns the adjacent nodes of the node """
        return self.adjacent[node].keys()

    def __repr__(self):
        ret = "- Graph\n"
        ret += "\tNode list:\n"
        for node in self:
            ret += "\t\tNode [%s] = %s\n" % (node, self.getNeighbors(node))
        return ret

