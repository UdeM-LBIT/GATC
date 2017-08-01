import copy
from ..TreeLib import TreeClass, TreeUtils, TreeFun, memorize
from ..reclkl import  computeQe, get_discr_size, computeProb, nodeLimitter
from scipy.stats import gamma
import numpy as np
import random

class EdgeParams:
    ISFixed = {'mu':False, 'sigma':False}
    rateLimit = [1e-6, 10, 1.0]

    def __init__(self, edgerates, ratelim=None):
        self.k = 0
        self.theta = 0
        self.mu = self._getparam(edgerates[0], 'mu')
        self.sigma = self._getparam(edgerates[1], 'sigma')        
        if ratelim and len(ratelim) == 3:
            EdgeParams.rateLimit = ratelim
        self.update()

    def update(self):
        self.k = (self.mu/self.sigma)**2
        self.theta = (self.sigma**2)/self.mu

    def _getparam(self, param, name):
        if isinstance(param, basestring):
            if param[-1] in ['f', 'F']:
                EdgeParams.ISFixed[name] =  True
                return float(param[:-1])
            else:
                EdgeParams.ISFixed[name] = False
        return float(param)

    def get_shape(self):
        return self.k

    def get_scale(self):
        return self.theta

    def clone(self):
        other = EdgeParams((self.mu, self.sigma))
        return other

    def is_mutable(self):
        return not all(EdgeParams.ISFixed.values())

    def mutate(self):
        other = self.clone()
        poss_choice = [k for k in EdgeParams.ISFixed.keys() if not EdgeParams.ISFixed[k]]
        if len(poss_choice) > 0:
            param = random.choice(poss_choice)
            setattr(other, param, min(EdgeParams.rateLimit[1], max(random.gauss(getattr(other, param), EdgeParams.rateLimit[-1]), EdgeParams.rateLimit[0])))
            other.update()
        return other

class DTLParams:
    ISFixed = {'dup':False, 'loss':False, 'trans':False}
    rateLimit = [1e-6, 10, 1.0]
    
    def __init__(self, dtl, parcim=True, ratelim=None):
        self.dup = self._getparam(dtl[0], 'dup')
        self.trans = self._getparam(dtl[1], 'trans')
        self.loss = self._getparam(dtl[-1], 'loss')
        self.parcim = parcim
        if self.parcim:
            for k, v in DTLParams.ISFixed.items():
                DTLParams.ISFixed[k] = True
        
        if ratelim and len(ratelim) == 3:
            DTLParams.rateLimit = ratelim

    def is_mutable(self):
        return not all(DTLParams.ISFixed.values())

    def clone(self):
        other = DTLParams((self.dup, self.trans, self.loss), parcim=self.parcim)
        return other

    def _getparam(self, param, name):
        if isinstance(param, basestring):
            if param[-1] in ['f', 'F']:
                DTLParams.ISFixed[name] = True
                return float(param[:-1])
            else: 
                DTLParams.ISFixed[name] = False
        return float(param)

    def getDup(self):
        return self.dup

    def getTrans(self):
        return self.trans

    def getLoss(self):    
        return self.loss

    def hastrans(self):
        return self.trans >0

    def getDTL(self):
        return (self.dup, self.trans, self.loss)

    def mutate(self):
        other = self.clone()
        poss_choice = [k for k in DTLParams.ISFixed.keys() if not DTLParams.ISFixed[k]]
        if len(poss_choice) > 0:
            param = random.choice(poss_choice)
            setattr(other, param, min(DTLParams.rateLimit[1], max(random.gauss(getattr(other, param), DTLParams.rateLimit[-1]), DTLParams.rateLimit[0])))
        return other

    def getHash(self):
        return hash(self.getDTL())


class ReconParams(object):
    """Little class to keep
    transfer parameter"""
    EVENT_LIST = ["SPR", "ROOT", "DTL", "EDGE"]
    def __init__(self, sptree, gtreesize, discrsize=10, parcim=False, stemlen=1.0, event_selector=[0.4, 0.2, 0.2, 0.2]):
        self.sptree = TreeClass(sptree)
        self.gtreesize = gtreesize
        self.discrsize =  discrsize
        self.stemlen = stemlen
        self.parcim = parcim
        self.data = {}
        self.default_event_selector = self._fixed_event_list(event_selector)
        self._spectree_preprocess()
          

    def _fixed_event_list(self, event_selector=None):
        if not event_selector:
            event_selector = self.default_event_selector
        tot_event = sum(event_selector)
        tot_event_at_zero = len([x for x in event_selector if x==0.0])
        diffprob = (tot_event - 1.0)/(4-tot_event_at_zero)
        return [(x - diffprob if x else x) for x in event_selector]


    def select_event(self, genome):
        if not genome.dtlrates.is_mutable():
            self.default_event_selector[2]=0.0
        if not genome.erates.is_mutable():
            self.default_event_selector[3]=0.0
        self.default_event_selector = self._fixed_event_list()
        choice = np.random.choice(ReconParams.EVENT_LIST, 1, p=self.default_event_selector)
        return choice

    def _spectree_preprocess(self):
        if self.parcim:
            self.sptree.label_internal_node()
            TreeUtils.lcaPreprocess(self.sptree)
        else:
            if not self.sptree.is_ultrametric():
                TreeFun.make_clock_like(self.sptree)

            self.sptree.compute_branches_length()            
            TreeFun.time_slice(self.sptree, timeframes=self.sptree.traverse(), single_nodes=True)
            TreeFun.scale_tree_height(self.sptree) # scale height to one
            slicelist, node_data = TreeFun.list_time_slice(self.sptree)
            self.data['slicelist'] = dict(slicelist)
            self.data['node_d'] = node_data
            self.data['leafslice'] = max(self.data['slicelist'].keys())
            self.discrsize = get_discr_size(self.gtreesize, self.discrsize, self.data['leafslice'])
            self.sptree.add_features(name2node=dict((x.name, x) for x in self.sptree.get_leaves()))

    def computeRecCost(self, gind, **kwargs):
        # get Qef ==> memorize or not
        #mu, sigma =  edgeparams.get_mu(), edgeparams.get_sigma()
        #k = (mu/sigma)**2
        #theta = (sigma**2)/mu
        dtlparams = gind.dtlrates
        gind.set_species()
        if self.parcim:
            if dtlparams.hastrans():
                return TreeUtils.computeDTLScore(gind.tree, self.sptree, dtlparams.getDup(), dtlparams.getTrans(), dtlparams.getLoss())
            else:
                lcamap = TreeUtils.lcaMapping(gind.tree, self.sptree)
                score = TreeUtils.computeDLScore(gind.tree, lcamap, dtlparams.getDup(), dtlparams.getLoss()) 
                return sum(score)
        else:
            k = 1.0
            theta = 1.0
            Qef = computeMat(dtlparams.getHash(), self, dtlparams)

            edgeparams = gind.erates#kwargs.get('edgeparams', None)
            if edgeparams:
                k = edgeparams.get_shape()
                theta = edgeparams.get_scale()
            ratedens = gamma(k, scale=theta)
            nodeLimitter(gind.tree, self.discrsize, self.data['leafslice'])
            prob_Ax = computeProb(ratedens, gind.tree, self.sptree.name2node, self.data['slicelist'], self.data['node_d'], self.discrsize, self.stemlen, dtlparams.getDup(), dtlparams.getTrans(), Qef)
            val = np.asarray(prob_Ax)[gind.tree.ind, self.sptree.edge_i, self.discrsize-1]
            #print val, -np.log(val)
            return -np.log(val)


@memorize
def computeMat(rec, dtlparams):
    return computeQe(rec.sptree, rec.data['slicelist'], rec.data['node_d'], rec.discrsize, dtlparams.getDup(), dtlparams.getTrans(), dtlparams.getLoss(), rec.stemlen)



