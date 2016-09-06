from pyevolve.GenomeBase import GenomeBase
from collections import Counter, defaultdict as ddict
import numpy as np
import hashlib
from ..TreeLib import TreeClass


class BranchRep:
    
    def __init__(self,part1, part2, btype):
        self.part1 =  part1
        self.part2 = part2
        self.btype =  btype
    
    def get_hash(self, partition):
        return hashlib.sha384(",".join(sorted(partition))).hexdigest()

    def link_leaf(self):
        return self.btype == 0
        
    def __hash__(self):
        return hash((self.get_hash(self.part1), self.get_hash(self.part2)))
    
    def __eq__(self, other):
        p1 = self.get_hash(self.part1) == other.get_hash(other.part1)
        p2 = self.get_hash(self.part2) == other.get_hash(other.part2)
        return self.btype == other.btype and p1 and p2 
        
class Utils:
    
    @staticmethod
    def shuffle_map(gmap):
        for k,v in gmap.items():
            np.random.shuffle(v)
        return gmap

    @staticmethod
    def initialize(genome, **args):
        gmap = Utils.shuffle_map(args["gmap"])
        pos = {}
        for node in genome.tree:
            node.add_features(species=node.name)
            pos[node.species] = pos.get(node.species, 0)
            node.name =  gmap[node.species][pos[node.species]]
            pos[node.species] += 1
    

    @staticmethod
    def iter_partition(tree):
        find_root = False
        for node in tree.iter_descendants("levelorder"):
            part1 = set(node.get_leaf_names())
            part2 = set(tree.get_leaf_names()) - part1
            nodes = [node]
            importance = 1
            if(not node.up.is_root()):
                if node.is_leaf():
                    importance = 0
                yield part1, part2, nodes, importance
            elif not find_root:
                find_root = True
                nodes.extend(node.get_sisters())
                yield part1, part2, nodes, importance

            
    @staticmethod
    def two_step_branch_selection(candbranch, prob):
        btype = np.random.choice([1,0], p=[prob, 1-prob])
        new_cands = [x for x in candbranch if x.btype == btype]
        # if it's empty, then the list only contains abs(btype-1)
        if not new_cands:
            new_cands = candbranch     
        #print new_cands 
        #print '------------------------'
        #print candbranch
        #print '\n'
        selected_branch = np.random.choice(new_cands)
        return selected_branch
    
    @staticmethod
    def find_and_swap(branch, g1, g2):
        
        if not branch.link_leaf():
            node1 = g1.tree.get_common_ancestor(branch.part1)
            node2 = g1.tree.get_common_ancestor(branch.part2)
            g1swap,  g1part = (node1,  branch.part1) if node2.is_root() else (node2, branch.part2)
            g2swap = g2.tree.get_common_ancestor(g1part)
            
        else:
            leaf_node = list((branch.part1 if len(branch.part1) == 1 else branch.part2))[0]
            g1swap = g1.tree&leaf_node
            g2swap = g2.tree&leaf_node
                    
                
        g1_parent = g1swap.up
        g2_parent = g2swap.up
        g2_parent.add_child(g1swap.detach())
        g1_parent.add_child(g2swap.detach())
        return g1, g2
        
        
    @staticmethod
    def permute_seq(genome, spec):
        nlist = genome.tree.search_nodes(species=spec)
        if len(nlist)>1:
            n1, n2 = np.random.choice(nlist, 2, replace=False)
            n1.name, n2.name = n2.name, n1.name            
        
    
    @staticmethod
    def crossover(genome, **args):
        gdad = args['dad']
        gmom = args['mom']
        genome1 = gdad.clone()
        genome2 = gmom.clone()
        prob = max(genome1.intbrnp, genome2.intbrnp)
        current_dict = ddict(int)
        genomes = [genome1, genome2]
        for g in genomes:
            for part1, part2, nodes, imp in Utils.iter_partition(g.tree):
                brnch = BranchRep(part1, part2, imp)
                current_dict[brnch] += 1
        
        candidates = [k for k,v in current_dict.items() if v > 1]
        selected_branch = Utils.two_step_branch_selection(candidates, prob)
        return Utils.find_and_swap(selected_branch, genome1, genome2)
        
        
    @staticmethod
    def mutate(genome, **args):
        if args["pmut"] <= 0.0:
            return 0
        nspec = genome.get_spec_len()
        av_mutation =  args["pmut"]*nspec
        spec_list = genome.spcount.keys()
        if av_mutation < 1.0:
            av_mutation = 0
            for spec in spec_list:
                if np.random.rand() < args["pmut"]:
                    Utils.permute_seq(genome, spec)
                    av_mutation += 1
        else:
            for i in range(int(round(av_mutation))):
                # choose a random 
                spec = np.random.choice(spec_list)
                Utils.permute_seq(genome, spec)
                
            
        return av_mutation

    @staticmethod
    def evaluate(genome, **args):
        raxmlmodel = genome.model
        score = raxmlmodel.optimize_model(genome.tree)
        return score

    
class GPolySolver(GenomeBase):
    
    gmap = {}

    def __init__(self, tree, model, intbrnp=0.95, gmap={}):
        GenomeBase.__init__(self)
        self.tree =  tree
        self.model = model
        self.intbrnp = intbrnp
        if gmap:
            self.setGeneMap(gmap)
        self.spcount = None
        try:
            self.spcount =  Counter(tree.get_leaf_names())
        except:
            pass
        self.initializator.set(Utils.initialize)
        self.mutator.set(Utils.mutate)
        self.evaluator.set(Utils.evaluate)
        self.crossover.set(Utils.crossover)

    def get_spec_len(self):
        return len(self.spcount.keys())
    
    @classmethod
    def setGeneMap(clc, val):
        clc.gmap = val

    def initialize(self, **args):
        """ Called to initialize genome
        :param args: this parameters will be passed to the initializator
        """
        args['gmap'] = GPolySolver.gmap
        for it in self.initializator.applyFunctions(self, **args):
            pass

    def copy(self, g):
        """ Copy the current GenomeBase to 'g'
        :param g: the destination genome
        .. note:: If you are planning to create a new chromosome representation, you
                **must** implement this method on your class.
        """
        GenomeBase.copy(self, g)
        g.tree =  self.tree.copy()
        g.spcount = self.spcount
        g.intbrnp = self.intbrnp
        g.model = self.model
        
    def clone(self):
        """ Clone this GenomeBase
        :rtype: the clone genome
        .. note:: If you are planning to create a new chromosome representation, you
            **must** implement this method on your class.
        """
        newcopy = GPolySolver(None, None)
        self.copy(newcopy)
        return newcopy
