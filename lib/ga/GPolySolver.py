from pyevolve.GenomeBase import GenomeBase
from collections import Counter, defaultdict as ddict
import numpy as np
import hashlib
from itertools import permutations, product
from ..TreeLib import TreeClass
from heapq import heappushpop, heappush


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
    def _enumerate_permutation(gmap):
        perm_keeper =  ddict(list)
        for k,v in gmap.items():
            perm_keeper[k] = [x for x in permutations(v)]
        return perm_keeper
    
    @staticmethod
    def best_tree_finder(trees, raxmlmod, gmap, ntrees, timelimit=None):
        best_trees = []
        perm_keeper = Utils._enumerate_permutation(gmap)
        keylist =  dict((x,y) for y, x in enumerate (perm_keeper.keys()))
        product_list = product(*perm_keeper.values())
        i = 0 
        start_time = time.time()
        stop = False
        for plabel in product_list:
            if not stop:
                for t in trees:
                    ind = [0 for ckey in keylist]
                    tcop = t.copy()
                    for node in tcop:
                        spec_pos = keylist[node.name]
                        node.name =  plabel[spec_pos][ind[spec_pos]]
                        ind[spec_pos] += 1
                    score = raxmlmod.optimize_model(tcop)
                    i += 1
                    print("Trees %d, score : %f" %(i, score))
                    if len(best_trees) <= ntrees:
                        heappush(best_trees, (score, tcop))
                    else:
                        # here push the new item then pop the tree with 
                        # smallest value from it 
                        heappushpop(best_trees, (score, tcop))
                        if timelimit and cur_time - start_time > timelimit: 
                            stop = True
                            break
        return best_trees

    
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
            leafset1 = set(node.get_leaves())
            leafset2 = set(tree.get_leaves()) - part1
            part1 =  [x.species for x in leafset1]
            part2 =  [x.species for x in leafset2]
            importance = 1
            if(not node.up.is_root()):
                if node.is_leaf():
                    importance = 0
                yield part1, part2, leafset1, importance
            elif not find_root:
                find_root = True
                yield part1, part2, leafset1, importance

            
    @staticmethod
    def two_step_branch_selection(candbranch, prob):
        btype = np.random.choice([1,0], p=[prob, 1-prob])
        new_cands = [x for x in candbranch if x[0].btype == btype]
        # if it's empty, then the list only contains abs(btype-1)
        if not new_cands:
            new_cands = candbranch     

        selected_branch = np.random.choice(new_cands)
        return selected_branch
    
    @staticmethod
    def find_and_swap(br, g1, g2):
        branch, nodelist = br
        g1_list = []
        g2_list = []
        for x  in nodelist:
            if x[-1] == 0:
                g1_list.append(x[0])
            else:
                g2_list.append(x[0])
            
        g1swap = g1.tree.get_common_ancestor(np.random.choice(g1_list))
        g2swap = g2.tree.get_common_ancestor(np.random.choice(g2_list))
       
        g1_parent = g1swap.up
        g2_parent = g2swap.up

        # fix problem due to same name multiple time
        g1_swapper = ddict(list)
        g2_swapper = ddict(list)
        not_in_g1_swapper = set(g1.tree) - set(g1_swapper)
        for n in not_in_g1_swapper:
            g1_swapper[n.species].append(n)

        not_in_g2_swapper = list(set(g2.tree) - set(g2_swapper))
        for n in not_in_g2_swapper:
            g2_swapper[n.species].append(n)


        swap_spec, g1_nname = zip(*[(x.species, x.name) for x in g1swap])
        g2_nname = [x.name for x in g2swap]
        for spec in set(swap_spec):
            cand_replace1, cand_replace2 = set(g1_nname) - set(g2_nname), set(g2_nname) - set(g1_nname) 
            for n in g1_swapper[spec]:
                if n.name in g2_nname:
                    n.name = set.pop(cand_replace1)
            for n in g2_swapper[spec]:
                if n.name in g1_nname:
                    n.name = set.pop(cand_replace2)     
       
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
    def branch_is_valid(branch, leafnodes):
        # branch should be present 
        if len(leafnodes) < 2:
            return False
        else:
            multbrnch = set([x[-1] for x in leafnodes]) 
            # branch should be present on both parent trees
            if len(multbrnch) < 2:
                return False
            # branch should not be a single duplication node
            ln = Counter([y for y in x[0] for x in leafnodes ])
            if len(set(branch.part1)) == 1 and len(set(ln.values()))< 2 and len(branch.part1)==2:
                return False
            return True 
         

    @staticmethod
    def crossover(genome, **args):
        gdad = args['dad']
        gmom = args['mom']
        genome1 = gdad.clone()
        genome2 = gmom.clone()
        prob = max(genome1.intbrnp, genome2.intbrnp)
        current_dict = ddict(list)
        genomes = [genome1, genome2]
        for (gind, g) in enumerate(genomes):
            for part1, part2, nodes, imp in Utils.iter_partition(g.tree):
                brnch = BranchRep(part1, part2, imp)
                current_dict[brnch].append((nodes,gind))
        
        candidates = [(k,v) for k,v in current_dict.items() if Utils.branch_is_valid(k,v)]
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

    def __init__(self, tree, model, intbrnp=0.95, gmap={}, is_init=False):
        GenomeBase.__init__(self)
        self.tree =  tree
        self.model = model
        self.intbrnp = intbrnp
        self.is_init = is_init
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
        if not self.is_init:
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
        g.is_init = self.is_init
        
    def clone(self):
        """ Clone this GenomeBase
        :rtype: the clone genome
        .. note:: If you are planning to create a new chromosome representation, you
            **must** implement this method on your class.
        """
        newcopy = GPolySolver(None, None)
        self.copy(newcopy)
        return newcopy
