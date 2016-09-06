
# python libraries
import sys, os, re
import subprocess
import glob

import tempfile

# import RAxML SWIG module
import raxml

from ..TreeLib import TreeClass
from Bio import AlignIO



def calculate_likelihood(cmd, title, basedir=os.getcwd()):
    cmd = cmd+ "-n %s -w %s" % (title, os.path.abspath(basedir))
    rst = executeCMD(cmd)
    infofiles = glob.glob("%s/RAxML*.%s" % (basedir,title))
    f = [x for x in infofiles if 'info' in x][0]
    likelihoods = extractRAXMLikelihood(f, 1)[0]
    for f in infofiles:
        os.remove(f)
    return likelihoods
        

def executeCMD(cmd, dispout=False):
    
    p = subprocess.Popen(
        cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = p.communicate()
    #print "STDERR\n---------\n", err
    if dispout:
        print("\nSTDOUT\n---------\n", out)
    return err

def extractRAXMLikelihood(filename, n):
    likelihoods = []
    with open(filename) as IN:
        patern = re.compile('Tree [0-9]+: ')
        for line in reversed(IN.readlines()):
            if (patern.match(line)):
                likelihoods.append(float(line.split(':')[1].strip()))
                n -= 1
            if(n <= 0):
                break

    return list(reversed(likelihoods))


class LklModel():
    """computes statitic using raxml command line"""
    def __init__(self, alignment, cmd="raxmlHPC-SSE3", model="GTRGAMMA", eps=2.0, title="test", extra_string=""):
        self.alignment = alignment
        self.cmd = cmd
        fd, self.alignment = tempfile.mkstemp('.align')
        os.close(fd)
        AlignIO.write(alignment, self.alignment, format="phylip-relaxed")
        self.model = model
        self.eps = eps
        self.title  = title
        self.extra = extra_string 
        self.currLH = 0

    def optimize_model(self, gtree, **args):
        """Optimizes the RAxML model"""
        fd, treefile = tempfile.mkstemp('.tree')
        os.close(fd)
        gtree.write(outfile=treefile)
        cmdline = "%s -f g -z %s -s %s -m %s %s"%(self.cmd, treefile, self.alignment, self.model, self.extra)
        self.currLH = calculate_likelihood(cmdline, self.title, basedir=os.getcwd())
        os.remove(treefile)
        return self.currLH

    def print_raxml_tree(self, *args, **kargs):
        """Draw raxml tr -- adef and tr must have been previously defined"""
        print(self.curr_LH)
          


class RAxMLModel():
    """Computes test statistics using RAxML site-wise likelihoods"""

    def __init__(self, alignment, model="GTRGAMMA", eps=2.0, title="test", extra_string=""):
        """Initializes the RAxML model"""
        self._raxml = RAxML()
        self.model = model
        self.alignment = alignment
        fd, self.alignment = tempfile.mkstemp('.align')
        os.close(fd)
        AlignIO.write(alignment, self.alignment, format="phylip-relaxed")

        self.eps = eps
        self.title = title
        self.extra = extra_string

    def __del__(self):
        """Cleans up the RAxML model"""
        del self._raxml
        os.remove(self.alignment)


    def optimize_model(self, gtree):
        """Optimizes the RAxML model"""
        fd, treefile = tempfile.mkstemp('.tree')
        os.close(fd)
        gtree.write(outfile=treefile)
        self._raxml.optimize_model(treefile, self.alignment,
                                   "-m %s -e %s -n %s %s" % (self.model, self.eps, self.title, self.extra))
        os.remove(treefile)
        return self._raxml.best_LH

    def compute_lik_test(self, gtree, stat="SH", alternative=None):
        """Computes the test statistic 'stat' using RAxML likelihoods"""
        return self._raxml.compute_lik_test(gtree, stat, alternative)


    def print_raxml_tree(self, *args, **kargs):
        """Draw raxml tr -- adef and tr must have been previously defined"""
        treestr = raxml.tree_to_string(self._raxml.tr, self._raxml.adef)
        #tree = TreeClass(treestr)
        print(tree)
        #print(treestr)
        print(self._raxml.best_LH)
          

class RAxML:
    """Wrapper for RAxML functions"""

    #=========================================
    # constructors/destructors

    def __init__(self):
        self.adef = raxml.new_analdef()
        raxml.init_adef(self.adef)
        self.tr = raxml.new_tree()
        self.optimal = False
        self.best_LH = None; self.weight_sum = None; self.best_vector = None

    def __del__(self):
        raxml.delete_analdef(self.adef)
        raxml.delete_tree(self.tr)
        if self.best_vector is not None:
            raxml.delete_best_vector(self.best_vector)

    #=========================================
    # utilities

    def read_tree(self, tree):
        """Read treelib tree to raxml tr"""
        r,w = os.pipe()
        fr,fw = os.fdopen(r, 'r'), os.fdopen(w, 'w')

        tree.write(outfile=fw)
        fw.close()

        raxml.read_tree(fr, self.tr, self.adef)
        fr.close()

    #=========================================
    # model optimization

    def optimize_model(self, treefile, seqfile, extra="-m GTRGAMMA -n test"):
        """Optimizes the RAxML model"""

        # default model to use is GTRGAMMA
        # initialize parameters based on input
        cmd = "raxmlHPC -t %s -s %s %s" %\
              (treefile, seqfile, extra)
        raxml.init_program(self.adef, self.tr, cmd.split(' '))

        # optimize
        raxml.optimize_model(self.adef, self.tr)

        # reset best LH
        if self.best_vector is not None:
            raxml.delete_best_vector(self.best_vector)
        self.best_vector, self.best_LH, self.weight_sum = raxml.compute_best_LH(self.tr)

        # set flags
        self.optimal = True

    #=========================================
    # test statistics

    def compute_lik_test(self, tree, test="SH", alternative=None):
        """Computes the test statistic, returning the pvalue and Dlnl"""
        ##use scipy.stats to determine whether zscore is significant
        ##sf = 1 - cdf, zprob = cdf
        ##>>> stats.norm.sf(2)*2      # two-sided
        ##0.045500263896358417
        ##>>> (1-stats.zprob(2))*2    # two-sided
        ##0.045500263896358417
        ##>>> stats.zprob(2)
        ##0.97724986805182079
        ##>>> stats.norm.cdf(2)
        ##0.97724986805182079

        if test == "SH":
            if not self.optimal:
                raise Exception("The model is not optimized: call optimize_model.\n")

            self.read_tree(tree)
            zscore, Dlnl = raxml.compute_LH(self.adef, self.tr,
                                            self.best_LH, self.weight_sum, self.best_vector)

            # note that RAxML uses a one-sided comparison with a two-sided threshold
            # that is, it determines whether z>z_thr, where z_thr corresponds to a significance level of alpha/2
            # this is equivalent to testing sf(zscore)*2
            pval = sf(zscore)
        else:
            raise Exception("%s test statistic not implemented" % test)

        return pval, Dlnl
