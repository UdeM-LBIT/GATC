
# python libraries
import sys, os, re
import subprocess
import glob

import tempfile

# import RAxML SWIG module
import raxml

from ..TreeLib import TreeClass
from Bio import AlignIO


def get_rid_of(listfile):
    for f in listfile:
        os.remove(f)
    

def calculate_likelihood(cmd, title, ext="", basedir=os.path.abspath(os.getcwd()), size=1):
    cmd = cmd+ "-n %s -w %s" % (title+ext, basedir)
    rst = executeCMD(cmd)
    infofiles = glob.glob("%s/RAxML*.%s" % (basedir,title+ext))
    f = [x for x in infofiles if 'info' in x][0]
    likelihoods = extractRAXMLikelihood(f, size)
    get_rid_of(infofiles)
    return likelihoods

def compute_sh_test(cmd, title, basedir=os.path.abspath(os.getcwd())):
    cmd = cmd+ "-n %s -w %s" % (title, basedir)
    rst = executeCMD(cmd)
    infofiles = glob.glob("%s/RAxML*.%s" % (basedir,title))
    f = [x for x in infofiles if 'info' in x][0]
    results = extractSHTest(f)  
    get_rid_of(infofiles)
    return results

def run_consel(inputfile, type, sort=9, basedir=os.getcwd(), basename="RAxML_perSiteLLs"):
    makermtcmd = "makermt --%s %s" % (type, inputfile)
    fname = os.path.join(basedir, basename)
    conselcmd = "consel %s" % fname
    catpvcmd = "catpv %s > %s-pv.txt" % (fname, fname)
    if sort:
        catpvcmd += " -s %s" % sort
    executeCMD(makermtcmd)
    executeCMD(conselcmd)
    executeCMD(catpvcmd)
    conselOut = fname + "-pv.txt"
    return parseConselOutfile(conselOut, sort)

def parseConselOutfile(file, sort):
    title = []
    content = []
    first = True
    with open(file, 'r') as CONSELOUT:
        for line in CONSELOUT:
            if(line.startswith('#')):
                if first and 'reading' not in line and ("%s+" % sort not in line):
                    title.extend(
                        line.replace('#', '').replace('|', '').strip().split())
                    first = False
                elif not first:
                    values = line.replace('#', '').replace(
                        '|', '').strip().split()
                    content.append([float(x) for x in values])
    return dict(zip(title, zip(*content)))


def consel(cmd, title='consel', basedir=os.getcwd(), sort=9):
    cmd = cmd+ "-n %s -w %s" % (title, basedir)
    # run raxml
    executeCMD(cmd)
    infofiles = glob.glob("%s/RAxML*.%s" % (basedir,title))
    cons_input = glob.glob("%s/RAxML_perSiteLLs.%s" % (basedir,title))[0]
    # run consel
    consel_output = run_consel(cons_input, 'puzzle', sort=sort, basedir=basedir)
    infofiles.extend(glob.glob("%s/RAxML_perSiteLLs*" % (basedir)))
    #get_rid_of(infofiles)
    print consel_output
    sys.exit()
    return consel_output


def executeCMD(cmd, dispout=True):
    p = subprocess.Popen(
        cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = p.communicate()
    #print "STDERR\n---------\n", err
    if dispout:
        #print("\nSTDOUT\n---------\n", out)
        print("\nSTDERR\n---------\n")
        print(err)
        print(out)
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


def extractSHTest(filename):
    bestLK, treeLK = 0, 0
    alpha5, alpha2, alpha1 = False, False, False
    
    with open(filename) as IN:
        patern = re.compile('Model optimization, best Tree:')
        prevpatern = re.compile('Likelihood:')
        previous_line = ""
        for i, line in enumerate(reversed(IN.readlines())):
            if (patern.match(line)):
                bestLK = float(line.split(':')[-1].strip())
                if not prevpatern.search(previous_line):
                    raise IndexError("Cannot parse raxml info file")
                lline = previous_line.split(":")
                treeLK = float(lline[2].rstrip("D(LH)").strip())
                alpha5, alpha2, alpha1 = [not x.strip().upper().startswith("YES") for x in lline[5].split(",")]
                break
            previous_line = line
    return bestLK, treeLK, alpha5, alpha2, alpha1                        
    

class LklModel():
    """computes statitic using raxml command line"""
    def __init__(self, alignment, cmd="raxmlHPC-SSE3", model="GTRGAMMA", eps=2.0, title="test", extra_string="", reestimate=False):
        self.alignment = alignment
        self.cmd = cmd
        fd, self.alignment = tempfile.mkstemp('.align')
        os.close(fd)
        AlignIO.write(alignment, self.alignment, format="phylip-relaxed")
        self.model = model
        self.reestimate = reestimate
        self.eps = eps
        self.title  = title
        self.extra = extra_string 
        self.currLH = 0
        self.wdir = os.path.join(os.path.abspath(os.getcwd()), '_raxmltmp')
        if not os.path.exists(self.wdir):
            os.makedirs(self.wdir)

        get_rid_of(glob.glob("%s/RAxML*.%s" % (self.wdir,title)))

    def _build_cmd_line(self, treefile):
        bcmd = "-f g"
        if self.reestimate:
            bcmd = "-f G"
        cmdline = "%s %s -z %s -s %s -m %s %s"%(self.cmd, bcmd, treefile, self.alignment, self.model, self.extra)
        return cmdline

    def _build_treecomp_line(self, besttree, othertree):
        cmdline = "%s -f h -t %s -z %s -s %s -m %s %s"%(self.cmd, besttree, othertree, self.alignment, self.model, self.extra)
        return cmdline

    def optimize_model(self, gtree, **args):
        """Optimizes the RAxML model"""
        fd, treefile = tempfile.mkstemp('.tree')
        os.close(fd)
        size = 1
        if isinstance(gtree, list):
            size = len(gtree)
            with open(treefile, 'w') as GOUT:
                for gt in gtree:
                    GOUT.write(gt.write()+"\n")
        else:
            gtree.write(outfile=treefile)
        cmdline = self._build_cmd_line(treefile)
        self.currLH = calculate_likelihood(cmdline, self.title, ext=args.get("ext", ""), basedir=self.wdir, size=size)
        #print treefile
        os.remove(treefile)
        return self.currLH

    def print_raxml_tree(self, *args, **kargs):
        """Draw raxml tr -- adef and tr must have been previously defined"""
        print(self.curr_LH)
    
    def compute_consel_test(self, *args, **kwargs):
        """Compute consel output for a bunch of trees in argument"""
        # start by building a file with th two trees
        fd, treefile = tempfile.mkstemp('.trees')
        alpha = kwargs.get('alpha', 0.05)
        # get query tree position
        # best tree at second position by default
        querypos = kwargs.get('querypos', 2)
        os.close(fd)
        trees = []
        for t in args:
            trees.append(t.write())
        with open(treefile, 'w') as IN:
            IN.write("\n".join(trees))
        cmdline = self._build_cmd_line(treefile)
        consel_output = consel(cmdline, "consel", basedir=self.wdir)
        os.remove(treefile)
        item_pos = [int(x) for x in consel_output['item']].index(querypos)
        return consel_output['au'][item_pos] < alpha

    
    def compute_lik_test(self, besttree, tree, test="SH", alpha=0.05):
        """Compute sh test between two trees"""
        fd, besttreefile = tempfile.mkstemp('.btree')
        os.close(fd)
        fd, curtreefile = tempfile.mkstemp('.tree')
        os.close(fd)
        besttree.write(besttreefile)
        tree.write(curtreefile)
        cmdline  =  self._build_treecomp_line(besttreefile, curtreefile)
        if test == 'SH':
            bestlk, treelk, p5, p2, p1 = compute_sh_test(cmdline, self.title, basedir=self.wdir)
        else:
            raise NotImplementedError("%s test statistic not implemented" % test)
        os.remove(curtreefile)
        os.remove(besttreefile)
        p = {0.05:p5, 0.02:p2, 0.01:p1}
        return bestlk, treelk, p.get(alpha, False)    

    def __eq__(self, other):
        return (self.alignment == other.alignment) and (self.model == other.model) and (self.cmd == other.cmd)
        
      

class RAxMLModel():
    """Computes test statistics using RAxML site-wise likelihoods"""

    def __init__(self, alignment, model="GTRGAMMA", eps=2.0, title="test", extra_string=""):
        """Initializes the RAxML model"""
        self._raxml = RAxML()
        self.model = model
        self.alignment = alignment
        self.reestimate = True
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


    def optimize_model(self, gtree, **args):
        """Optimizes the RAxML model"""
        fd, treefile = tempfile.mkstemp('.tree')
        os.close(fd)
        gtree.write(outfile=treefile)
        self._raxml.optimize_model(treefile, self.alignment,
                                   "-m %s -e %s -n %s %s" % (self.model, self.eps, self.title+args.get("ext", ""), self.extra))
        os.remove(treefile)
        return self._raxml.best_LH

    def compute_lik_test(self, besttree, tree, test="SH", alpha=0.05):
        """Computes the test statistic 'stat' using RAxML likelihoods"""
        bestlk = self.optimize_model(besttree)
        pval, dnl = self._raxml.compute_lik_test(tree, test)
        return bestlk, None, pval>alpha 

    def compute_consel_test(self, *args, **kwargs):
        """Compute consel output for a bunch of trees in argument"""
        # start by building a file with th two trees
        fd, treefile = tempfile.mkstemp('.trees')
        alpha = kwargs.get('alpha', 0.05)
        os.close(fd)
        trees = []
        for t in args:
            print t
            trees.append(t.write())
        with open(treefile, 'w') as IN:
            IN.write("\n".join(trees))
        cmdline = self._build_cmd_line(treefile)
        consel_output = consel(cmdline, "consel", basedir=self.wdir)
        os.remove(treefile)
        return consel_output


    def print_raxml_tree(self, *args, **kargs):
        """Draw raxml tr -- adef and tr must have been previously defined"""
        #treestr = raxml.tree_to_string(self._raxml.tr, self._raxml.adef)
        #tree = TreeClass(treestr)
        #print(tree)
        #print(treestr)
        print(self._raxml.best_LH)
    

    def __eq__(self, other):
        return self.alignment == other.alignment and self.model == other.model


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
            raise NotImplementedError("%s test statistic not implemented" % test)

        return pval, Dlnl
