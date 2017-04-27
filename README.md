# GAPol
Genetic Algorithm for leaves labelling on a minimal reconciliation cost  gene tree.

GApol find the best tree (leaf labeling) from a list of unlabeled trees obtained using polytomysolver. The algorithm proceed by using a genetic algorithm to find the best candidates which are scored according to their likelihood. The binary is name 'gaperm'. 

gaperm --help


__optional arguments:__
  
    -h, --help            show this help message and exit
    -v, --version         show program's version number and exit
    --input TREES, -t TREES
                          file containing the trees in profileNJ output format
    --aln ALIGN, -a ALIGN
                          A sequence alignment file
    --alnfmt {fasta,stockholm,clustal,nexus,maf,phylip}, -f {fasta,stockholm,clustal,nexus,maf,phylip}
                          The file format of the sequence alignment. The
                          alignment is assumed to be in fasta format by default
    --smap SMAP, -S SMAP  Gene to species map. Use the standard format.
    --sep GENESEP         Gene-Specie separator for each leaf name in the
                          genetree. This is an alternative for the --smap option
    --spos SPOS           The position of the specie name according to the
                          separator.
    --cap                 Capitalize the species name of the genetree leaves to
                          match each species. Almost all functions are case
                          sensitive.
    --ignoreduptop        Check for trees with same topology if leaves are not
                          labeled, then ignore them. If this is not used
                          overrepresenation will be consider done on purpose.
    --output OUTPUT, -o OUTPUT
                          Output file in which best trees should be saved
    --model RAXMLMODEL, -m RAXMLMODEL
                          Raxml model to use. If you do not provide this, it
                          will guess your sequence type and use either GTRGAMMA
                          or PROTGAMMALG
    --lklreestimate       Ensure that LKL model parameters will be re-estimate
                          each time - this will be longuer
    --extras RAXMLEXTRA   Raxml extra arguments
    --timelim [TIMELIM]   Set time limit in minutes
    --maxrcost MAXRCOST   Filtering input based on cost, if output is profileNJ
                          like
    --parallel [PARALLEL]
                          Set parallel mode for tree evaluation and
                          mutation/crossover
    --ignoreleaf          Ignore label in input tree (Transform profileNJ to
                          polytomysolver)
    --plot_lkl            Plot best ind likelihood for each generation



__Genetic algo:__  Use genetic algorithm to find the

    --gen NGEN            Number of generations for the G.A.
    --popsize POPSIZE     Genetic population size, will fill or remove initial
                          pop till this is met.
    --freqrep FREQREP     Frequency of report (each freqrep) generations
    --mutrate MUTRATE     Mutation rate
    --fastconv            Set the ga engine in accelerated mode, where each
                          generation divergence is increased. Can rapidly be
                          stuck in a local optimum
    --crossrate CROSSRATE
                          Crossover rate
    --elitism [ELITISM]   Number of elitist to bring to next generation
                          (according to their lkl only). If 0, elitism is
                          disabled
    --selector {roulette,tournament,rank,uniform}
                          Selector at each generation
    --besttree MLTREE     Path to known best tree, given the sequences data
    --verbose             Print GA status at each step
    --stopcrit {CONVD,FC,WC,SH}
                          Stopping criterion
    --alpha ALPHA         Threshold for WC and FC stopping criterion (0.95) and
                          alpha for SH criterion (0.05), Should be a float in
                          ]0,1[. Please avid FC and WC for small pop size
    --sloop SLOOP         Number of iteration for WC and FC stopping criterion
    --deltalkl DELTALKL   Maximum difference of score for CONV stopping
                          criterion


__All permutation:__ Test all possibility for all input trees

    --allsearch [ALLSEARCH]  Perform all permutation on all input and return best trees.
                          You can specifiy the number of trees you want in the output.
                          
                          
