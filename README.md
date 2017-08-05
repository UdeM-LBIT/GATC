# GATC
Genetic Algorithm for gene Trees Construction/Correction.

GATC find the best tree according to both sequence data and species tree reconciliation using a genetic algorithm.

A typical run of the program starting with a population of trees to correct is :


```
gatc correct -t example/gtree.trees  --aln example/aln.fasta -S example/smap  -f 'fasta' --gen 50 --plot_lkl --output example/output  --timelim 120 --freqrep 1 --crossrate 0.8 --mutrate 0.6  --crit 'AU' --besttree example/mltree.nw  --parallel  --rectype par --dtlrate 2f 3f 1f -s example/sptree.nw --alpha 0.05 --nout 5
```

In this example, the dtl rates are fixed, by appending `"f"` to their values. GATC is used to correct the initial set of tree with the MOOP framework.


    usage: GATC [-h] [-v] [--plot_lkl] [--sample_space] [--nout NOUT]
                [--output OUTPUT] [--verbose] [--aln ALIGN]
                [--alnfmt {fasta,stockholm,clustal,nexus,maf,phylip}]
                [--seqmodel RAXMLMODEL] [--extras RAXMLEXTRA] [--gen NGEN]
                [--popsize POPSIZE] [--freqrep FREQREP] [--mutrate MUTRATE]
                [--crossrate CROSSRATE] [--elitism [ELITISM]]
                [--selector {roulette,tournament,rank,uniform}]
                [--parallel [PARALLEL]] [--smap SMAP] [--sep GENESEP]
                [--spos SPOS] [--sprconstr SPRCONSTR] [--use_weight aln recon]
                [--use_sigmoid] [--norec] [--rectype {par,lkl}] [--sptree SPTREE]
                [--dtlrate dup hgt loss] [--raterange low high sigma]
                [--edgerate mean sigma] [--discrsize DISCRSIZE]
                [--eventselector spr reroot dtl edge] [--stemlen STEMLEN]
                [--crit {CONV,FC,WC,SH,AU}] [--besttree MLTREE] [--alpha ALPHA]
                [--sloop SLOOP] [--deltalkl DELTALKL] [--timelim [TIMELIM]]
                [--allsearch [ALLSEARCH]]
                {correct,construct} ...

    GATC

    positional arguments_:
      {correct,construct}
        correct             Find best tree using a list of input tree
        construct           Construct tree from scratch

    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit
      --plot_lkl            Plot best ind likelihood for each generation
      --sample_space        Return all individuals in the population for each
                            generation. Will slow the algorithm
      --nout NOUT           Number of output trees to save
      --output OUTPUT, -o OUTPUT
                            Output file in which best trees should be saved
      --verbose             Print GA status at each step

    Sequence likelihood:
      Sequence likelihood parameters

      --aln ALIGN, -a ALIGN
                            A sequence alignment file
      --alnfmt {fasta,stockholm,clustal,nexus,maf,phylip}, -f {fasta,stockholm,clustal,nexus,maf,phylip}
                            The file format of the sequence alignment. The
                            alignment is assumed to be in fasta format by default
      --seqmodel RAXMLMODEL, -m RAXMLMODEL
                            Raxml model to use. If you do not provide this, it
                            will guess your sequence type and use either GTRGAMMA
                            or PROTGAMMALG
      --extras RAXMLEXTRA   Raxml extra arguments

    Genetic algo:
      Use genetic algorithm to find the best permutation

      --gen NGEN            Number of generations for the G.A.
      --popsize POPSIZE     Genetic population size. Genomes will be cloned or
                            removed from the initial population until the size
                            requirement is met
      --freqrep FREQREP     Frequency of report (each freqrep) generations
      --mutrate MUTRATE     Mutation rate
      --crossrate CROSSRATE
                            Crossover rate
      --elitism [ELITISM]   Number of elitist to bring to next generation
                            (according to their lkl only). If 0, elitism is
                            disabled
      --selector {roulette,tournament,rank,uniform}
                            Selector at each generation
      --parallel [PARALLEL]
                            Set parallel mode for tree evaluation and
                            mutation/crossover
      --smap SMAP, -S SMAP  Gene to species map. Use the standard format.
      --sep GENESEP         Gene-Specie separator for each leaf name in the
                            genetree. This is an alternative for the --smap option
      --spos SPOS           The position of the specie name according to the
                            separator.
      --sprconstr SPRCONSTR
                            Constraints for the SPR moves given the species tree,
                            *** Not implemented
      --use_weight aln recon
                            Use weight values for sequence likelihood and
                            reconciliation likelihood when computing fitness. The
                            MOOP framework will always be used for parsimony
                            reconciliation. Unit values will be used if not
                            provided.
      --use_sigmoid         Use sigmoid function with the provided weight. Sigmoid
                            will scale both function to 0-1.
      --norec               Only perform mutation/crossover that preserve
                            reconciliation cost, which mean that only sequence lkl
                            should be optimized
      --rectype {par,lkl}   Type of reconciliation to use: parcimonie or
                            likelihood. Parcimonie is faster

    Reconciliation parameters:
      Reconciliation parameters to set for algorithm

      --sptree SPTREE, -s SPTREE
                            HGT rate in a genome
      --dtlrate dup hgt loss
                            DTL cost to use, order: Dup, Trans, Loss
      --raterange low high sigma
                            Rates range limit when the rate is not fixed. New
                            values are obtained from a truncated normal
                            distribution with mean being the current value
      --edgerate mean sigma
                            Edge rate variation parameter under a iid gamma
      --discrsize DISCRSIZE
                            Interval discretisation size for the species tree
      --eventselector spr reroot dtl edge
                            Probability of selection for the following mutation:
                            SPR, REROOT, DTL, EDGE. Events probabilities are
                            redistributed when some hyperparameters are fixed.
      --stemlen STEMLEN     Stem length for the species tree

    Stopping criteria:
      Criteria to stop genetic algorithm

      --crit {CONV,FC,WC,SH,AU}
                            Stopping criterion
      --besttree MLTREE     Path to known best tree, given the sequences data
      --alpha ALPHA         Threshold for WC and FC stopping criterion (1- alpha) and
                            alpha for SH criterion (alpha = 0.05), Should be a float in
                            ]0,1[. Avoid FC and WC for small pop size
      --sloop SLOOP         Number of iteration for WC and FC stopping criterion
      --deltalkl DELTALKL   Maximum difference of score for CONV stopping
                            criterion
      --timelim [TIMELIM]   Set time limit in minutes

    All permutation:
      Test all possibility for all input trees

      --allsearch [ALLSEARCH]
                            Perform all permutation on all input and return best
                            trees. You can specifiy the number of trees you want.
                            If you set it to a negative value, all trees will be
                            returned
