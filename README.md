# GAPol
Genetic Algorithm for leaves labelling on a minimal reconciliation cost  gene tree.

GApol find the best tree (leaf labeling) from a list of unlabeled trees obtained using polytomysolver. The algorithm proceed by using a genetic algorithm to find the best candidates which are scored according to their likelihood. The binary is name 'gaperm'. 

gaperm --help


__optional arguments:__
  
    -h, --help      show this help message and exit
  
    -v, --version   show program's version number and exit
  
    --input TREES, -t TREES  file containing the trees in profileNJ output format
  
    --aln ALIGN, -s ALIGN   A sequence alignment file
  
    --alnfmt, -f {fasta,stockholm,clustal,nexus,maf,phylip}  The file format of the sequence alignment. 
                                                 The alignment is assumed to be in fasta format by default
  
    --smap SMAP, -S SMAP  Gene to species map. Use the standard format.
  
    --sep GENESEP     Gene-Specie separator for each leaf name in the genetree. 
                      This is an alternative for the --smap option
    
    --spos SPOS   The position of the specie name according to the separator. 
                  Supported options are 'prefix' and 'postfix'
    
    --cap         Capitalize the species name of the genetree leaves to  match each species.
                  Almost all functions are case sensitive.
    
    --corleaf     Change tree name to species
    
    --ignoredup   Check for duplicated trees and ignore them. This is useless when option --corleaf is used.
  
    --output OUTPUT, -o OUTPUT   Output file in which best trees should be saved
    
    --model, -m RAXMLMODEL,  Raxml model to use. If you do not provide this, 
                      it will guess your sequence type and use either GTRGAMMA or PROTGAMMALG
                      
    --eps RAXMLEPS        Raxml eps to use for likelihood computation
    
    --extras RAXMLEXTRA   Raxml extra arguments
    
    --timelim [TIMELIM]   Set time limit in minutes

__Genetic algo:__  Use genetic algorithm to find the

    --gen NGEN            Number of generations for the G.A.
  
    --popsize POPSIZE     Genetic population size
    
    --freqrep FREQREP     Number of generation before stat report
    
    --selector {roulette,tournament,rank,uniform}  Selector at each generation
    
    --besttree MLTREE     Path to known best tree, given the sequences data
    
    --verbose             Print GA status at each step
    
    --plot_lkl            Plot best ind likelihood for each generation

__All permutation:__ Test all possibility for all input trees

    --allsearch [ALLSEARCH]  Perform all permutation on all input and return best trees.
                          You can specifiy the number of trees you want in the output.
                          
                          
