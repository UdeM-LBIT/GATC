/* File: raxml.i */
%module raxml

%{
#define SWIG_FILE_WITH_INIT
#include "../../src/raxml/axml.h"
#include "../../src/raxml/globalVariables.h"
%}

%include typemaps.i
%include tmaps.i  // additional typemaps

%inline %{
/* struct helper functions */
analdef *new_analdef()
{
    return (analdef *)malloc(sizeof(analdef));
}

void delete_analdef(analdef *adef)
{
    free(adef);
}

tree *new_tree()
{
    return (tree *)malloc(sizeof(tree));
}

void delete_tree(tree *tr)
{
    free(tr);
}

/* wrapper for tree input and output */
void read_tree(FILE *fp, tree *tr, analdef *adef)
{
    treeReadLen(fp, tr, adef);
}

char *tree_to_string(tree *tr, analdef *adef)
{
    Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, FALSE, adef, SUMMARIZE_LH);
    return tr->tree_string;
}
%}

%inline %{
/* raxml axml.c: main */
void init_adef(analdef *adef)
{
    initAdef(adef);
}

/* raxml axml.c: main */
void init_program(analdef *adef, tree *tr,
                  int argc, char **argv)
{
    get_args(argc, argv, adef, tr);
}

/* raxml axml.c: main -> TREE_EVALUATION -> likelihood test */
void optimize_model(analdef *adef, tree *tr)
{
    rawdata *rdta = (rawdata *)malloc(sizeof(rawdata));
    cruncheddata *cdta = (cruncheddata *)malloc(sizeof(cruncheddata));

    if(adef->model == M_PROTCAT || adef->model == M_GTRCAT) {
        tr->rateHetModel = CAT;
    }
    else {
        tr->rateHetModel = GAMMA;
    }

    if(adef->useInvariant && adef->likelihoodEpsilon > 0.001)
        adef->likelihoodEpsilon = 0.001;

    readData(adef, rdta, cdta, tr);

    checkOutgroups(tr, adef);

    //makeFileNames();

    checkSequences(tr, rdta, adef);

    makeweights(adef, rdta, cdta, tr);
    makevalues(rdta, cdta, tr, adef);

    initModel(tr, rdta, cdta, adef);

    // case TREE_EVALUATION -> likelihood test
    getStartingTree(tr, adef);

    modOpt(tr, adef);
}
%}

%newobject compute_best_LH;
%delobject delete_best_vector;
%apply double *OUTPUT { double *bestLH, double *weightSum };
%inline %{
/* raxml axml.c: computeLH */
double *compute_best_LH(tree *tr, double *bestLH, double *weightSum)
{
    *bestLH = tr->likelihood;

    int i;
    *weightSum = 0.0;
    for(i = 0; i < tr->cdta->endsite; i++)
        *weightSum += (double)(tr->cdta->aliaswgt[i]);

    double *bestVector = (double*)malloc(sizeof(double) * tr->cdta->endsite);
    evaluateGenericInitrav(tr, tr->start);
    evaluateGenericVector(tr, tr->start, bestVector);
    /*
    double bestsum = 0;
    for(i = 0; i < tr->cdta->endsite; i++){
        printf("bestvector[%d] = %f\n", i, bestVector[i]);
        bestsum += bestVector[i];
    }
        printf("bestsum ==> %f\n",bestsum);
    */
    return bestVector;
}

void delete_best_vector(double *bestVector)
{
    free(bestVector);
}
%}
%clear double *bestLH, double *weightSum;

%apply double *OUTPUT { double *zscore, double *Dlnl };
%inline %{
/* raxml axml.c: computeLHTest */
void compute_LH(analdef *adef, tree *tr,
                double bestLH, double weightSum, double *bestVector,
                double *zscore, double *Dlnl)
{
    double currentLH;
    double *otherVector = (double*)malloc(sizeof(double) * tr->cdta->endsite);

    treeEvaluate(tr, 2);
    tr->start = tr->nodep[1];

    evaluateGenericInitrav(tr, tr->start);

    currentLH = tr->likelihood;
    if(currentLH > bestLH) {
        //printf("Better tree found at %f\n", currentLH);
        /*exit(1);*/
    }

    evaluateGenericVector(tr, tr->start, otherVector);

    {
        int j;
        double temp, wtemp, sum, sum2, sd;

        sum = 0.0;
        sum2 = 0.0;

        for (j = 0; j < tr->cdta->endsite; j++) {
            temp  = bestVector[j] - otherVector[j];
            wtemp = tr->cdta->aliaswgt[j] * temp;
            sum  += wtemp;
            sum2 += wtemp * temp;
        }

        sd = sqrt( weightSum * (sum2 - sum*sum / weightSum)
                   / (weightSum - 1) );

        //printf("Tree Likelihood: %f D(LH): %f SD: %f Significantly Worse: %s\n", currentLH, currentLH - bestLH, sd, (sum > 1.95996 * sd) ? "Yes" : " No");
        *zscore = sum/sd;
        *Dlnl = bestLH - currentLH;
    }

    free(otherVector);
}
%}
%clear double *zscore, double *Dlnl;

