/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *  
 *  and 
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 * 
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#ifndef GLOBALVARIABLES_H
#define GLOBALVARIABLES_H

#ifdef PARALLEL
extern int numOfWorkers;
#endif

extern int processID;
extern infoList iList;
extern FILE   *INFILE;

extern int Thorough;

extern char run_id[128], 
  workdir[1024], 
  seq_file[1024], 
  tree_file[1024], 
  weightFileName[1024], 
  modelFileName[1024], 
  excludeFileName[1024],
  bootStrapFile[1024], 
  permFileName[1024], 
  resultFileName[1024], 
  logFileName[1024], 
  checkpointFileName[1024], 
  infoFileName[1024], 
  randomFileName[1024],   
  bootstrapFileName[1024], 
  bipartitionsFileName[1024],
  ratesFileName[1024], 
  perSiteLLsFileName[1024], 
  lengthFileName[1024], 
  lengthFileNameModel[1024],
  proteinModelFileName[1024];

extern char *likelihood_key,
  *ntaxa_key,
  *smoothed_key;

extern char inverseMeaningDNA[16];
extern char inverseMeaningPROT[23]; 
extern double masterTime;

extern const int protTipParsimonyValue[23]; 

extern int partCount;

extern int optimizeRatesInvocations;  
extern int optimizeRateCategoryInvocations;
extern int optimizeAlphaInvocations;
extern int optimizeInvarInvocations;

#ifdef _USE_OMP
extern volatile int             NumberOfThreads;
#endif


#ifdef _USE_PTHREADS
extern volatile int             jobCycle;
extern volatile int             threadJob;
extern volatile int             NumberOfThreads;
extern volatile double          *reductionBuffer;
extern volatile double          *reductionBufferTwo;
extern volatile int             *reductionBufferParsimony;
extern volatile int             *barrierBuffer;
#endif

#endif // GLOBALVARIABLES_H
