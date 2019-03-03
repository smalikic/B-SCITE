/*
*
* bulkScoring.h  --- defines functions for scoring tree using bulk data and input provided by single-cell steps 
*
* Created by Salem Malikic on Dec 11 2016
*
*/
#include <string.h>
#include "Mutation.h"
#include "variousFunctions.h"
#include "matrices.h"
#include "output.h"

double bulkScoreTree(bool** anc_matrix, Mutation* bulkMutations, int n, double w, string pathOptimalBulkSolutionFile);
double bulkScoreTree(bool** anc_matrix, Mutation* bulkMutations, int n, double w);
double getConstPartOfObjective(Mutation* bulkMutations);
