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
#include "CombinedScoresStruct.h"

double getVarianceCoeffInBulkScoring(Mutation mut);
double bulkScoreTree(bool** anc_matrix, Mutation* muts, int n, double* optimal_x, double* optimal_y, int* numCellsMutPresent, double w);
Mutation* readBulkInput(string bulkFileLocation);
double absFunctionLocal(double x);
void writeOptimalVAFToFile(int n, Mutation* bulkMutations, double* optimal_x, double* optimal_y, string outFilenamePrefix);
void writeOptimalMatricesToFile(int n, CombinedScoresStruct& optimalCombinedScore, CombinedScoresStruct& optimalSCScore, CombinedScoresStruct& optimalBulkScore, string outFilenamePrefix);

