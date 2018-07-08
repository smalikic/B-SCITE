#include<vector>
#include<string>
#include<sstream>
#include "Mutation.h"
#ifndef VARIOUSFUNCTIONS_H_
#define VARIOUSFUNCTIONS_H_

using namespace std;

vector<string> split(string s);  // analogous to python split() function 
double getBulkVariance(Mutation* bulkMutations, int n);
int getNumZerosInD(int** dataMatrix, int m, int n);
double getSCScoreScalingCoeff(int** dataMatrix, int m, int n, double alpha, double beta);
double getBulkScoreScalingCoeff(Mutation* bulkMutations, int n);
string doubleToString(double x, int numDecimals);
#endif
