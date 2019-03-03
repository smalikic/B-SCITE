#include<vector>
#include<string>
#include<sstream>
#include<fstream>
#include<algorithm>
#include<stdlib.h>
#include<assert.h>
#include "Mutation.h"

#ifndef VARIOUSFUNCTIONS_H_
#define VARIOUSFUNCTIONS_H_

using namespace std;

vector<string> split(string s);  // analogous to python split() function 
double getBulkVariance(Mutation* bulkMutations, int n);
int getNumZerosInBinaryMatrix(int** dataMatrix, int m, int n);
string doubleToString(double x, int numDecimals);
bool assertFileExists(string pathToFile);
#endif
