#include "variousFunctions.h"
#include "bulkScoring.h"
#include <cmath>

vector<string> split(string s){
	stringstream ss(s);
	string buf;
	vector<string> tokens;
	while (ss >> buf)
		tokens.push_back(buf);

	return tokens;
}

double getBulkVariance(Mutation* bulkMutations, int n){
	double VAFSum = 0.0;
	for(int i=0; i<n; i++)
		VAFSum += bulkMutations[i].getVAF();
	double averageVAF = VAFSum/n;

	double sumSquaredDistancesFromAverage = 0.0;
	for(int i=0; i<n; i++)
		sumSquaredDistancesFromAverage += (bulkMutations[i].getVAF()-averageVAF)*(bulkMutations[i].getVAF()-averageVAF);
	
	return sumSquaredDistancesFromAverage/n;
}


int getNumZerosInD(int**dataMatrix, int m, int n){
	int numZerosInD = 0;
	for(int i=0; i<m; i++)
		for(int j=0; j<n; j++)
			if(dataMatrix[i][j] == 0)
				numZerosInD += 1;

	return numZerosInD;
}

double getSCScoreScalingCoeff(int** dataMatrix, int m, int n, double alpha, double beta){
	int numZerosInD = getNumZerosInD(dataMatrix, m, n);
	return numZerosInD * (log(1-alpha)-log(beta));

}


double getBulkScoreScalingCoeff(Mutation* bulkMutations, int n){
	double averageVAF = 0.0;
	for(int i=0; i<n; i++)
		averageVAF += bulkMutations[i].getVAF();
	averageVAF = averageVAF/n;

	double result = 0;
	for(int i=0; i<n; i++){
		double differenceFromAverage = bulkMutations[i].getVAF() - averageVAF; 
		result += getVarianceCoeffInBulkScoring(bulkMutations[i]) * differenceFromAverage * differenceFromAverage;
	}
	
	return result;
}


/*
double getBulkScoreScalingCoeff(Mutation* bulkMutations, int n){
	double maxVAF = 0.0;
	for(int i=0; i<n; i++){
		if(bulkMutations[i].getVAF() > maxVAF){maxVAF = bulkMutations[i].getVAF();}
	}

	double result = 0;
	for(int i=0; i<n; i++){
		double differenceFromMaxVAF = bulkMutations[i].getVAF() - maxVAF; 
		result += getVarianceCoeffInBulkScoring(bulkMutations[i]) * differenceFromMaxVAF * differenceFromMaxVAF;
	}
	
	return result;
}
*/

string doubleToString(double x, int numDecimals){
	std::stringstream ss;
	ss << std::fixed << setprecision(numDecimals) << x;
	return ss.str();
}
