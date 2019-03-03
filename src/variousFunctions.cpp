#include "variousFunctions.h"
#include "bulkScoring.h"
#include <cmath>


/*
vector<string> split(string s){
	stringstream ss(s);
	string buf;
	vector<string> tokens;
	while (ss >> buf)
		tokens.push_back(buf);

	return tokens;
}
*/



int getNumZerosInBinaryMatrix(int**dataMatrix, int m, int n){
	int numZerosInD = 0;
	for(int i=0; i<m; i++)
		for(int j=0; j<n; j++)
			if(dataMatrix[i][j] == 0)
				numZerosInD += 1;

	return numZerosInD;
}



string doubleToString(double x, int numDecimals){
	std::stringstream ss;
	ss << std::fixed << setprecision(numDecimals) << x;
	return ss.str();
}


bool assertFileExists(string filePath){
	ifstream fin(filePath);
	if(fin.good()){
		fin.close();
		return true;
	}
	else{
		cout << "ERROR. There does not exist file " << filePath;
		assert(false);
	}	
	return false;
}
