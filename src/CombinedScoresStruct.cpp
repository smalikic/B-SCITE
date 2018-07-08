#include "CombinedScoresStruct.h"
#include <cmath>
#include <utility>
#include <assert.h>

double calcCombinedSCBulkScore(double SCScore, double bulkScore, double w ){
       return 2*(w*SCScore + (1-w)*bulkScore);
}

string roundDoubleToString(double argument, int numDecimals){
	string strOfArgument = to_string(argument);
	int dotOccurence = strOfArgument.find(".");
	if(dotOccurence == string::npos)
		return strOfArgument;
	else
		return strOfArgument.substr(0, dotOccurence+3);
}	

 
