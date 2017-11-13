#include "CombinedScoresStruct.h"
#include <cmath>
#include <utility>
#include <assert.h>


double getSCContributionInCombinedScore(double SCScore, double bestIndependentSCScore, double w){
	return w * (SCScore-bestIndependentSCScore);
}

double getBulkContributionInCombinedScore(double bulkScore, double SCScoreScalingCoeff, double bulkScoreScalingCoeff, double w){	
	assert(bulkScoreScalingCoeff != 0);
	return (1-w) * bulkScore; // * SCScoreScalingCoeff)/bulkScoreScalingCoeff;
}

double calcCombinedSCBulkScore(double SCScore, double bulkScore, double bestIndependentSCScore, double SCScoreScalingCoeff, double bulkScoreScalingCoeff, double w ){
       return getSCContributionInCombinedScore(SCScore, bestIndependentSCScore, w) + getBulkContributionInCombinedScore(bulkScore, SCScoreScalingCoeff, bulkScoreScalingCoeff, w); 
}


