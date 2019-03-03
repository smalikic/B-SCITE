#include "CombinedScoresStruct.h"
#include <cmath>
#include <utility>
#include <assert.h>

double calcCombinedSCBulkScore(double SCScore, double bulkScore, double w ){
       return (w*SCScore + (1-w)*bulkScore);
}
