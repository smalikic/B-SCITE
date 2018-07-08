/*
 * mcmc.h
 *
 *  Created on: Mar 7, 2017
 *      Author: Salem Malikic
 */
#include "variousFunctions.h"
#include "Mutation.h"
#include "matrices.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <utility>
#ifndef COMBINEDSCORESSTRUCT_H_
#define COMBINEDSCORESSTRUCT_H_

double calcCombinedSCBulkScore(double SCScore, double bulkScore, double w);
string roundDoubleToString(double argument, int numdecimals); // first rounds given double to numdecimals and then returns string representation

struct CombinedScoresStruct{
	double SCScore, bulkScore; // current values of SCITE and CITUP scores, as well as beta
	int n,m;
	bool** ancMatrix; // ancMatrix for which the optimal Combined Score is achieved
	double w; // Constant
	
	CombinedScoresStruct(int n, int m, double w){
		this->SCScore   = -100000.0;
		this->bulkScore = -100000.0;
		this->n = n;
		this->m = m;
		this->w = w;
		ancMatrix = allocate_boolMatrix(n, n);
	}

	double getScore(){
		return calcCombinedSCBulkScore(SCScore, bulkScore, w);
	}
	
	// If we observe that the new optimal solution has been found then the above function is used to print message about it to summaryFile
	// rep --- repetition,  it --- iteration
	void updateScore(int rep, int it, double propSCScore, double propBulkScore, bool** propAncMatrix, ofstream& summaryFile)
	{
		double propCombinedScore = calcCombinedSCBulkScore(propSCScore, propBulkScore, w);
		if(propCombinedScore <= getScore())
			return;
	
		string scoreType  = "COMBINED_SCORE";
		if(w == 1)
			scoreType = "SC_SCORE";
		if(w == 0)
			scoreType = "BULK_SCORE";	
		
		SCScore = propSCScore;
		bulkScore = propBulkScore;
		valuesCopy_boolMatrix(propAncMatrix, ancMatrix, n, n);
		// Rep It CominedScore ScaledSCScore ScaledBulkScore NonScaledSCore NonScaledBulkScore w
		summaryFile << setw(3)  << left << rep << setw(9) << left << it;
		summaryFile << setw(18) << left << scoreType;
		summaryFile << setw(10) << left << roundDoubleToString(-propCombinedScore, 2);
		summaryFile << setw(10) << left << roundDoubleToString(-w*SCScore, 2);
		summaryFile << setw(10) << left << roundDoubleToString(-(1-w)*bulkScore,  2);
		summaryFile << setw(10) << left << roundDoubleToString(-SCScore, 2);
		summaryFile << setw(10) << left << roundDoubleToString(-bulkScore, 2);
		summaryFile << setw(8)  << left << roundDoubleToString(w, 2);
		summaryFile << endl;
		
		return;
	}
	
	void printSummaryOfScores()
	{
		cout << "OVERALL OPTIMAL combined score  is:  " << getScore() << endl;
		cout << "	         where LOG score is:  " << SCScore << endl;
		cout << "		 and  BULK score is:  " << bulkScore << endl;
	}
};
		

#endif
