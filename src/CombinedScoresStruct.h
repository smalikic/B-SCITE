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

struct CombinedScoresStruct{
	double SCScore, bulkScore; // current values of SCITE and CITUP scores, as well as beta
	int n,m;
	bool** bestAncMatrix; // ancMatrix for which the optimal Combined Score is achieved
	bool** secondBestAncMatrix;
	double w; // Constant
	
	CombinedScoresStruct(int n, int m, double w){
		this->SCScore   = -100000.0;
		this->bulkScore = -100000.0;
		this->n = n;
		this->m = m;
		this->w = w;
		bestAncMatrix = allocate_boolMatrix(n, n);
		secondBestAncMatrix = allocate_boolMatrix(n,n);
	}

	double getScore(){
		return calcCombinedSCBulkScore(SCScore, bulkScore, w);
	}
	
	// If we observe that the new optimal solution has been found then the above function is used to print message about it to summaryFile
	// rep --- repetition,  it --- iteration
	void updateScore(int rep, int it, double propSCScore, double propBulkScore, bool** propAncMatrix, ofstream& summaryFile)
	{
		double propCombinedScore = calcCombinedSCBulkScore(propSCScore, propBulkScore, w);
		if(propCombinedScore < getScore())
			return;
	
		bool propTreeEqualToCurrentBestTree = true;
		for(int i=0; i<n;i++){
			if(propTreeEqualToCurrentBestTree == false){
				break;
			}
			for(int j=0; j<n; j++){
				if(propAncMatrix[i][j] != bestAncMatrix[i][j]){
					propTreeEqualToCurrentBestTree = false;
					break;
				}
			}
		}
		if(propTreeEqualToCurrentBestTree){
			return;
		}

		string scoreType  = "COMBINED_SCORE";
		if(w == 1)
			scoreType = "SC_SCORE";
		if(w == 0)
			scoreType = "BULK_SCORE";	
		
		SCScore = propSCScore;
		bulkScore = propBulkScore;
		valuesCopy_boolMatrix(bestAncMatrix, secondBestAncMatrix, n, n);
		valuesCopy_boolMatrix(propAncMatrix, bestAncMatrix, n, n);
		// Rep It CominedScore ScaledSCScore ScaledBulkScore NonScaledSCore NonScaledBulkScore w
		summaryFile << setw(3)  << left << rep << setw(9) << left << it;
		summaryFile << setw(18) << left << scoreType;
		summaryFile << setw(10) << left << doubleToString(-propCombinedScore, 2);
		summaryFile << setw(10) << left << doubleToString(-w*SCScore, 2);
		summaryFile << setw(10) << left << doubleToString(-(1-w)*bulkScore,  2);
		summaryFile << setw(10) << left << doubleToString(-SCScore, 2);
		summaryFile << setw(10) << left << doubleToString(-bulkScore, 2);
		summaryFile << setw(8)  << left << doubleToString(w, 2);
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
