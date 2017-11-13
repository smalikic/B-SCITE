/*
 * mcmcDoublet.cpp
 *
 *  Created on: May 18, 2016
 *      Author: jahnka
 */

#include <stdbool.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <time.h>
#include "matrices.h"
#include "treelist.h"
#include "trees.h"
#include "mcmc.h"
#include "scoreTree.h"
#include "scoreBinTree.h"
#include "rand.h"
#include "limits.h"
#include "output.h"
#include "mcmcBinTreeMove.h"
#include "mcmcTreeMove.h"
#include "doublets.h"
#include "bulkScoring.h"
#include "Mutation.h"
#include "variousFunctions.h"
//#include "doubleMut.h"
using namespace std;




#define XOR(first, second) ((first && !second) || (!first && second)) ? true : false
#define XNOR(first, second) !XOR(first, second)
#define AND(first, second) (first && second) ? true : (false)

class Beta_Distr {
	public:
    double mean, sd, alpha, beta;
    bool uniform;
  	public:
    Beta_Distr (double priorMean, double priorSD, bool uni);
    double logBetaPDF (double x) {return lgamma(alpha+beta)+(alpha-1)*log(x)+(beta-1)*log(1-x)-lgamma(alpha)-lgamma(beta);}
};


Beta_Distr::Beta_Distr (double priorMean, double priorSD, bool uni) {
	mean = priorMean;
	sd = priorSD;
	uniform = uni;
	if(uniform){
		alpha = 1.0;                                    // setting alpha and beta of the beta distribution for alpha to 1 gives a uniform prior
		beta = 1.0;
	}
	else{
		alpha = ((1-mean)*mean*mean/(sd*sd)) - mean;     // <-10.13585344 turn the mean and sd into parameters of the beta distribution
	  	beta = alpha*((1/mean)-1);                                                                 //<-13.38666556
	}
}

class MCMC_theta {
	public:
	double alpha;
	double beta;
	double alphaLogScore;
	double betaLogScore;
	double thetaLogScore;
	MCMC_theta (double newAlpha, double newBeta, double newAlphaLogScore, double newBetaLogScore, int m, int doublet_totalzeros);
	void update (double newAlpha, double newBeta, double newAlphaLogScore, double newBetaLogScore, int m, int doublet_totalzeros){

	}
};

MCMC_theta::MCMC_theta (double newAlpha, double newBeta, double newAlphaLogScore, double newBetaLogScore, int m, int doublet_totalzeros) {
	alpha = newAlpha;
	beta = newBeta;
	alphaLogScore = newAlphaLogScore;
	betaLogScore = newBetaLogScore;
	thetaLogScore = betaLogScore + alphaLogScore + m*doublet_totalzeros*log(1-alpha);
}

#define DOUBLET_UNIFORM_BETA_DISTR_beta      false
//#define UNIFORM_BETA_DISTR_beta      true    // uses a uniform beta distribution for the prior of beta, sets bpriora_beta=1 and bpriorb_beta=1
#define DOUBLET_UNIFORM_BETA_DISTR_alpha      false
//#define UNIFORM_BETA_DISTR_alpha      true  // uses a uniform beta distribution for the prior of alpha, sets bpriora_alpha=1 and bpriorb_alpha=1


double doublet_totalzeros = 100000000;           //  the number of sequenced sites observed to be mutated in no sample
double doublet_alphaBetaMoveRatio = 0.5;              // percentage of beta moves
bool   fixedDoubletProb = false;
double bestAvgSingletScoreSum = -DBL_MAX;
double bestAvgRelevantDoubletScoreSum = -DBL_MAX;
double propAvgSingletScoreSum = -DBL_MAX;
double propAvgRelevantDoubletScoreSum = -DBL_MAX;
double currAvgSingletScoreSum = -DBL_MAX;
double currAvgRelevantDoubletScoreSum = -DBL_MAX;

double precision = 1e-6;    // the precision of doubles basically
//default_random_engine generator;

double frac_w_Equals_1 = 0.0;


//COMMENT SALEM: can't we remove betterLogScore from the calculations below ?
/* combined probability for singlet and doublet for given p: (1-p)*L_0 + p * L_1 */
double getCombinedProb(double L_0, double L_1, double p){

	double betterLogScore = max(L_0, L_1);
	return log((1-p)*exp(L_0-betterLogScore) + p*exp(L_1-betterLogScore)) + betterLogScore;
}

/* gets the probability of sample being a doublet given p */
double getDoubletProb(double L_0, double L_1, double p){
	double betterLogScore = max(L_0, L_1);
	//if(p<=0.0){
//		p = 0.00005;                      // doublet prob is 0
//	}
//	else if(p >= 1.0){
//		p = 0.9999;                      // doublet prob is 1
//	}

		//cout << L_0 << "  " << L_1 << "  computing new score: " << p*exp(L_1) << " / (" << (1-p)*exp(L_0) << " + " << p*exp(L_1) << ")\n";
		//getchar();
		// COMMENT SALEM: Can't we remove betterScore ?
		 return p*exp(L_1-betterLogScore)/ ((1-p)*exp(L_0-betterLogScore) + p*exp(L_1-betterLogScore));
		 //return p*exp(L_1)/ ((1-p)*exp(L_0) + p*exp(L_1));

}


/* optimize doublet probability p given the tree and sample attachment probabilities as singlets or doublets (EM type optimization)*/
double optimizeDoubletProb(vector<double> L_0, vector<double> L_1, double p, int m){


//	if(p>=1){
//		p = 0.9999;
//	}
//	else if(p <= 0){
//		p = 0.0001;
//	}
	p = 0.5;
	double done = false;        // optimal p not yet found

	vector<double> sampleDoubletProb(m);
	int z=0;
	double currProbSum = 0.0;
	while(!done){
		double changed = false;
		double newProbSum = 0.0;


		z++; //cout << "round " << z++ << "\n";
		for(int sample=0; sample<m; sample++){
			double newProb = getDoubletProb(L_0[sample], L_1[sample], p);  // get the doublet probability for this sample
			newProbSum += newProb;


			sampleDoubletProb[sample] = newProb;


		}
		//cout.precision(17);
		//cout << "old/new p: " << p << "     " << newProbSum/m << "\n";
		if(abs(newProbSum-currProbSum) > precision ){
		//if(newProbSum != currProbSum){    // p has not yet converged

			changed = true;
			currProbSum = newProbSum;
		}
		if(!changed){       // p has converged
			break;
		}
		p = newProbSum/m;     // take the average probability over all samples
		//getchar();
	}
	//cout  << "new p: " << p << " after " << z  << " rounds\n";
	//getchar();
	return p;    // this is the best p for this tree when averaging over sample attachments
}


/* This runs the MCMC for learning the tree and beta, or only the tree with a fixed beta, it samples from the posterior and/or records the optimal trees/beta */
double scoreTreeFastWithDoublets(int*parent, int n, int m, double** logScores, int** dataMatrix, bool** ancMatrix, double p, double& best_p, double& relDoubletProb, int doubleMut){

	int treeSize = n+1;                   // for root node
	if(doubleMut>=0){ treeSize++; }       // in case of double mutation
	int parVecSize = treeSize-1;
	int* bft = getBreadthFirstTraversal(parent, parVecSize);   // get breadth first traversal for tree
	vector<double> singletSumScores(m);                          // the sum of likelihoods for single attachments for each sample
	vector<double> doubletSumScores(m);                          // the sum of likelihood for doublet attachments for each sample
	vector<double> relDoubletScores(m);                       // the sum of likelihoods of the relevant doublets
	vector<double> L_0(m);                        // the average likelihood for singlet attachment
	vector<double> L_1(m);                        // the average likelihood for doublet attachment

	double relDoubletScore = 0.0;
	for(int sample=0; sample<m; sample++){    // get singlet and doublet probabilities for all samples

		double* singletScores;                         // attachment scores of sample as singlet to each node
		vector<vector<double> > doubletScores;   // attachment scores of sample as doublet to each node pair

		
        	singletScores = getAttachmentScoresFast(parent, parVecSize, logScores, dataMatrix[sample], bft);
        	doubletScores = getDoubletAttachmentScoreMatrixFast(singletScores, parVecSize, bft, parent, logScores, dataMatrix[sample]);


		singletSumScores[sample] = getSampleSingletScoreFast(singletScores, parVecSize, bft);	// sum over all attachment points                                                                                                                  // there are no relevant doublets in a linear tree
		doubletSumScores[sample] = getSampleDoubletScoreFast(doubletScores, parVecSize, ancMatrix, relDoubletScore); // sum over all pairs of attachment points
		relDoubletScores[sample] = relDoubletScore;

		L_0[sample] = singletSumScores[sample] - log(treeSize);        // average over all attachment points
		L_1[sample] = doubletSumScores[sample] - 2*log(treeSize);      // average over all pairs of attachment points

		delete singletScores;
	}

	// optimize the doublet probability if the option was chosen
	if(fixedDoubletProb){
		best_p = p;
	}
	else{
		best_p = optimizeDoubletProb(L_0, L_1, p, m);    // computes the best p for this tree when averaging over sample attachments
	}


	double sumScore = 0.0;
	for(int sample=0; sample<m; sample++){
		relDoubletScore +=getDoubletProb(L_0[sample], L_1[sample], best_p) * relDoubletScores[sample];
		sumScore += getCombinedProb(L_0[sample], L_1[sample], best_p);
	}
	relDoubletProb = relDoubletScore/m;
	//cout << "sum score after:  " << sumScore << "\n";
	delete [] bft;
	//cout << "best p: " << best_p << "  best rel doublet rate: " << relDoubletScore << "\n";
	//getchar();
	return sumScore;
}





std::string runMCMCbetaDoublet(vector<struct treeBeta>& bestTrees, double* errorRates, int noOfReps, int noOfLoops, double gamma, vector<double> moveProbs, int n, int m, int** dataMatrix, char scoreType, int* trueParentVec, int step, bool sample, double chi, double priorSd_beta, double priorSd_alpha, bool useTreeList, char treeType, double doubletProb, int doubleMut, Mutation* bulkMutations, double w, string outFilenamePrefix, double* optimal_x, double* optimal_y, int* numCellsMutPresent, double bestIndependentSCScore){

	cout << endl << " ------ Starting MCMC with " << noOfReps << " repeats and " << noOfLoops << " loops ------" << endl;
	ofstream SCORES_SummaryFile;
//	SCORES_SummaryFile.open((outFilenamePrefix + ".scores").c_str(), ios::out);
//	SCORES_SummaryFile << "Rep\tIt\tScoreType\tScoreValue\tSCContribution\tbulkContribution\trawSCScore\trawBulkScore\tw" << endl;

	double burnInPhase = 0.25;                   // first quarter of steps are burn in phase
	unsigned int optStatesAfterBurnIn = 0;
	int burnIn = noOfLoops*burnInPhase;
	stringstream sampleOutput;

	//double changeDoubletProb = 0.0;     // probability that MCMC move changes beta
	double doubletBetaPriorMean = doubletProb;  // the prior probability for having a doublet
	//double doubletBetaPriorSd = 0.05;                                     //  prior sd for AD error rate
	//double doublet_jumpSd = doubletBetaPriorSd/chi;

	Beta_Distr alpha (errorRates[0], priorSd_alpha, DOUBLET_UNIFORM_BETA_DISTR_alpha);                  // prior distribution alpha (FD)
	Beta_Distr beta (errorRates[1] + errorRates[2], priorSd_beta, DOUBLET_UNIFORM_BETA_DISTR_beta);     // prior distribution beta (AD1 + AD2)


	int parentVectorSize = n;
	if(doubleMut >= 0){ parentVectorSize++; }

	if(treeType=='t'){parentVectorSize = (2*m)-2;}                     // transposed case: binary tree, m leafs and m-1 inner nodes, root has no parent


	double jumpSd_beta = beta.sd/chi;             // chi: scaling of the known error rate for the MH jump; resulting jump sd
	double jumpSd_alpha = alpha.sd/chi;           // chi: scaling of the known error rate for the MH jump; resulting jump sd


	int minDistToTrueTree = INT_MAX;             // smallest distance between an optimal tree and the true (if given)
	double bestTreeLogScore = -DBL_MAX;          // log score of T in best (T,beta)
	double bestScore = -DBL_MAX;                 // log score of best combination (T, beta)
	double bestBeta = beta.mean;
	double bestAlpha = alpha.mean;
	double bestRelDoubletRate = -1.0;
	double bestDoubletRate = doubletProb;



	double SCScoreScalingCoeff   = getSCScoreScalingCoeff(dataMatrix, m, n, bestAlpha, bestBeta);
	double bulkScoreScalingCoeff = getBulkScoreScalingCoeff(bulkMutations, n);
	CombinedScoresStruct optimalBulkScore     = CombinedScoresStruct(n, m, dataMatrix, bulkMutations, 0.0, bestAlpha, bestBeta, bestIndependentSCScore);
	CombinedScoresStruct optimalSCScore       = CombinedScoresStruct(n, m, dataMatrix, bulkMutations, 1.0, bestAlpha, bestBeta, bestIndependentSCScore);
	CombinedScoresStruct optimalCombinedScore = CombinedScoresStruct(n, m, dataMatrix, bulkMutations, w, bestAlpha, bestBeta, bestIndependentSCScore);

	clock_t lastTimeStamp = clock();
	for(int r=0; r<noOfReps; r++){   // repeat the MCMC, start over with random tree each time, only best score and list of best trees is kept between repetitions
		double topologySearch_w = 1.0; // w value used when deciding about move from current to next tree
		cout << endl << "REPEAT NUMBER: " << r << endl;

		int*   currTreeParentVec;
		if(treeType=='m'){
			if(doubleMut<0){
				currTreeParentVec = getRandParentVec(parentVectorSize);       // start MCMC with random tree
			}
			else{
				bool invalidDoubleMutTree = true;
				while(invalidDoubleMutTree==true){
					currTreeParentVec = getRandParentVec(parentVectorSize);
					invalidDoubleMutTree = isInvalidDoubleMutTree(currTreeParentVec, doubleMut, n);
				}
			}
		}
		else{             currTreeParentVec = getRandomBinaryTree(m);}                                                 // transposed case: random binary tree

		bool** currTreeAncMatrix =  parentVector2ancMatrix(currTreeParentVec,parentVectorSize);
		double** currLogScores = getLogScores(errorRates[0], errorRates[1], errorRates[2], errorRates[3]);           // compute logScores of conditional probabilities
		double currBeta = beta.mean;                                                                                  // the current AD rate
		double currAlpha = alpha.mean;
		double currDoubletProb = doubletBetaPriorMean;
		double currTreeLogScore;
		double propDoubletProb = currDoubletProb;    // this is optimized for the proposed tree in every MCMC loop
		double currRelDoubletProb = 0.0;
		double propRelDoubletProb = 0.0;
		double currTreeBulkScore;
		double currTreeCombinedScore;
		if(treeType=='m')
		{ 
			currTreeLogScore = scoreTreeFastWithDoublets(currTreeParentVec, n, m, currLogScores, dataMatrix, currTreeAncMatrix, currDoubletProb, propDoubletProb, propRelDoubletProb, doubleMut);
			currTreeBulkScore = bulkScoreTree(currTreeAncMatrix, bulkMutations, n, optimal_x, optimal_y, numCellsMutPresent, w);
			
			optimalSCScore.updateScore(r, -1, currTreeLogScore, currTreeBulkScore, currTreeAncMatrix, SCORES_SummaryFile);
			optimalBulkScore.updateScore(r, -1, currTreeLogScore, currTreeBulkScore, currTreeAncMatrix, SCORES_SummaryFile);
			optimalCombinedScore.updateScore(r, -1, currTreeLogScore, currTreeBulkScore, currTreeAncMatrix, SCORES_SummaryFile);
		        //sumScoreTreeFastSinglet(n, m, currLogScores, dataMatrix, currTreeParentVec);
		}
		else{              currTreeLogScore = getBinTreeScore(dataMatrix, n, m, currLogScores, currTreeParentVec);
						   // TODO: scoring of trees with doublets for the transposed case
		}

		double currBetaLogScore = (moveProbs[0]==0) ? 0.0 : beta.logBetaPDF(currBeta);            // zero if beta is fixed
		double currAlphaLogScore = (moveProbs[0]==0) ? 0.0 : alpha.logBetaPDF(currAlpha);          // zero if alpha is fixed
		double currThetaLogScore = currBetaLogScore + currAlphaLogScore + m*doublet_totalzeros*log(1-currAlpha);
		double currScore = currTreeLogScore + currThetaLogScore;
		//cout << "start score: " << currScore << "\n";
		
		currTreeCombinedScore = calcCombinedSCBulkScore(currTreeLogScore, currTreeBulkScore, bestIndependentSCScore, SCScoreScalingCoeff, bulkScoreScalingCoeff, topologySearch_w);
		int switch_w_It = int(frac_w_Equals_1 * noOfLoops); // defines index of iteration where w = 1 is replaced with actual w

		for(int it=0; it<noOfLoops; it++)
		{
        		
			if(it == switch_w_It)
			{
				topologySearch_w = w;
				currTreeCombinedScore = calcCombinedSCBulkScore(currTreeLogScore, currTreeBulkScore, bestIndependentSCScore, SCScoreScalingCoeff, bulkScoreScalingCoeff, topologySearch_w);
			}
			
			int reportWhenDivisibleBy = 1000;
			//int reportWhenDivisibleBy = noOfLoops/10;
			if(it % reportWhenDivisibleBy == 0){
				cout.precision(5);
				cout << "At MCMC (repeat,iteration) = (" << r << "," << it << ")\t and best score is " << optimalCombinedScore.getScore() << endl;
				clock_t currentTimeStamp = clock();
			//	//cout << "CLOCKS_PER SEC\t" << CLOCKS_PER_SEC << endl;
			//	cout << "Total time (in seconds) between latest reported (r,it) pair\t" << (double(currentTimeStamp-lastTimeStamp))/CLOCKS_PER_SEC << endl << endl; 
				lastTimeStamp = currentTimeStamp;
    			}

        		bool moveAccepted = false;                                           // Is the MCMC move accepted?
        		//bool moveChangesDoubletProb = false;

        		//moveChangesDoubletProb = changeBeta(changeDoubletProb);            // true if move changes the probability of having a doublet
        		//if(moveChangesDoubletProb==false){
        		bool moveChangesTheta = changeBeta(moveProbs[0]);                     // true if this move changes beta, not the tree


			/* Commented by SALEM, should be uncommented if we decide to learn alpha and beta
			if(moveChangesTheta)                                                 // new beta is proposed, log scores change tree is copy of current tree
			{
				double propBeta = proposeNewBeta(currBeta, jumpSd_beta);
				double propAlpha = proposeNewAlpha(currAlpha, jumpSd_alpha);

				if (sample_0_1() < doublet_alphaBetaMoveRatio){
					propBeta = currBeta;
				}
				else {
				    propAlpha = currAlpha;
				}

				//cout.precision(17);
				//cout << "new theta proposed: (" << propBeta << "," << propAlpha << ")\n";

				double** propLogScores = deepCopy_doubleMatrix(currLogScores, 4, 2);
				updateLogScoresAlphaBeta(propLogScores, propBeta, propAlpha);
				double propBetaLogScore = beta.logBetaPDF(propBeta);
					double propAlphaLogScore = alpha.logBetaPDF(propAlpha);
				double propThetaLogScore =  m*doublet_totalzeros*log(1-propAlpha) + propBetaLogScore + propAlphaLogScore;
				double propTreeLogScore;

				if(treeType=='m'){ propTreeLogScore = scoreTreeFastWithDoublets(currTreeParentVec, n, m, propLogScores, dataMatrix, currTreeAncMatrix, currDoubletProb, propDoubletProb, propRelDoubletProb, doubleMut);   // compute the new tree score for new beta
									//cout << "proposed tree score diff:  " << propTreeLogScore- currTreeLogScore << "     prop score: " << propTreeLogScore << "  curr score: " << currTreeLogScore << "\n";
									//cout << "proposed theta score diff: " << propThetaLogScore- currThetaLogScore << "      prop score: " << propThetaLogScore << "  curr score: " << currThetaLogScore << "\n";
				}
				else{              propTreeLogScore = getBinTreeScore(dataMatrix, n, m, propLogScores, currTreeParentVec);}

				if (sample_0_1() < exp((propTreeLogScore+propThetaLogScore-currTreeLogScore-currThetaLogScore)*gamma)){               // the proposed move is accepted
					moveAccepted = true;
					free_doubleMatrix(currLogScores);
				    currTreeLogScore  = propTreeLogScore;                                       // update score of current tree
				    currBeta = propBeta;                                                        // the current AD rate
				    currAlpha = propAlpha;
				    currBetaLogScore = propBetaLogScore;
				    currAlphaLogScore = propAlphaLogScore;
				    currThetaLogScore = propThetaLogScore;
				    currScore = currTreeLogScore+currThetaLogScore;                          // combined score of current tree and current beta
				    currLogScores = propLogScores;
				    currAvgSingletScoreSum = propAvgSingletScoreSum;
				    currAvgRelevantDoubletScoreSum	= propAvgRelevantDoubletScoreSum;
				    currDoubletProb = propDoubletProb;
				}
				else{
					delete [] propLogScores[0];
					delete [] propLogScores;
				}
			}
			else*/
			{                                   // move changed tree
				double nbhcorrection = 1.0;
        			int* propTreeParVec;
        			double propTreeLogScore;
				double propTreeBulkScore;
				double propTreeCombinedScore;
				bool** propTreeAncMatrix;
		
        			if(treeType=='m')
				{

        				propTreeParVec = proposeNewTree(moveProbs, parentVectorSize, currTreeAncMatrix, currTreeParentVec, nbhcorrection);  // propose new tree and
			
        				if(doubleMut>=0)
					{
        					while(isInvalidDoubleMutTree(propTreeParVec, doubleMut, n))
						{
        						delete [] propTreeParVec;
        						propTreeParVec = proposeNewTree(moveProbs, parentVectorSize, currTreeAncMatrix, currTreeParentVec, nbhcorrection);
        					}
        				}
		
        				propTreeAncMatrix = parentVector2ancMatrix(propTreeParVec, parentVectorSize);
        				propTreeLogScore = scoreTreeFastWithDoublets(propTreeParVec, n, m, currLogScores, dataMatrix, propTreeAncMatrix, currDoubletProb, propDoubletProb, propRelDoubletProb, doubleMut);
					propTreeBulkScore = bulkScoreTree(propTreeAncMatrix, bulkMutations, n, optimal_x, optimal_y, numCellsMutPresent, w);
					
					optimalSCScore.updateScore(r, it, propTreeLogScore, propTreeBulkScore, propTreeAncMatrix, SCORES_SummaryFile);
					optimalBulkScore.updateScore(r, it, propTreeLogScore, propTreeBulkScore, propTreeAncMatrix, SCORES_SummaryFile);
					optimalCombinedScore.updateScore(r, it, propTreeLogScore, propTreeBulkScore, propTreeAncMatrix, SCORES_SummaryFile);
        				//cout << "proposed score: " << propTreeLogScore << " after tree move, before tree score: " << currTreeLogScore << "\n";
        				free_boolMatrix(propTreeAncMatrix);
        			}
        			else
				{              
					propTreeParVec = proposeNextBinTree(moveProbs, m, currTreeParentVec, currTreeAncMatrix);
        		                propTreeLogScore = getBinTreeScore(dataMatrix, n, m, currLogScores, propTreeParVec);
				}

				propTreeCombinedScore = calcCombinedSCBulkScore(propTreeLogScore, propTreeBulkScore, bestIndependentSCScore, SCScoreScalingCoeff, bulkScoreScalingCoeff, topologySearch_w);
        			//if (sample_0_1() < nbhcorrection*exp((propTreeLogScore-currTreeLogScore)*gamma)) // commented by SALEM
				// the proposed tree is accepted
				double randomNumber_0_1 = sample_0_1();
				if(((propTreeCombinedScore - currTreeCombinedScore)>20) || (sample_0_1() < nbhcorrection*exp((propTreeCombinedScore - currTreeCombinedScore)*gamma)))
				{
        				moveAccepted = true;
        				free_boolMatrix(currTreeAncMatrix);                                            // discard outdated tree
        				delete[] currTreeParentVec;
        				currTreeAncMatrix = parentVector2ancMatrix(propTreeParVec,parentVectorSize); // update matrix of current tree
        				currTreeParentVec = propTreeParVec;                                         // update parent vector of current tree
        				currTreeLogScore  = propTreeLogScore;                                       // update score of current tree
        				currTreeBulkScore = propTreeBulkScore;
					currTreeCombinedScore = propTreeCombinedScore; 
					currScore = currTreeLogScore + currThetaLogScore;
        				currAvgSingletScoreSum = propAvgSingletScoreSum;
        				currAvgRelevantDoubletScoreSum	= propAvgRelevantDoubletScoreSum;
        				currDoubletProb = propDoubletProb;
        				currRelDoubletProb = propRelDoubletProb;
        				//cout  << "curr p = " << currDoubletProb << " current relevant doublet prob: " << currRelDoubletProb << "\n";
        			}
        			else
				{
        				delete [] propTreeParVec;            // discard proposed tree
        			}
        		}

        		/* If the true tree is given update the smallest distance between a currently best tree and the true tree */
        		if(trueParentVec)
			{
        			minDistToTrueTree = updateMinDistToTrueTree(trueParentVec, currTreeParentVec, parentVectorSize, minDistToTrueTree, currScore, bestScore);
        		}

        		/* If the list of optimal trees is used, update it */
        		if(useTreeList)
			{
        			updateTreeList(bestTrees, currTreeParentVec, parentVectorSize, currScore, bestScore, currBeta);
        		}

        		/* Sample from the posterior if required and past the burn-in phase */
        		if(sample && it>=burnIn && it % step == 0)
			{
        			sampleOutput << sampleFromPosterior_doublets(currTreeLogScore, parentVectorSize, currTreeParentVec, moveProbs[0], currBeta, currScore, currAlpha, currRelDoubletProb, currDoubletProb);
        		}


        		/* Update best tree in case we have found a new best one */
			// Below is commented by SALEM because we do not search for alpha, beta, doublet rate now
			/*  
     		 	if(currScore > bestScore)
			{
        			optStatesAfterBurnIn = 0;                    // new opt state found, discard old count
        			bestTreeLogScore = currTreeLogScore;
        			bestScore = currScore;                 // log score of best combination (T, beta)
        			bestBeta = currBeta;
        			bestAlpha = currAlpha;
        			bestRelDoubletRate = currRelDoubletProb;
        			bestDoubletRate = currDoubletProb;
        			//cout << "new best score: " << bestTreeLogScore << "   " << currDoubletProb << " current relevant doublet prob: " << currRelDoubletProb << "   all doublets rate: " << currDoubletProb << " beta: " << currBeta << "  alpha: " << currAlpha << "\n";
        		}
			

        		// Update the number of MCMC steps we spent in an optimal state 
	    		if(currScore == bestScore && it>=burnIn)
			{
        			optStatesAfterBurnIn++;
        		}
			*/			
			if(moveAccepted)
			{
				optStatesAfterBurnIn = 0;
			}

			if(currTreeCombinedScore == optimalCombinedScore.getScore() && it >= burnIn)
			{
				optStatesAfterBurnIn++;
			}

       		 } // the end of iteration
     
	   	delete [] currTreeParentVec;
        	free_doubleMatrix(currLogScores);
       		free_boolMatrix(currTreeAncMatrix);
	}                                              // last repetition of MCMC done



	unsigned int noStepsAfterBurnin = noOfReps*(noOfLoops-burnIn);
	cout.precision(17);
	cout << endl << endl;
	
	cout << "Best score found:\t" << optimalCombinedScore.getScore() <<  "\n";
	cout << "#optimal steps after burn-in:\t" << optStatesAfterBurnIn << "\n";
	cout << "total #steps after burn-in:\t" << noStepsAfterBurnin << "\n";
	cout << "#optimal steps after burn-in:\t" << (1.0*optStatesAfterBurnIn)/noStepsAfterBurnin << "\n";
	if(moveProbs[0]!=0.0){
		cout << "best value for beta:\t" << bestBeta << "\n";
		cout << "best value for alpha:\t" << bestAlpha << "\n";
		cout << "best relevant doublet rate:\t" << bestRelDoubletRate << "\n";
		cout << "best doublet rate:\t" << bestDoubletRate << "\n";
		cout << "best log score for (T, theta):\t" << bestScore << "\n";
	}

	optimalCombinedScore.printSummaryOfScores();
	//writeOptimalVAFToFile(n, bulkMutations, optimal_x, optimal_y, outFilenamePrefix);
	writeOptimalMatricesToFile(n, optimalCombinedScore, optimalSCScore, optimalBulkScore, outFilenamePrefix);
//	SCORES_SummaryFile.close();
	

	return sampleOutput.str();
}


bool isInvalidDoubleMutTree(int* parent, int doubleMut,int copy){
	if(parent[doubleMut] == parent[copy]){
		return true;
	}
	if(parent[doubleMut]==copy){
		return true;
	}
	if(parent[copy]==doubleMut){
		return true;
	}
	return false;
}

vector<vector<double> > getDoubletAttachmentScoreMatrixFast(double* singletAttachmentScore, int n, int* bft, int* parent, double** logScores, int* dataVector){

	vector<vector<double> > doubletAttachmentScore(n+1,vector<double>(n+1));
	int root = n;
	doubletAttachmentScore[root][root] = singletAttachmentScore[root];
	for(int i=1; i<=n; i++){                                      // score all attachment points in combination with root attachment and attachment pairs to the same node
		int node = bft[i];
		doubletAttachmentScore[node][root] = singletAttachmentScore[node];       // second cell attaches to root -> doublet score = singlet root score
		doubletAttachmentScore[root][node] = singletAttachmentScore[node];       // first cell attaches to root
		doubletAttachmentScore[node][node] = singletAttachmentScore[node];     // second cell attaches to same node
	}

	for(int i=1; i<=n; i++){                       // try all attachment points for the first cell of the doublet
		int first = bft[i];
		for(int j=i+1; j<=n; j++){                                                 // try all attachment points for the second cell (wlogs with larger bft-index)
			int second = bft[j];
			doubletAttachmentScore[first][second] = doubletAttachmentScore[first][parent[second]]; // doublet score for attaching second cell at parent of other node is known
			doubletAttachmentScore[first][second] -= logScores[dataVector[second]][0];            //  expected mutation state of the doublet for the mutation at second node
			doubletAttachmentScore[first][second] += logScores[dataVector[second]][1];            //  changes from 0 to 1; it cannot have already been 1 due to the bft and j>i
			doubletAttachmentScore[second][first] = doubletAttachmentScore[first][second];
		}
	}
	return doubletAttachmentScore;
}

vector<vector<double> > getDoubleMutDoubletAttachmentScoreMatrixFast(double* singletAttachmentScore, int n, int* bft, bool** anc, int* parent, double** logScores, int* obs, int doubleMut){

	vector<vector<double> > doubletAttachmentScore(n+2,vector<double>(n+2));
	int root = n+1;
	int smallerBftCopy = doubleMut;
	int biggerBftCopy = n;

	for(int i=0; i<=n; i++){
		if(bft[i] == n){
			smallerBftCopy = n;                        // the instance of the double mutation that has the higher bft-index
			biggerBftCopy = doubleMut;                  // the instance of the double mutation that has the lower bft-index
			break;
		}
		if(bft[i] == doubleMut){
			break;
		}
	}
//	if(bft[n] < bft[doubleMut]){
//		smallerBftCopy = n;                        // the instance of the double mutation that has the higher bft-index
//		biggerBftCopy = doubleMut;                  // the instance of the double mutation that has the lower bft-index
//	}

	bool backmutation = false;
	if(anc[smallerBftCopy][biggerBftCopy]){ backmutation = true; }     // the two copies of the double mutation are in one linage

	doubletAttachmentScore[root][root] = singletAttachmentScore[root];

	for(int i=1; i<=n+1; i++){                                          // score all attachment points in combination with root attachment and attachment pairs to the same node
		int node = bft[i];
		doubletAttachmentScore[node][root]  = singletAttachmentScore[node];     // second cell attaches to root -> doublet score = singlet root score
		doubletAttachmentScore[root][node]  = singletAttachmentScore[node];     // first cell attaches to root
		doubletAttachmentScore[node][node]  = singletAttachmentScore[node];     // second cell attaches to same node
	}

	for(int i=1; i<=n+1; i++){                 // try all non-root attachment points for the first cell of the doublet
		int smallerBftAtt = bft[i];

		if(anc[smallerBftCopy][smallerBftAtt] != anc[biggerBftCopy][smallerBftAtt]){  // double mutation expected to be present due to placement of 1st copy of doublet

			for(int j=i+1; j<=n+1; j++){                                             // loop through all attachment points for the 2nd copy (wlog with larger bft-index)
				int biggerBftAtt = bft[j];                   // attachment point for 2nd copy

				if(biggerBftAtt == biggerBftCopy){  // attaching to the second copy of the double mutation, no change compared to attachment to parent, as mutation already expected due to first attachment
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] = doubletAttachmentScore[smallerBftAtt][parent[biggerBftAtt]];
					doubletAttachmentScore[biggerBftAtt][smallerBftAtt] = doubletAttachmentScore[smallerBftAtt][biggerBftAtt];
				}
				else{                          // attaching to a node other than the second copy of the double mutation -> default update
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] = doubletAttachmentScore[smallerBftAtt][parent[biggerBftAtt]];

					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] -= logScores[obs[biggerBftAtt]][0];             //  expected mutation state of the doublet for the mutation at second node
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] += logScores[obs[biggerBftAtt]][1];             // changes from 0 to 1; it cannot have already been 1 due to the bft and j>i
					doubletAttachmentScore[biggerBftAtt][smallerBftAtt] = doubletAttachmentScore[smallerBftAtt][biggerBftAtt];
				}
			}
		}
		else{                                             // placement of first copy of doublet not enough to determine whether doublet mutation is expected to be present
			for(int j=i+1; j<=n+1; j++){                                 // loop through all attachment points for the 2nd copy (wlog with larger bft-index)
				int biggerBftAtt = bft[j];                            // attachment point for 2nd copy of doublet

				if( biggerBftAtt == biggerBftCopy && anc[smallerBftCopy][biggerBftAtt]){                                       // backmutation
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] = doubletAttachmentScore[smallerBftAtt][parent[biggerBftAtt]];
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] -= logScores[obs[doubleMut]][1];                 //  expected mutation state of the doublet for the mutation at second node
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] += logScores[obs[doubleMut]][0];                 // changes from 1 to 0 because backmutation
					doubletAttachmentScore[biggerBftAtt][smallerBftAtt] = doubletAttachmentScore[smallerBftAtt][biggerBftAtt];
				}
				else if(biggerBftAtt == smallerBftCopy || biggerBftAtt == biggerBftCopy){           // attaching to a parallel double mutation
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] = doubletAttachmentScore[smallerBftAtt][parent[biggerBftAtt]];
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] -= logScores[obs[doubleMut]][0];             //  expected mutation state of the doublet for the mutation at second node
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] += logScores[obs[doubleMut]][1];             // changes from 0 to 1; it cannot have already been 1 due to the bft and j>i
					doubletAttachmentScore[biggerBftAtt][smallerBftAtt] = doubletAttachmentScore[smallerBftAtt][biggerBftAtt];
				}
				else{                // attaching to a node other than a copy of the double mutation -> default update
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] = doubletAttachmentScore[smallerBftAtt][parent[biggerBftAtt]];
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] -= logScores[obs[biggerBftAtt]][0];             //  expected mutation state of the doublet for the mutation at second node
					doubletAttachmentScore[smallerBftAtt][biggerBftAtt] += logScores[obs[biggerBftAtt]][1];             // changes from 0 to 1; it cannot have already been 1 due to the bft and j>i
					doubletAttachmentScore[biggerBftAtt][smallerBftAtt] = doubletAttachmentScore[smallerBftAtt][biggerBftAtt];
				}
			}

		}
	}
	return doubletAttachmentScore;
}

/* computes L_1, the average likelihood over all doublet attachments of a sample */
/* also computes L_2, the by updating relDoubletScore, */
double getSampleDoubletScoreFast(vector<vector<double> > doubletAttachmentScore, int n, bool** ancMatrix, double& relDoubletScore){

	vector<double> doubletScores;            // list of the scores of all attachment pairs
	vector<double> relevantDoubletScores;    // list of the scores only for the relevant pairs

	for(int i=0; i<=n; i++){                                             // try all attachment points for the first cell of the doublet
		doubletScores.push_back(doubletAttachmentScore[i][i]);
		for(int j=i+1; j<=n; j++){                                       // try all attachment points for the second cell
			doubletScores.push_back(doubletAttachmentScore[i][j]);        // add score to list
			doubletScores.push_back(doubletAttachmentScore[j][i]);
			if(i<n && j<n && ancMatrix[i][j]==0 && ancMatrix[j][i]==0){                        // if the doublet is relevant, add to list
				relevantDoubletScores.push_back(doubletAttachmentScore[i][j]);
				relevantDoubletScores.push_back(doubletAttachmentScore[j][i]);
			}
		}
	}

	double doubletScoreSum  = 0.0;
	double relevantDoubletScoreSum = 0.0;

	doubletScoreSum = getSumOfVectorElementsExpLog(doubletScores);                    // sum over all doublet scores
	relevantDoubletScoreSum = getSumOfVectorElementsExpLog(relevantDoubletScores);   // sum over all relevant doublet scores
	relDoubletScore = exp(relevantDoubletScoreSum-doubletScoreSum);                   // L_1,rel / L_1
	return doubletScoreSum;    // return average doublet score
}


/* computes the sum of likelihoods over all singlet attachments of a sample */
double getSampleSingletScoreFast(double* singletScores, int n, int*bft){
	double bestSingletScore = getMaxEntry(singletScores, n+1);                 // best of the scores (used to compute with score differences rather than scores)
	double sumSingletScore = 0.0;

	for(int i=0; i<=n; i++){                                               // sum over all attachment scores, exp is necessary as scores are actually log scores
		sumSingletScore += exp(singletScores[bft[i]]-bestSingletScore);    // subtraction of best score to calculate with score differences
		//cout << "for node " << bft[i] << ": " << singletScores[bft[i]] << "    " << exp(singletScores[bft[i]]) << "\n";
	}

	return log(sumSingletScore)+bestSingletScore;          // transform back to log scores and change from score differences to actual scores

}

/* returns true if the tree is single chain, and false elsewise */
bool isLinearTree(int* parent, int* bft, int n){
	for(int i=0; i<n; i++){
		if(parent[bft[i]]==parent[bft[i+1]]){
			return false;
		}
	}
	return true;
}

/* computes an approximate scoring for a tree summing the score over all attachment points per sample */
/* this is basically the old scoring, and just used for comparison with the new scoring */
double sumScoreTreeFastSinglet(int n, int m, double** logScores, int** dataMatrix, int* parent){

	double sumTreeScore = 0.0;
	int* bft = getBreadthFirstTraversal(parent, n);   // get breadth first traversal for tree

	for(int sample=0; sample<m; sample++){
		double* scores = getAttachmentScoresFast(parent, n, logScores, dataMatrix[sample], bft); // attachments scores of sample to each node
		double bestMaxTreeScore = getMaxEntry(scores, n+1);                                     // best of the scores (used to compute with score differences rather than scores)

		double sumScore = 0.0;
		for(int i=0; i<=n; i++){                                                 // sum over all attachment scores, exp is necessary as scores are actually log scores
			sumScore += exp(scores[bft[i]]-bestMaxTreeScore);                   // subtraction of best score to calculate with score differences (smaller values)
		}
		delete [] scores;
		sumTreeScore += log(sumScore)+bestMaxTreeScore;                     // transform back to log scores and change from score differences to actual scores

	}
	//cout << "                                                      " << sumTreeScore << "\n";
	//getchar();
	delete [] bft;
	return sumTreeScore;
}

/* prints out the current tree and beta to sample from the posterior distribution */
string sampleFromPosterior_doublets(double currTreeLogScore, int n, int* currTreeParentVec, double thetaProb, double currBeta, double currScore, double currAlpha, double currRelDoubletRate, double currDoubletRate){

	std::stringstream content;
	content << currTreeLogScore  << "\t";                 // logscore of current tree
	content << countBranches(currTreeParentVec, n);       // number of branches in current tree
	if(thetaProb>0.0){
		content << "\t" << currBeta;                      // current beta
		content << "\t" << currAlpha;                      // current alpha
		content << "\t" << currScore;                     // current combined logscore for tree and theta
		content << "\t" << currRelDoubletRate;                // current relevant doublet rate
		content << "\t" << currDoubletRate;                // current doublet rate
	}
	content << "\t";
	for(int i=0; i<n; i++){
		content << currTreeParentVec[i] << " ";
	}
	content << "\n";
	return content.str();
}
