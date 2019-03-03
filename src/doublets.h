/*
 * doublet.h
 *
 *  Created on: May 25, 2016
 *      Author: jahnka
 */

#ifndef DOUBLETS_H_
#define DOUBLETS_H_

std::string runMCMCbetaDoublet(std::vector<struct treeBeta>& bestTrees, double* errorRates, int noOfReps, int noOfLoops, double gamma, std::vector<double> moveProbs, int n, int m, int** dataMatrix, char scoreType, int* trueParentVec, int step, bool sample, double chi, double priorSd_beta, double priorSd_alpha, double doublet_alphaBetaMoveRatio, bool useTreeList, char treeType, double doubletProb, int doubleMut, Mutation* bulkMutations, double w, string outFilenamePrefix);

double getSampleSingletScoreFast(double* singletScores, int n, int*bft);
double sumScoreTreeFastSinglet(int n, int m, double** logScores, int** dataMatrix, int* parent);
double getCombinedProb(double L_0, double L_1, double p);
std::string sampleFromPosterior_doublets(double currTreeLogScore, int n, int* currTreeParentVec, double thetaProb, double currBeta, double currScore, double currAlpha, double currRelDoubletRate, double currDoubletRate);
std::vector<std::vector<double> > getDoubletAttachmentScoreMatrixFast(double* singletAttachmentScore, int n, int* bft, int* parent, double** logScores, int* dataVector);
std::vector<std::vector<double> > getDoubleMutDoubletAttachmentScoreMatrixFast(double* singletAttachmentScore, int n, int* bft, bool** anc, int* parent, double** logScores, int* obs, int doubleMut);
double getSampleDoubletScoreFast(std::vector<std::vector<double> > doubletAttachmentScore, int n, bool** ancMatrix, double& relDoubletScore);
bool isInvalidDoubleMutTree(int* parent, int doubleMut,int copy);
#endif /* DOUBLETS_H_ */
