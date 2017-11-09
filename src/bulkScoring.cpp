#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <ilconcert/iloenv.h>
#include <ilcplex/cplex.h>
#include <math.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <float.h>
#include <vector>
#include "bulkScoring.h"
#include "CombinedScoresStruct.h"

const double pi= 3.1415926535897;
const double min_x = 0.05;
using namespace std; 



double getVarianceCoeffInBulkScoring(Mutation mut){
	double halfVAF = mut.getVAF()/2;
	double totalReads = mut.mutReads + mut.refReads;
	if(halfVAF == 1)
		halfVAF = 0.99;
	if(halfVAF == 0)
		halfVAF = 0.01;
	if(halfVAF < 0 || halfVAF > 1){
		cout << "ERROR in bulk mutation frequency. Half VAF equals " << halfVAF << " but was set to 0.01." << endl;
		halfVAF = 0.01;
	} 
	double varianceCoefficient = totalReads/(8.0*halfVAF*(1-halfVAF));
	
	return varianceCoefficient;
}


double bulkScoreTree(bool** anc_matrix,  Mutation* mutations, int n, double* optimal_x, double* optimal_y, int* numCellsMutPresent, double w){
	if(w==1.00)
		return 0;

	bool quadraticError = true;
	
	int nodeOfMut[n]; // These variables are in fact the other representation of delta's
	for(int i=0; i<n; i++){
		nodeOfMut[i] = i; 
	};
 
	try
	{
		IloEnv env;
		IloModel model(env);
	
		IloFloatVarArray x = IloFloatVarArray(env, n, 0, 1);
		//IloSemiContVarArray x = IloSemiContVarArray(env, n, 0.05, 1);
		IloExpr sum_x_Constraint(env);
		for(int i=0; i<n; i++)
		{
			sum_x_Constraint += x[i];
		}
		model.add(sum_x_Constraint <= 1); /* Root is always kept for normal cells hence the purity is 1 - sum of all x */

		IloExpr objective(env);
		IloFloatVarArray absValues(env, n, 0, 1);	
		for(int i=0; i<n; i++)
		{
			if (numCellsMutPresent[i] == 0) continue;

			int nodeOfCurrentMut = nodeOfMut[i]; // take the node to which the mutation is already assigned through single cell steps
			
			IloExpr errorOfCurrentMut(env);  // equals to (f_observed_in_bulk - f_observed_in_tree)^2
			errorOfCurrentMut +=  mutations[i].getVAF();	
			for(int j=0; j<n; j++)
				if(anc_matrix[nodeOfCurrentMut][j])
					errorOfCurrentMut -= x[j];
	
	
			if(quadraticError)
			{
				objective += getVarianceCoeffInBulkScoring(mutations[i]) * errorOfCurrentMut * errorOfCurrentMut; 
			}
			else // !!!!! Analyze coefficient in this case
			{
				model.add(absValues[i] >= errorOfCurrentMut);
				model.add(absValues[i] >= -errorOfCurrentMut);
				objective += absValues[i];
			}		
		}

		IloObjective objectiveExpression = IloMinimize(env, objective);

		model.add(objectiveExpression);

		IloCplex cplex(model);
		cplex.setParam(IloCplex::TiLim, 3 * 3600); //3 hours time limit
		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());		
		cplex.setParam(IloCplex::EpGap, 0.002);
		cplex.setParam(IloCplex::Threads, 8);
		double startTime = cplex.getCplexTime();
		cplex.solve();
		double endTime = cplex.getCplexTime();
		double cplex_seconds = endTime - startTime;
		IloCplex::CplexStatus status = cplex.getCplexStatus();
		

		for(int i=0; i<n; i++)
		{
			optimal_x[i] = cplex.getValue(x[i]);
		}	
		
		for(int i=0; i<n; i++)
		{
			optimal_y[i]=0;
			for(int j=0; j<n; j++)
				if(anc_matrix[i][j])
					optimal_y[i] += optimal_x[j];
		}


		double objValue = cplex.getObjValue();
		env.end();

	/*
		 //The score below is the full score that accounts the approximation of binomial distribution using normal distribution
		double fullFinalScore = 0.0; // this score accounts for all summands (including constants) existing when using Normal approximation for Binomial 
		fullFinalScore += n*log(1.0/sqrt(2*pi));
		double half_VAF, variance_coefficient, totalCount;
		for(int i=0; i<n; i++)
		{
			half_VAF = mutations[i].getVAF()/2;
			totalCount = mutations[i].mutReads + mutations[i].refReads;
			double varianceCoeff = 1/(totalCount*(half_VAF)*(1-half_VAF));
			fullFinalScore += ((0.5) * log(varianceCoeff));
		}
		fullFinalScore -= objValue;
		return fullFinalScore;
		// if we want full Final score to be returned then return(fullFinalScore); statement should come here
	*/
		return -objValue;  

	} 

	catch (IloAlgorithm::CannotExtractException& e) { 
		cout << "ERROR_SALEM" << endl;
		std::cerr << "CannoExtractException: " << e << std::endl; 
		IloExtractableArray failed = e.getExtractables(); 
		for (IloInt i = 0; i < failed.getSize(); ++i) 
			std::cerr << "\t" << failed[i] << std::endl; 
 
	}
 
	catch(const IloException& e){
		cout << "ERROR_SALEM" << endl;
		cerr << e;
	}
}



Mutation* readBulkInput(string bulkFileLocation)
{
	ifstream bulkFile(bulkFileLocation, ifstream::in);
	if(!bulkFile.is_open()){
		cout << "There probably does not exist file " << bulkFileLocation << "\t EXITING !!!" << endl;
		assert(false);
	};
	

	string line;	
	vector<string> bulkFileLines;
	while(getline(bulkFile, line)){
		bulkFileLines.push_back(line);
	};
	bulkFile.close();
	
	bulkFileLines.erase(bulkFileLines.begin()); // remove header line
	int numMutations = bulkFileLines.size();
	Mutation* mutations = new Mutation[numMutations];
	for(int mutIndex = 0; mutIndex < numMutations; mutIndex++){
		vector<string> mutFields = split(bulkFileLines[mutIndex]);
		string ID = mutFields[0];
		string chromosome = mutFields[1];
		int position = stoi(mutFields[2]);
		int mutReads = stoi(mutFields[3]);
		int refReads = stoi(mutFields[4]);
		string INFO = mutFields[5];
		mutations[mutIndex] = Mutation(ID, chromosome, position, mutReads, refReads, INFO);
	 	//cout << mutIndex << "\t" << ID << "\t" << chromosome << "\t" << position << "\t" << mutReads << "\t" << refReads << "\t" << INFO << endl;	
	}	
		
	return mutations;
}


double absFunctionLocal(double x){
	if(x>0)
		return x;
	else
		return -x;
}

void writeOptimalVAFToFile(int n, Mutation* bulkMutations, double* optimal_x, double* optimal_y, string outFilenamePrefix){
	int width = 12;
	double min_significant_x   = 0.03;
	double accurracy_tolerance = 0.02;
	ofstream VAFSummaryFile;
	VAFSummaryFile.open((outFilenamePrefix + ".VAF").c_str(), ios::out);
	VAFSummaryFile << setw(width) << left << "mutID";
	VAFSummaryFile << setw(width) << left << "optimal_x";
	VAFSummaryFile << setw(width) << left <<  "optimal_y";
	VAFSummaryFile << setw(width) << left <<  "trueVAF";
	VAFSummaryFile << setw(18) << left << "ACCURRACY";
	VAFSummaryFile << setw(15) << left << ("optimal_x>" + doubleToString(min_significant_x, 2));
	VAFSummaryFile << endl;
	
	for(int i=0; i<n; i++)
	{
		VAFSummaryFile << setw(width) << left << bulkMutations[i].getID();
		VAFSummaryFile << setw(width) << setprecision(2) << left <<  optimal_x[i];
		VAFSummaryFile << setw(width) << setprecision(2) << left << optimal_y[i];
		VAFSummaryFile << setw(width) << setprecision(2) << left << bulkMutations[i].getVAF();

		if(absFunctionLocal(bulkMutations[i].getVAF() - optimal_y[i]) > accurracy_tolerance)
			VAFSummaryFile << setw(18) << left << "INACCURRATE";
		else
			VAFSummaryFile << setw(18) << left << "ACC";

		if(optimal_x[i] > min_significant_x)
			VAFSummaryFile << setw(15) <<  setprecision(2) << left << optimal_x[i];

		VAFSummaryFile << endl;	
	}
	VAFSummaryFile.close();	
}



void writeOptimalMatricesToFile(int n, CombinedScoresStruct& optimalCombinedScore, CombinedScoresStruct& optimalSCScore, CombinedScoresStruct& optimalBulkScore, string outFilenamePrefix){
	int parentVectorSize = n;
	ofstream MATRICES_SummaryFile;
	MATRICES_SummaryFile.open((outFilenamePrefix + "." + "matrices").c_str(), ios::out);
	MATRICES_SummaryFile << "NUM_MUTATIONS\t" << n << endl;

	MATRICES_SummaryFile << "ANCESTRY_MATRIX_OPTIMAL_COMBINED_SCORE" << endl;
	MATRICES_SummaryFile << ancMatrixToString(optimalCombinedScore.ancMatrix, parentVectorSize);
	MATRICES_SummaryFile << "PARENT_VECTOR_OPTIMAL_COMBINED_SCORE\t";
	int* parVecOptCombScoreTree = ancMatrixToParVector(optimalCombinedScore.ancMatrix, parentVectorSize);
	MATRICES_SummaryFile << parVectorToString(parVecOptCombScoreTree, parentVectorSize);
	delete [] parVecOptCombScoreTree;

	MATRICES_SummaryFile << "ANCESTRY_MATRIX_OPTIMAL_SC_SCORE" << endl;
	MATRICES_SummaryFile << ancMatrixToString(optimalSCScore.ancMatrix, parentVectorSize);
	MATRICES_SummaryFile << "PARENT_VECTOR_OPTIMAL_SC_SCORE\t";
	int* parVecOptSCScoreTree = ancMatrixToParVector(optimalSCScore.ancMatrix, parentVectorSize);
	MATRICES_SummaryFile << parVectorToString(parVecOptSCScoreTree, parentVectorSize);
	delete [] parVecOptSCScoreTree;

	MATRICES_SummaryFile << "ANCESTRY_MATRIX_OPTIMAL_BULK_SCORE" << endl;
	MATRICES_SummaryFile << ancMatrixToString(optimalBulkScore.ancMatrix, parentVectorSize);
	MATRICES_SummaryFile << "PARENT_VECTOR_OPTIMAL_BULK_SCORE\t";
	int* parVecOptBulkScoreTree = ancMatrixToParVector(optimalBulkScore.ancMatrix, parentVectorSize);
	MATRICES_SummaryFile << parVectorToString(parVecOptBulkScoreTree, parentVectorSize);
	delete [] parVecOptBulkScoreTree;

	MATRICES_SummaryFile.close();
}
