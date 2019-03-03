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

const double pi= 3.1415926535897;
using namespace std; 

// pathOptimalBulkSolutionFile is file where values of optimal x (phi in the manuscript) and optimal y are stored. Additional INFO about mutation is also stored.
double bulkScoreTree(bool** anc_matrix,  Mutation* bulkMutations, int n, double w, string pathOptimalBulkSolutionFile){

	if((w==1.00) && (pathOptimalBulkSolutionFile == ""))
	{
		return 0;
	}

	int nodeOfMut[n]; // These variables are in fact the other representation of delta's
	for(int i=0; i<n; i++){nodeOfMut[i] = i;}
 
	try
	{
		// Cplex environment 
		IloEnv env;
		IloModel model(env);

		// Number of available bulk samples
		int numSamples = bulkMutations[0].getNumSamples();

		IloArray<IloFloatVarArray> x(env, n);
		for(int i=0; i<n; i++)
		{
			x[i] = IloFloatVarArray(env, numSamples, 0, 1);
		}

		for(int s=0; s<numSamples; s++)
		{
			IloExpr sum_x_Constraint(env);
			for(int i=0; i<n; i++)
			{
				sum_x_Constraint += x[i][s];
			}
			model.add(sum_x_Constraint <= 1);
		}


		IloExpr objective(env);	
		for(int i=0; i<n; i++)
		{
			int nodeOfCurrentMut = nodeOfMut[i]; // take the node to which the mutation is already assigned through single cell steps
			for(int s=0; s<numSamples; s++)
			{
				IloExpr errorOfCurrentMut(env);  // equals to (f_observed_in_bulk - f_observed_in_tree)^2
				errorOfCurrentMut +=  bulkMutations[i].getVAFinSample(s);
				for(int j=0; j<n; j++)
				{
					if(anc_matrix[nodeOfCurrentMut][j])
					{
						errorOfCurrentMut -= x[j][s];
					}
				}
				objective += bulkMutations[i].getQuadraticTermCoeffInBulkObjective(s) * errorOfCurrentMut * errorOfCurrentMut; 
			}		
		}

		model.add(IloMaximize(env, objective));


		IloCplex cplex(model);
		cplex.setParam(IloCplex::TiLim, 3 * 3600); // 3 hours time limit
		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());		
		cplex.setParam(IloCplex::EpGap, 0.002);
		cplex.setParam(IloCplex::Threads, 8);

		double startTime = cplex.getCplexTime();

		cplex.solve();

		double endTime = cplex.getCplexTime();
		double cplex_seconds = endTime - startTime;
		IloCplex::CplexStatus status = cplex.getCplexStatus();
		
		if(pathOptimalBulkSolutionFile != "")
		{
			double** optimal_x = allocate_doubleMatrix(n, numSamples);
			for(int i=0; i<n; i++)
			{	
				for(int s=0; s<numSamples; s++)
				{
					optimal_x[i][s] = cplex.getValue(x[i][s]);
				}
			}

			double** optimal_y = allocate_doubleMatrix(n, numSamples);	
			for(int i=0; i<n; i++)
			{
				for(int s=0; s<numSamples; s++)
				{
					optimal_y[i][s] = 0;
					for(int j=0; j<n; j++)
					{
						if(anc_matrix[i][j])
						{
							optimal_y[i][s] += optimal_x[j][s];
						}
					}
				}
			}

			ofstream optimalBulkSolutionFile;
			optimalBulkSolutionFile.open(pathOptimalBulkSolutionFile.c_str());
			optimalBulkSolutionFile << "MUT_ID\tTRUE_VAF\tINFERRED_VAF\t|TRUE-INFERRED|_VAF" << endl;
			for(int i=0; i<n; i++)
			{
				
				optimalBulkSolutionFile << bulkMutations[i].getID();
				optimalBulkSolutionFile << "\t";
				for(int s=0; s<numSamples; s++)
				{
					if (bulkMutations[i].getTrueVAFinSample(s) == "NA")
					{
						 optimalBulkSolutionFile << "NA" << ";";
					}
					else
					{ 
						optimalBulkSolutionFile << doubleToString(atof(bulkMutations[i].getTrueVAFinSample(s).c_str()), 5) << ";";
					}
				}
				optimalBulkSolutionFile << "\t";
				for(int s=0; s<numSamples; s++)
				{
					optimalBulkSolutionFile << doubleToString(optimal_y[i][s], 5) << ";";
				}
				optimalBulkSolutionFile << "\t";
				for(int s=0; s<numSamples; s++)
				{
					if(bulkMutations[i].getTrueVAFinSample(s) == "NA")
					{
						optimalBulkSolutionFile << "NA";
					}
					else
					{
						optimalBulkSolutionFile << doubleToString(abs(atof(bulkMutations[i].getTrueVAFinSample(s).c_str()) - optimal_y[i][s]), 5);
					}
					optimalBulkSolutionFile << ";";
				}
				optimalBulkSolutionFile << endl;
			}
			optimalBulkSolutionFile.close();
			free_doubleMatrix(optimal_x);
			free_doubleMatrix(optimal_y);
		}

		double objValue = cplex.getObjValue();
		env.end();

	return objValue;  

	} 

	catch (IloAlgorithm::CannotExtractException& e) { 
		cout << "Error in bulkScoreTree function. Problem with CPLEX.\n" << endl;
		cerr << "CannoExtractException: " << e << std::endl; 
		IloExtractableArray failed = e.getExtractables(); 
		for (IloInt i = 0; i < failed.getSize(); ++i) 
		{
			cerr << "\t" << failed[i] << endl; 
		}
	}
 
	catch(const IloException& e){
		cout << "ERROR! Cplex error in bulkScoreTree function." << endl << endl;
		cerr << e.getMessage();
	}
}


double bulkScoreTree(bool** anc_matrix,  Mutation* bulkMutations, int n, double w){
	return bulkScoreTree(anc_matrix, bulkMutations, n, w, "");
}
