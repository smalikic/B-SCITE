/*
*
* Created on Dec 11 2016 by Salem Malikic
*
* In this file definition  
*
*/
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <stdlib.h>
#include <iomanip>
#include <vector>
#ifndef MUTATION_H_
#define MUTATION_H_

using namespace std;

vector<string> splitString(string inputString, string delimitersString);

struct Mutation{
		string ID;
		string chromosome;
		int position;
		vector<int> mutReads;  // number of variant reads supporting mutation (size of array = number of samples)
		vector<int> refReads;  // number of reference reads spanning genomic position of mutation (size of array = number of samples)
		string INFO;
		// fields below are kept in order to spee-up model build time in cplex etc.
		int numSamples;        // number of bulk samples available = size of array mutReads (= size of array refReads)
		double* VAFs;
		double* bulkObjectiveQuadraticCoeffs;

		Mutation(){ }

		Mutation(string ID, string chromosome, int position, vector<int> mutReads, vector<int> refReads, string INFO){
			this->ID = ID;
			this->chromosome = chromosome;
			this->position = position;
			assert((mutReads.size() == refReads.size()) && "ERROR. Assertion in Mutation constructor failed. Unequal lengths of mutReads and refReads vectors");
			numSamples = mutReads.size();
			for(int i=0; i<numSamples; i++){
				this->mutReads.push_back(mutReads[i]);
				this->refReads.push_back(refReads[i]);
			}
			this->INFO = INFO;

			VAFs = new double[numSamples];
			for(int i=0; i<numSamples; i++){
				int totalReads = this->mutReads[i] + this->refReads[i];
				if(totalReads == 0){
					VAFs[i] = 0.0;
				}
				else{ 
					VAFs[i] = (2.0 * this->mutReads[i])/totalReads;}
				}
			
			bulkObjectiveQuadraticCoeffs = new double[numSamples];
			for(int i=0; i<numSamples; i++){
				double t = this->mutReads[i] + this->refReads[i];
				double z = 2 * (double(this->mutReads[i] + 0.5))/(t + 1.0);
				if( (0<z) && (z<2)){
					bulkObjectiveQuadraticCoeffs[i] = -t/(8*(z/2)*(1-z/2));
				}
				else if(z==0){
					bulkObjectiveQuadraticCoeffs[i] = 0;
				}
				else if(z==2){
					bulkObjectiveQuadraticCoeffs[i] = 0;
				}
				else{
					cout << "Assertion 0 <= z <= 2 failed in the constructor of Mutation class";
					assert(false);
				}
			}
		}
	
		
		// VAF used in the code is actually 2*vaf, where vaf is standard definition of variant allele frequency
		double getVAFinSample(int sampleIndex){ // 0-based indexing for samples
			// assert(sampleIndex < numSamples); -- commented to improve running time performance
			return VAFs[sampleIndex];
		}

		// although it is counterintuitive, we return string because in some cases (e.g. real data) we do not have information about true VAF
		string getTrueVAFinSample(int sampleIndex){
			vector<string> INFO_entries = splitString(INFO,";");
			for(int i=0; i<INFO_entries.size(); i++)
			{
				string currentEntry = INFO_entries[i];
				if(currentEntry.find_first_of("=")==string::npos){
					cout << "INFO field of bulk mutation must be in the form key1=value1;key2=value2;key3=value3;... Problem with " << currentEntry;
					assert(false);
				}
				string key   = splitString(currentEntry, "=")[0];
				string value = splitString(currentEntry, "=")[1];
				if(key == "trueVAF")
				{ 
					// ideally one needs to add here some code to verify that value represents float - not implemented as is only applicable to sim. data
					return value;
				}
			}

			return "NA";
		}
		
		int getNumSamples(){
			return numSamples;
		}
		
		string getID(){
			return ID;
		}
	
		double getQuadraticTermCoeffInBulkObjective(int sampleIndex){
			// assert(sampleIndex < numSamples); -- commented to improve running time performance
			return bulkObjectiveQuadraticCoeffs[sampleIndex];
		}

		string toString(){
			string stringRepresentation = "";
			stringRepresentation += ID + "\t";
			stringRepresentation += chromosome + "\t";
			stringRepresentation += to_string(position) + "\t";
			for(int sampleIndex=0; sampleIndex<numSamples-1; sampleIndex++){ stringRepresentation += to_string(mutReads[sampleIndex]) + ";";}
			stringRepresentation += to_string(mutReads[numSamples-1]) + "\t";
			for(int sampleIndex=0; sampleIndex<numSamples-1; sampleIndex++){ stringRepresentation += to_string(refReads[sampleIndex]) + ";";}
			stringRepresentation += to_string(refReads[numSamples-1]) + "\t";
			stringRepresentation += INFO;
			return stringRepresentation;
		}


		friend ostream& operator<<(ostream& os, Mutation& mut){
			os << mut.toString();
			return os;
		}
};

#endif
