/*
*
* Created on Dec 11 2016 by Salem Malikic
*
* In this file definition  
*
*/

#include <string>

#ifndef MUTATION_H_
#define MUTATION_H_

using namespace std;

struct Mutation{
		string ID;
		string chromosome;
		int position;
		int mutReads, refReads; /** VAF is calculated as 2*var_reads/(var_reads + ref_reads) */
		string INFO;

		Mutation(){ }

		Mutation(string ID, string chromosome, int position, int mutReads, int refReads, string INFO){
			this->ID = ID;
			this->chromosome = chromosome;
			this->position = position;
			this->mutReads = mutReads;
			this->refReads = refReads;
			this->INFO = INFO;
		}
	
		double getVAF(){
			return (2.0 * mutReads)/(mutReads + refReads);
		}
		
		string getID(){
			return ID;
		}	
};
#endif
