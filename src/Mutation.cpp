#include "Mutation.h"

/*
	Function below serves as .split() function from Python. An example of output is:
	splitString("My name is Salem Malikic", " yk")
	output is vector consisting of strings: 
	["M", "name", "is", "Salem", "Mali", "ic"]
*/
vector<string> splitString(string inputString, string delimitersString){
	if(delimitersString.size() == 0){
		cout << endl << "ERROR in function splitString. Delimiters string can not be empty" << endl;}
	assert(delimitersString.size() > 0);	
	vector<string> delimiters;
	for(int i=0; i<delimitersString.size(); i++){
		delimiters.push_back(delimitersString.substr(i, 1));}

	vector<string> entries; // results of the split -- entry represents each of the elements of vector of strings returned as result
	string currentEntry = "";
	for(int i=0; i<inputString.size(); i++){
		string currentCharacter = inputString.substr(i,1);
		if(std::find(delimiters.begin(), delimiters.end(), currentCharacter) != delimiters.end()){
			if(currentEntry.size()>0){entries.push_back(currentEntry);}
			currentEntry = "";
		}
		else{
			currentEntry += currentCharacter;
		}
	}
	// collecting possibly the last element
	if(currentEntry.size()>0){
		entries.push_back(currentEntry);
		currentEntry = "";
	}
	return entries;
}


