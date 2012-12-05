#ifndef __QUAD_CC__
#define __QUAD_CC__

#include "mcmc.h"


using namespace std;

QuadHandler::QuadHandler(parameterMap *parmap, MCMC * mcmc_in){
	mcmc = mcmc_in;
	// cout << "QuadHandler: Constructor Start" << endl;
	
	QuadScriptHome = parameter::getS(*parmap, "QuadFILEPATH", "./");
	//Observables = parameter::getS(mcmc->parmap, "OBSERVABLES", "Observables not specified");
	//Observables = mcmc->Observables;

	// cout << "The Quad is located at " << QuadScriptHome << endl;
	EmOutputFile = QuadScriptHome + "/src/QuadOutput.txt";
	
	fstream f;
	
	//create input/output files. I think this works, but it seems crude to me (fix later)
	f.open(EmOutputFile.c_str(), ios::out);
	f.flush();
	f.close();
	
	// cout << "QuadHandler: Constructor Done." << endl;
}

QuadHandler::~QuadHandler(){
	if(remove(EmOutputFile.c_str()) != 0){
		cerr << "Warning: ~QuadHandler: Unable to erase input/output/error files." << endl;
	}
	
	//cd back to original directory.
	string command = "cd " + string(pPath);
	int temp = system(command.c_str());
}

void QuadHandler::QueryQuad(ParameterSet Theta,vector<double> &Means, vector<double> &Errors){
	// cout << "Querying Quad." << endl;
	ofstream outputfile;
	ifstream inputfile;
	string command;
	string currentline;
	char * token;
	int NumDataRows = 1;
	EmulatedParams = mcmc->EmulatorParams;
	
	command = QuadScriptHome + "/src/quad ";
	for(int i = 0; i < Theta.Values.size(); i++){
		if(EmulatedParams.find(Theta.Names[i]) != string::npos){
			stringstream ss;
			ss << Theta.Values[i];
			command = command + ss.str() + " ";
		}
		else{
			cout << "Warning: Observable " << Theta.Names[i] << " is not an emulated observable." << endl;
		}
	}
	command = command + "> "+ EmOutputFile;
	
	int result = system(command.c_str());
	
	if(result != 0){
		cerr << "Error: Unable to run Quad." << endl;
		exit(1);
	}
	
	inputfile.open(EmOutputFile.c_str());
	if(inputfile){
		while(!inputfile.eof()){
			getline(inputfile, currentline, '\n');
			// cout << currentline << endl;
			if(currentline.compare(0,1,"#") != 0 && !currentline.empty() ){ //Comments have a # character
				if(currentline.compare(0,5,"Error") ==0){ //pass Quad errors to cerr stream
					cerr << "Error during emulation: Check " << EmOutputFile << endl;
					cerr << currentline << endl;
					exit(1);
				}
				
				stringstream ss;
				
				double tempnum;
				
				ss << currentline;
				while(ss >> tempnum){
					Means.push_back(tempnum);
					//Errors.push_back(1); //The errors are always 1
				}
				NumDataRows++;
			}
		}
		if(Means.size()==17){
			float myints[] = {23.2096,28.6014,45.2516,62.9723,0.390475,0.36656,0.383316,10.566,28.7311,45.3162,62.7765,0.00482302,0.00455732,0.00534056,0.326125,0.31816,0.326674};
			Errors.assign (myints,myints+17);
		}
		if(Means.size()==6){
			float myints[] = {1,1,1,1,1,1};
			Errors.assign (myints,myints+6);
		}
		
		inputfile.close();
	}
	else{
		cerr << "Unable to open Quad output file." << endl;
		exit(1);
	}
	if(Means.size() != Errors.size()){
		cerr << "Error: Quad output size mismatch. Error in reading Quad output in." << endl;
		cerr << "Mean vector size: " << Means.size() << " Error vector size: " << Errors.size() << endl;
		exit(1);
	}
	
}
#endif