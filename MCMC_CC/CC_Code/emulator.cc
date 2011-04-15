#ifndef __EMULATOR_CC__
#define __EMULATOR_CC__

#include "mcmc.h"


using namespace std;

EmulatorHandler::EmulatorHandler(MCMC *mcmc_in){
	mcmc = mcmc_in;
	
	EmulatorScriptHome = parameter::getS(mcmc->parmap, "EMULATOR_FILEPATH", "/Users/kevinnovak/Research/RHIC_Research/madai-analysis");
	EmInputFile = EmulatorScriptHome + "/src/InputPts.txt";
	EmOutputFile = EmulatorScriptHome + "/src/EmulatorOutput.txt";
	
	fstream f;
	
	
	//create input/output files. I think this works, but it seems crude to me (fix later)
	f.open(EmInputFile.c_str(), ios::out);
	f << flush;
	f.close();
	
	f.open(EmOutputFile.c_str(), ios::out);
	f << flush;
	f.close();
	
	
	//check if emulator has been run yet.
	string checkfilename = EmulatorScriptHome + "/model-data/" + mcmc->dir_name + "/theta-tables.dat";
	f.open(checkfilename.c_str());
	if(f){
		f.close();
	}else{
		cerr << "EmulatorHander: Emulator doesn't exist for this project yet!" << endl;
		exit(1);
	}
	
	//get original PWD. Have to cd to the emulator directory, this is so we can cd back.
	pPath = getenv ("PATH");
	
	int result = system(("cd " + EmulatorScriptHome + "/src").c_str());
	
	if(result != 0){
		cerr << "Error: Unable to cd to emulator directory. The given directory is: " << EmulatorScriptHome << endl;
		exit(1);
	}
	
}

EmulatorHandler::~EmulatorHandler(){
	if(remove(EmInputFile.c_str()) != 0 || remove(EmOutputFile.c_str()) != 0){
		cerr << "Warning: ~EmulatorHandler: Unable to erase input/output files." << endl;
	}
	
	//cd back to original directory.
	string command = "cd " + string(pPath);
	int temp = system(command.c_str());
}

void EmulatorHandler::QueryEmulator(ParameterSet Theta,vector<double> &Means, vector<double> &Errors){
	ofstream outputfile;
	ifstream inputfile;
	string command;
	string currentline;
	char * token;
	int NumDataRows = 1;
	
	outputfile.open(EmInputFile.c_str());
	
	if(outputfile){
		for(int i = 0; i < Theta.Values.size(); i++){
			outputfile << Theta.Values[i];
			if(i != Theta.Values.size()-1){
				outputfile << " ";
			}
		}
		outputfile << endl;
		outputfile.close();
	}else{
		cerr << "Error: Unable to open input file to generate emulator input points." << endl;
		exit(1);
	}
	
	command = "cat " + EmInputFile + " | ./computePoints.sh > "+ EmOutputFile;
	int result = system(command.c_str());
	
	if(result != 0){
		cerr << "Error: Unable to run emulator." << endl;
		exit(1);
	}
	
	inputfile.open(EmOutputFile.c_str());
	if(inputfile){
		while(!inputfile.eof()){
			getline(inputfile, currentline, '\n');
			if(currentline.compare(0,1,"#") != 0 && !currentline.empty() ){ //Comments have a # character
				if(currentline.compare(0,5,"Error") ==0){ //pass emulator errors to cerr stream
					cerr << "Error during emulation: Check " << EmOutputFile << endl;
					cerr << currentline << endl;
					exit(1);
				}
				
				//not sure if this works.
				char * templine;
				strcpy(templine, currentline.c_str());
				token = strtok(templine, " ");
				//Since the emulator returns alternating rows of values and errors, if the row number is even its a row of error values
				if(NumDataRows % 2 == 0){
					while(token != NULL){
						Errors.push_back(atof(token));
						token = strtok(NULL," ");
					}
				}
				else{
					while(token != NULL){
						Means.push_back(atof(token));
						token = strtok(NULL," ");
					}
				}
			}
		}

		inputfile.close();
	}
	else{
		cerr << "Unable to open emulator output file." << endl;
		exit(1);
	}
	if(Means.size() != Errors.size()){
		cerr << "Error: Emulator output size mismatch. Error in reading emulator output in." << endl;
		exit(1);
	}
}
#endif