#ifndef __EMULATOR_CC__
#define __EMULATOR_CC__

#include "mcmc.h"


using namespace std;

EmulatorHandler::EmulatorHandler(parameterMap *parmap, MCMCConfiguration * mcmc_in){
	mcmc = mcmc_in;
	// cout << "EmulatorHandler: Constructor Start" << endl;
	
	EmulatorScriptHome = parameter::getS(*parmap, "EMULATORFILEPATH", "Couldn't find filepath");
	Observables = parameter::getS(mcmc->parmap, "OBSERVABLES", "Observables not specified");
	Cent_Range = parameter::getS(mcmc->parmap, "CENTRALITY_RANGE" , "Centrality range not specified");

	// cout << "The emulator is located at " << EmulatorScriptHome << endl;
	EmInputFile = EmulatorScriptHome + "/src/InputPts.txt";
	EmOutputFile = EmulatorScriptHome + "/src/EmulatorOutput.txt";
	EmErrorFile = EmulatorScriptHome + "/src/MCMCEmulatorError.txt";
	
	fstream f;
	
	//create input/output files. I think this works, but it seems crude to me (fix later)
	f.open(EmInputFile.c_str(), ios::out);
	f.flush();
	f.close();
	
	f.open(EmOutputFile.c_str(), ios::out);
	f.flush();
	f.close();
	
	f.open(EmErrorFile.c_str(), ios::out);
	f.flush();
	f.close();
	
	//check if emulator has been run yet.
	string checkfilename = mcmc->dir_name + "/" + Observables + "-" + Cent_Range + "-thetas.txt";
	// cout << checkfilename << endl;
	
	f.open(checkfilename.c_str());
	if(f){
		f.close();
		// cout << "Emulator exists." << endl;
	}else{ //If it can't find it:
		f.close();
		string checkfilename = "./" + Observables + "-" + Cent_Range + "-thetas.txt"; //Check the base directory
		f.open(checkfilename.c_str());
		if(f){
			f.close();
			// cout << "Emulator exists." << endl;
		}else{
			cerr << "EmulatorHandler: Emulator doesn't exist for this project yet!" << endl;
			exit(1);
		}
	}

	// cout << "EmulatorHandler: Constructor Done." << endl;
}

EmulatorHandler::~EmulatorHandler(){
	if(remove(EmInputFile.c_str()) != 0 || remove(EmOutputFile.c_str()) != 0 || remove(EmErrorFile.c_str()) != 0){
		cerr << "Warning: ~EmulatorHandler: Unable to erase input/output/error files." << endl;
	}
	
	//cd back to original directory.
	string command = "cd " + string(pPath);
	int temp = system(command.c_str());
}

void EmulatorHandler::QueryEmulator(ParameterSet Theta,vector<double> &Means, vector<double> &Errors){
	// cout << "Querying emulator." << endl;
	ofstream outputfile;
	ifstream inputfile;
	string command;
	string currentline;
	char * token;
	int NumDataRows = 1;
	EmulatedParams = mcmc->EmulatorParams;
	outputfile.open(EmInputFile.c_str());

	if(outputfile){
		if(EmulatedParams.find(Theta.Names[0]) != string::npos){
			outputfile << Theta.Values[0];
			// cout << Theta.Values[0];
			// cout << Theta.Names[0] << endl;
		}
		else{
			cout << "Warning: Parameter " << Theta.Names[0] << " is not an emulated parameter." << endl;
		}
		for(int i = 1; i < Theta.Values.size(); i++){
			if(EmulatedParams.find(Theta.Names[i]) != string::npos){
				outputfile << " " << Theta.Values[i];
				// cout << " " << Theta.Values[i];
				// cout << Theta.Names[i] << endl;
			}
			else{
				cout << "Warning: Parameter " << Theta.Names[i] << " is not an emulated parameter." << endl;
			}
		}
		outputfile << endl;
		// cout << endl;
		outputfile.close();
	}else{
		cerr << "Error: Unable to open input file to generate emulator input points." << endl;
		cout << "Filename: " << EmInputFile << endl;
		exit(1);
	}
	
	// command = "cat " + EmInputFile + " | " + EmulatorScriptHome + "/src/computePoints.sh  "\
	// + mcmc->dir_name + " > "+ EmOutputFile + " 2> " + EmErrorFile;
	
	command = EmulatorScriptHome + "/src/computePoints.sh " + mcmc->dir_name + " "\
	+ mcmc->dir_name + "/fn-data-" + Observables + "-" + Cent_Range + ".dat < " + EmInputFile + " > "+ EmOutputFile + " 2> " + EmErrorFile;

	// cout << command << endl;
	
	int result = system(command.c_str());
	
	if(result != 0){
		cerr << "Error: Unable to run emulator." << endl;
		cerr << "Check " << EmErrorFile << " for more details." << endl;
		exit(1);
	}
	
	inputfile.open(EmOutputFile.c_str());
	if(inputfile){
		while(!inputfile.eof()){
			getline(inputfile, currentline, '\n');
			// cout << currentline << endl;
			if(currentline.compare(0,1,"#") != 0 && !currentline.empty() ){ //Comments have a # character
				if(currentline.compare(0,5,"Error") ==0){ //pass emulator errors to cerr stream
					cerr << "Error during emulation: Check " << EmOutputFile << " and " << EmErrorFile << endl;
					cerr << currentline << endl;
					exit(1);
				}
				
				stringstream ss;
				
				double tempnum;
				
				ss << currentline;
				if(NumDataRows % 2 == 0){
					while(ss >> tempnum){
						Errors.push_back(tempnum);
					}
				}else{
					while(ss >> tempnum){
						Means.push_back(tempnum);
					}
				}
				NumDataRows++;
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