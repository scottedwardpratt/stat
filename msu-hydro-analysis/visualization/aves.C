// This code is for calcultaing the weighted averages and sigmas of the mcmc trace data.

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include <math.h>

using namespace std;

int main(int argc, char *argv[]){
	if(argc < 2){
		printf("Usage: ./visualisualization/aves model_directory trace_file_1 trace_file_2 ...\n");
		printf("Example: ./visualisualization/aves model_directory/Aprin2012 model_directory/April2012/trace/default\n");
		exit(-1);
	}

	int NUM_OF_PARS = 6;
	int NUM_OF_OBS = 19;
	int WRITEOUT = 150; // This number controlls how often computerPoints.sh is called.
	int x;
	ifstream trace,inputfile,tempfile;
	ofstream outputfile;
	string line,filename,command,EmInputFile,EmOutputFile,temp,tempfilename;
	int index;
	double param[WRITEOUT][NUM_OF_PARS],average[NUM_OF_PARS],variance[NUM_OF_PARS],totalweight;
	double likelihood,data[NUM_OF_OBS],dataerror[NUM_OF_OBS],observable[WRITEOUT][NUM_OF_OBS],observablerror[WRITEOUT][NUM_OF_OBS];
	char buffer[50];

	// Zero things
	for(int i=0;i<NUM_OF_PARS;i++){
		average[i]=0;
	}

	// Read in the real data
	ifstream datafile;
	string exp_data;
	string model_directory(argv[1]);
	exp_data = model_directory + "exp_data/results_trimmed.dat"; // At the moment there is no reading of parameter lists to make sure things line up, so make a trimmed list with only the used observables in order
	datafile.open(exp_data.c_str());
	getline(datafile,line,'\n');
	index = 0;
	while(!datafile.eof()){
		stringstream ss;
		ss << line;
		ss >> temp >> temp >> data[index] >> dataerror[index];
		//cout << data[index] << " " << dataerror[index] << endl;
		getline(datafile,line,'\n');
		index++;
	}

	index=0;

	EmInputFile = "input.txt";
	EmOutputFile = "output.txt";
	tempfilename = "filelisttemp.txt";
	fstream f;
	f.open(EmInputFile.c_str(), ios::out);
	f.flush();
	f.close();
	f.open(EmOutputFile.c_str(), ios::out);
	f.flush();
	f.close();
	f.open(tempfilename.c_str());
	f.flush();
	f.close();

	x=0;

	for(int i=2;i<argc;i++){
		string dirname(argv[i]);
		//cout << dirname << endl;
		string command = "ls " + dirname + " > " + tempfilename;
		int hold = system(command.c_str());
		tempfile.open(tempfilename.c_str());
		getline(tempfile,filename,'\n');

		string lastfile;
		lastfile = "trace.dat"; //This file doesn't have any new info in it

		while(!(filename==lastfile)){ // Read each file
			filename = dirname + "/" + filename;
			//cout << filename << endl;
			trace.open(filename.c_str());
			getline(trace,line,'\n'); // The first line is the names of the parameters and Theta0
			getline(trace,line,'\n');
			//cout << line;
			while(!trace.eof()){ // Read each line of each file
				stringstream ss;
				ss << line;
				ss.getline(buffer,6,',');
				index=atoi(buffer);
				for(int j=0;j<NUM_OF_PARS;j++){
					ss.getline(buffer,10,',');
					//parameter[i][j]=atoi(buffer);
					param[x][j]=atof(buffer);
					//cout << param[x][j] << " ";
				}
				//cout << x << " ";
				x++;
				if(x==WRITEOUT){ //After reading in WRITEOUT sets of params we will cacluate the likelihoods and averages
					//cout << "Getting Observables from Emulator" << endl;
					outputfile.open(EmInputFile.c_str()); // Now that you've read the params, calculate the observable
					outputfile << "#Not sure, does this line matter?" << endl;
					for(int i=0;i<WRITEOUT;i++){
						for(int j=0;j<NUM_OF_PARS;j++){
							outputfile << " " << param[i][j];
						}
						outputfile << endl;
					}
					outputfile.close();
					command = "./src/computePoints.sh " + model_directory + " " + model_directory + "/fn-data-standardlong.dat < " + EmInputFile + " > " + EmOutputFile;
					int result = system(command.c_str()); // Caculate the observables from the emulator

					command = "rm " +EmInputFile;
					result = system(command.c_str());

					// Get the calculated observables into memory so we can calculate the likelihood.
					inputfile.open(EmOutputFile.c_str());
					getline(inputfile,line,'\n');
					int linecount=0;
					while(!inputfile.eof()){
						if(line.compare(0,1,"#") != 0 && !line.empty() ){ // Skip comments and empty lines
							stringstream is;
							is << line;
							double obsbuff;
							//cout << linecount << " ";
							if((linecount % 2) == 0){ // Every other line is the errors
								for(int j=0;j<NUM_OF_OBS;j++){ // Read untill there is nothing more to read
									is >> obsbuff;
									observable[int(linecount/2)][j]=obsbuff;
									//cout << observable[int(linecount/2)][j] << " ";
									//cout << observable[int(linecount/2)][j] << " " << int(linecount/2) <<  " " << linecount << " " << j << endl;
								}
							}else{
								for(int j=0;j<NUM_OF_OBS;j++){ // Read untill there is nothing more to read
									is >> obsbuff;
									observablerror[int((linecount-1)/2)][j]=obsbuff;
								}
							}
							linecount++;
						}
						getline(inputfile,line,'\n');
					}
					inputfile.close();

					command = "rm " +EmOutputFile;
					result = system(command.c_str());

					for(int i=0;i<WRITEOUT;i++){ // Calculate the likelihood
						likelihood=1;
						for(int j=0;j<NUM_OF_OBS;j++){
							if(observablerror[i][j]<1) {observablerror[i][j]=1;}
							likelihood=likelihood*(1/sqrt(2*M_PI*pow(observablerror[i][j],2)))*exp(-pow((observable[i][j]-data[j]),2)/(2*pow(observablerror[i][j],2)));
							//cout << M_PI <<  " " << observable[i][j] << " " << data[j] << " " << observablerror[i][j] << endl;
						}
						//likelihood=likelihood/NUM_OF_OBS;
						//cout << log(likelihood) << endl;
						// And now that you have the likelihood, calculate the weighted average
						//likelihood=1;
						for(int j=0;j<NUM_OF_PARS;j++){
							average[j]+=param[i][j]*likelihood;
						}
						totalweight+=double(likelihood);
					}
					x=0;
					//cout << "================="<< endl;
				}
				getline(trace,line,'\n');
			}
			trace.close();
			filename.clear();
			getline(tempfile,filename,'\n');
			for(int i=0;i<x;i++){ // Calculate the likelihood
				likelihood=1;
				for(int j=0;j<NUM_OF_OBS;j++){
					if(observablerror[i][j]<1) {observablerror[i][j]=1;}
					likelihood=likelihood*(1/sqrt(2*M_PI*pow(observablerror[i][j],2)))*exp(-pow((observable[i][j]-data[j]),2)/(2*pow(observablerror[i][j],2)));
					//cout << M_PI <<  " " << observable[i][j] << " " << data[j] << " " << observablerror[i][j] << endl;
				}
				//likelihood=likelihood/NUM_OF_OBS;
				//cout << likelihood << endl;
				// And now that you have the likelihood, calculate the weighted average
				//likelihood=1;
				for(int j=0;j<NUM_OF_PARS;j++){
					average[j]+=param[i][j]*likelihood;
				}
				totalweight+=double(likelihood);
			}
		}
		tempfile.close();
	}

	cout << "Averages:" << endl;
	for(int j=0;j<NUM_OF_PARS;j++){
		average[j]=average[j]/totalweight;
		cout << average[j] << endl;
	}

	//Now that we have done everything once and have the average values for the parameters, we calculate the weighted variances
	for(int i=2;i<argc;i++){
		string dirname(argv[i]);
		//cout << dirname << endl;;
		string command = "ls " + dirname + " > " + tempfilename;
		int hold = system(command.c_str());
		tempfile.open(tempfilename.c_str());
		getline(tempfile,filename,'\n');

		string lastfile;
		lastfile = "trace.dat"; //This file doesn't have any new info in it

		while(!(filename==lastfile)){ // Read each file
			filename = dirname + "/" + filename;
			//cout << filename << endl;
			trace.open(filename.c_str());
			getline(trace,line,'\n'); // The first line is the names of the parameters and Theta0
			getline(trace,line,'\n');
			//cout << line;
			while(!trace.eof()){ // Read each line of each file
				stringstream ss;
				ss << line;
				ss.getline(buffer,6,',');
				index=atoi(buffer);
				for(int j=0;j<NUM_OF_PARS;j++){
					ss.getline(buffer,10,',');
					//parameter[i][j]=atoi(buffer);
					param[x][j]=atof(buffer);
				}
				x++;
				if(x==WRITEOUT){ //After reading in WRITEOUT sets of params we will cacluate the likelihoods and averages
					outputfile.open(EmInputFile.c_str()); // Now that you've read the params, calculate the observable
					outputfile << "#Not sure, does this line matter?" << endl;
					for(x=0;x<WRITEOUT;x++){
						for(int j=0;j<NUM_OF_PARS;j++){
							outputfile << " " << param[x][j];
						}
						outputfile << endl;
					}
					outputfile.close();
					command = "./src/computePoints.sh " + model_directory + " " + model_directory + "/fn-data-standardlong.dat < " + EmInputFile + " > " + EmOutputFile;
					int result = system(command.c_str()); // Caculate the observables from the emulator

					// Get the calculated observables into memory so we can calculate the likelihood.
					inputfile.open(EmOutputFile.c_str());
					getline(inputfile,line,'\n');
					int linecount=0;
					while(!inputfile.eof()){
						if(line.compare(0,1,"#") != 0 && !line.empty() ){ // Skip comments and empty lines
							stringstream is;
							is << line;
							double obsbuff;
							//cout << linecount << " ";
							if((linecount % 2) == 0){ // Every other line is the errors
								for(int j=0;j<NUM_OF_OBS;j++){ // Read untill there is nothing more to read
									is >> obsbuff;
									//observable[int(linecount/2)][j]=obsbuff;
									//cout << observable[int(linecount/2)][j] << " " << int(linecount/2) <<  " " << linecount << " " << j << endl;
								}
							}else{
								for(int j=0;j<NUM_OF_OBS;j++){ // Read untill there is nothing more to read
									is >> obsbuff;
									observablerror[int((linecount-1)/2)][j]=obsbuff;
								}
							}
							linecount++;
						}
						getline(inputfile,line,'\n');
					}
					inputfile.close();

					for(x=0;x<WRITEOUT;x++){ // Calculate the likelihood
						likelihood=1;
						for(int j=0;j<NUM_OF_OBS;j++){
							if(observablerror[i][j]<1) {observablerror[i][j]=1;}
							likelihood=likelihood*(1/sqrt(2*M_PI*pow(observablerror[x][j],2)))*exp(-pow((observable[x][j]-data[j]),2)/(2*pow(observablerror[x][j],2)));
							//cout << M_PI <<  " " << observable[x][j] << " " << data[j] << " " << observablerror[x][j] << endl;
						}
						//likelihood=likelihood/NUM_OF_OBS;
						//cout << likelihood << endl;
						// And now that you have the likelihood, calculate the weighted variance
						//likelihood=1;
						for(int j=0;j<NUM_OF_PARS;j++){
							double diff=(average[j]-param[x][j]);
							//cout << diff << " ";
							//if(diff<0){ diff=-diff; }
							variance[j]+=diff*diff*likelihood;
						}
					}
					x=0;
				}
				getline(trace,line,'\n');
			}
			trace.close();
			filename.clear();
			getline(tempfile,filename,'\n');
			for(int i=0;i<x;i++){ // Calculate the likelihood
				likelihood=1;
				for(int j=0;j<NUM_OF_OBS;j++){
					if(observablerror[i][j]<1) {observablerror[i][j]=1;}
					likelihood=likelihood*(1/sqrt(2*M_PI*pow(observablerror[i][j],2)))*exp(-pow((observable[i][j]-data[j]),2)/(2*pow(observablerror[i][j],2)));
					//cout << M_PI <<  " " << observable[i][j] << " " << data[j] << " " << observablerror[i][j] << endl;
				}
				//likelihood=likelihood/NUM_OF_OBS;
				//cout << likelihood << endl;
				// And now that you have the likelihood, calculate the weighted average
				//likelihood=1;
				for(int j=0;j<NUM_OF_PARS;j++){
					double diff=(average[j]-param[i][j]);
					//cout << diff << " ";
					//if(diff<0){ diff=-diff; }
					variance[j]+=diff*diff*likelihood;
				}
				totalweight+=double(likelihood);
			}
		}
	}
	cout << "Variances:" << endl;
	for(int j=0;j<NUM_OF_PARS;j++){
		variance[j]=variance[j]/totalweight;
		cout << variance[j] << endl;
	}

	//Clean up the mess
	command = "rm " + EmOutputFile + " " + EmInputFile + " " + tempfilename;
	int result = system(command.c_str());
}