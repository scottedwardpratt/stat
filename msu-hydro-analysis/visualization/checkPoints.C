// This code is for comparing multiple emulator points with test points

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include <math.h>
#include <vector>

using namespace std;

int main(int argc, char *argv[]){
	if(argc < 2){
		printf("Usage: ./viz/checkPoints model_directory\n");
		printf("Example: ./viz/checkPoints model_directory/Aprin2012\n");
		printf("the code expects that the paramters will be in model_directory/parameters/run*\n");
		printf("and that the model values will be in model_directory/model_results/run*\n");
		exit(-1);
	}

	int NUM_OF_PARS = 6;
	int NUM_OF_OBS = 30;
	ifstream trace,inputfile,tempfile;
	ofstream outputfile;
	string line,filename,EmInputFile,EmOutputFile,temp,temp2,tempfilename;
	int index;
	double param[NUM_OF_PARS],average[NUM_OF_PARS],variance[NUM_OF_PARS];
	double data[NUM_OF_OBS],dataerror[NUM_OF_OBS],observable[NUM_OF_OBS],observablerror[NUM_OF_OBS];
	vector<string> datanames;
	char buffer[50];
	vector<string> names;

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

	string model_directory(argv[1]);
	//cout << dirname << endl;
	string command = "ls " + model_directory + "/model_results/ > " + tempfilename;
	int hold = system(command.c_str());
	tempfile.open(tempfilename.c_str());
	getline(tempfile,filename,'\n');
	getline(tempfile,filename,'\n');

	bool firsttime=true;

	//while(0<1){
	while(strcmp(filename.c_str(),"")!=0){
		cout << filename << endl;
	//while(!(filename==lastfile)){ // Read each file

		// First, read in the real data =========================================
		if(firsttime){
			ifstream datafile;
			string exp_data;
				exp_data = model_directory + "/model_results/" + filename + "/results.dat"; // At the moment there is no reading of parameter lists to make sure things line up, so make a trimmed list with only the used observables in order
				datafile.open(exp_data.c_str());
				getline(datafile,line,'\n');
				index = 0;

				while(!datafile.eof()){
					stringstream ss;
					ss << line;
					ss >> temp >> temp2 >> data[index] >> dataerror[index];
					//cout << data[index] << " " << dataerror[index] << endl;
					getline(datafile,line,'\n');
					if(firsttime){
						datanames.push_back(temp2);
						//cout << temp2 << "|" << endl;
					}
					index++;
				}
				datafile.close();

				string obser_names_file;
				obser_names_file = model_directory + "/standardlong.dat";
				datafile.open(obser_names_file.c_str());
				getline(datafile,line,'\n');
				bool done=false;
				while(!done){
					if(line.compare(0,1,"#") != 0 && !line.empty()){
						done=true;
						stringstream namestream;
						namestream << line;
						string name;
						while(!namestream.eof()){
							getline(namestream,name,' ');
							names.push_back(name);
							//cout << name << "|" << endl;
						}
						names.pop_back();
					}
					getline(datafile,line,'\n');
				}
			}

		// Now we read in the parameters ========================================
		filename = model_directory + "/parameters/" + filename + "/stats.param";
		cout << filename << endl;
		trace.open(filename.c_str());
		getline(trace,line,'\n');

		index=0;
		while(!trace.eof()){
			stringstream ss;
			ss << line;
			ss >> temp >> temp >> param[index] >> temp;
			//cout << line << endl;
			//cout << param[index] << endl;
			//cout << data[index] << " " << dataerror[index] << endl;
			getline(trace,line,'\n');
			index++;
		}

		// Now that you've read the params, make your observables ================================
		outputfile.open(EmInputFile.c_str()); 
		outputfile << "#Not sure, does this line matter?" << endl;

		for(int j=0;j<names.size();j++){
			outputfile << " " << param[j];
		}
		outputfile << endl;
			
		outputfile.close();
		command = "./src/computePoints.sh " + model_directory + " " + model_directory + "/fn-data-standardlong.dat < " + EmInputFile + " > " + EmOutputFile;
		int result = system(command.c_str()); // Caculate the observables from the emulator

		command = "rm " +EmInputFile;
		result = system(command.c_str());

		// Get the calculated observables into memory so we can do stuff with them =============================
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
					for(int j=0;j<names.size();j++){ // Read untill there is nothing more to read
						is >> obsbuff;
						observable[j]=obsbuff;
						//cout << observable[int(linecount/2)][j] << " ";
						//cout << observable[int(linecount/2)][j] << " " << int(linecount/2) <<  " " << linecount << " " << j << endl;
					}
				}else{
					for(int j=0;j<names.size();j++){ // Read untill there is nothing more to read
						is >> obsbuff;
						observablerror[j]=obsbuff;
					}
				}
				linecount++;
			}
			getline(inputfile,line,'\n');
		}
		inputfile.close();

		command = "rm " +EmOutputFile;
		result = system(command.c_str());

		int above,below;
		above=0;
		below=0;

		// Now, do what you will with what you have made ==================================================
		for(int j=0;j<names.size();j++){
			//cout << names[j] << endl;
			for(int k=0;k<datanames.size();k++){
				//cout << "     " << datanames[k] << endl;
				//if(names[j].compare(0,1,datanames[k]) == 0){
				if(strcmp(names[j].c_str(),datanames[k].c_str())==0){
					cout <<  names[j] << " mod " << data[k] << " emu " << observable[j] <<  " emu error " << observablerror[j] << endl;
					cout << "difference " << sqrt((data[k]-observable[j])*(data[k]-observable[j])) << " diff/mod " << sqrt((data[k]-observable[j])*(data[k]-observable[j]))/data[k] << " diff/emu " << sqrt((data[k]-observable[j])*(data[k]-observable[j]))/observable[j] << " error/emu " << observablerror[j]/observable[j] << endl;
				if(sqrt((data[k]-observable[j])*(data[k]-observable[j])) < observablerror[j]) {above++;}
				if(sqrt((data[k]-observable[j])*(data[k]-observable[j])) > observablerror[j]) {below++;}
				}
			}
		}
		cout << "The emulator error was too large " << above << " times, and too small " << below << " times" << endl;
		
		//cout << "================="<< endl;
		getline(trace,line,'\n');
		
		trace.close();
		filename.clear();
		getline(tempfile,filename,'\n');
		firsttime=false;
	}
	tempfile.close();

	//Clean up the mess
	//command = "rm " + EmOutputFile + " " + EmInputFile + " " + tempfilename;
	//int result = system(command.c_str());
}