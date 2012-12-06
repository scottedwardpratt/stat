#ifndef __MCMC_CC__
#define __MCMC_CC__

#include "mcmc.h"

using namespace std;

MCMC::MCMC(string info_dir, string configuration){
	srand(time(NULL));
	configname = configuration;
	dir_name = info_dir;
	parameterfile = dir_name+"/defaultpars/";
	cout << "There should be something here ->" << parameterfile << "<-" << endl;
	cout << "In config: " << parameterfile << endl;
	parameter_file_name = parameterfile + "/mcmc.param";
	cout << "Reading in " << parameter_file_name << endl;
	
	parameter::ReadParsFromFile(parmap, parameter_file_name.c_str());
	LOGLIKE = parameter::getB(parmap, "LOGLIKE", true);
	LOGPRIOR = parameter::getB(parmap, "LOGPRIOR", true);
	LOGPROPOSAL = parameter::getB(parmap, "LOGPROPOSAL", true);
	CREATE_TRACE = parameter::getB(parmap, "CREATE_TRACE", true);
	APPEND_TRACE = parameter::getB(parmap, "APPEND_TRACE", false);
	RESCALED_TRACE = parameter::getB(parmap, "RESCALED_TRACE", false);
	SUPPRESS_ERRORS = parameter::getB(parmap, "SUPPRESS_ERRORS", false);
	MODEL = parameter::getS(parmap,"MODEL","NOMODEL");
	PRESCALED_PARAMS = parameter::getB(parmap, "PRESCALED_PARAMS", false);
	WAIT_TO_PRINT_DENSITY_PLOTS = parameter::getB(parmap, "WAIT_TO_PRINT_DENSITY_PLOTS", false);
	MOVING_WINDOW = parameter::getB(parmap, "MOVING_WINDOW", true);
	DENSITY_PLOTS = parameter::getB(parmap, "DENSITY_PLOTS", true);
	
	GetRanges();

	if(PRESCALED_PARAMS){
        cout << "WARNING: PRESCALED_PARAMS = True. This means all parameter ranges are forced to be [0,1]" << endl;
		for(int i=0; i < ParamNames.size(); i++){
			//They go from zero to 1
			Min_Ranges[i]=0;
			Max_Ranges[i]=1;
		}
	}

	LoadObservables();
	
	vector<string> temp_logparam = parameter::getVS(parmap, "LOG_PARAMETERS", "");
	for(int i =0; i<temp_logparam.size(); i++){
		if(strcmp(temp_logparam[i].c_str(), "true") == 0 || strcmp(temp_logparam[i].c_str(), "True") == 0){
			LogParam.push_back(true);
		}else if(strcmp(temp_logparam[i].c_str(), "false") == 0 || strcmp(temp_logparam[i].c_str(), "False") == 0){
			LogParam.push_back(false);
		}else{
			cout << "Unrecognized LogParam value " << temp_logparam[i] << endl;
			exit(1);
		}
	}

    /*cout << "_____________________" << endl;
	cout << "parameterfile: " << parameterfile << endl;
	cout << "MODEL: " << MODEL << endl;
	cout << "dir_name: " << dir_name << endl;
	cout << "parameterfile: " << parameterfile << endl;
	cout << "configname: " << configname << endl;
	cout << "parameter_file_name: " << parameter_file_name << endl;
	cout << "EmulatorParams: " << EmulatorParams << endl;
	cout << "---------------------" << endl;*/

	randnum = new CRandom(time(NULL));
	if(strcmp(MODEL.c_str(),"CosmoSurvey")==0){
		Likelihood = new LikelihoodDistribution_Cosmo(this);
		Prior = new PriorDistribution_Cosmo(this);
	}
	else if(strcmp(MODEL.c_str(),"RHIC")==0){
		Likelihood = new LikelihoodDistribution_RHIC(this);
		Prior = new PriorDistribution_RHIC(this);
	}
	else if(strcmp(MODEL.c_str(),"INTERPOLATOR")==0){
		Likelihood = new LikelihoodDistribution_Interpolator(this);
		Prior = new PriorDistribution_Interpolator(this);
		PRESCALED_PARAMS = false;
	}
	else if(strcmp(MODEL.c_str(),"RHIC_PCA")==0){
		Likelihood = new LikelihoodDistribution_RHIC_PCA(this);
		Prior = new PriorDistribution_RHIC_PCA(this);
	}
	else if(strcmp(MODEL.c_str(),"TEST")==0){
		Likelihood = new LikelihoodDistribution_Test(this);
		Prior = new PriorDistribution_Test(this);
	}
	else{
		printf("Must define parameter MODEL in parameter file, or yours is unrecognized\n");
		exit(1);
	}

	Proposal = new ProposalDistribution(this);

	cout << "Parameter Ranges:" << endl;
	for( int i =0; i< ParamNames.size(); i++){
		cout << Min_Ranges[i] << " " << Max_Ranges[i] << endl;
	}
	
}

MCMC::~MCMC(){

}

void MCMC::FirstPass(){
	tracedir = dir_name + "/trace/" + configname;
	
	MAXITERATIONS = parameter::getI(parmap, "MAX_ITERATIONS", 500);
	VIZOUT = parameter::getI(parmap, "VIZOUT", 100);
	WRITEOUT = parameter::getI(parmap, "WRITEOUT", 100);
	BURN_IN = parameter::getI(parmap, "BURN_IN", 0);
	RANDOM_THETA0 = parameter::getB(parmap, "RANDOM_THETA0", false);
	VIZTRACE = parameter::getB(parmap, "VISUALIZE_TRACE", true);
	QUIET = parameter::getB(parmap, "QUIET", false);

	WriteOutCounter = 0;
	VizWriteOutCounter = 0;
	
	if(RANDOM_THETA0){
		cout << "We are using random theta0 values. They are:" << endl;
		for(int i=0;i<ParamNames.size();i++){
			Theta.push_back((double(rand())/double(RAND_MAX))*(Max_Ranges[i]-Min_Ranges[i])+Min_Ranges[i]);
			cout << ParamNames[i] << " " << Theta[i] << endl;
		}
		ThetaList.push_back(Theta);
	}
	else{
		parameterMap theta_parmap;
		string theta0_filename = parameterfile + "/theta0.param";
		cout << "We are reading theta0 from " << theta0_filename << " The values are: " << endl;
		parameter::ReadParsFromFile(theta_parmap, theta0_filename);
		vector<string> temp_names = parameter::getVS(theta_parmap, "NAMES", "");
		vector<double> temp_values = parameter::getV(theta_parmap, "VALUES", "");
		for(int i=0;i<temp_names.size();i++){
			if(strcmp(ParamNames[i].c_str(),temp_names[i].c_str())==0){
				if((temp_values[i]>Max_Ranges[i])||(temp_values[i]<Min_Ranges[i])){
					cout << "The theta0 value given for " << ParamNames[i] << " in " << theta0_filename << " is not within the ranges specified for this parameter. " << endl;
					cout << "The range specifed was: [" << Min_Ranges[i] << "," << Max_Ranges[i] << "] but the theta0 value given was: " << temp_values[i] << endl;
					exit(-1);
				}
				Theta.push_back(temp_values[i]);
				cout << ParamNames[i] << " " << temp_values[i] << endl;
			}
		}
		if(Theta.size()!=ParamNames.size()){
			cout << "The theta0 vector specifed in " << theta0_filename << " is not the right length, or the names are incorrect." << endl;
			cout << "Only " << Theta.size() << " values were read in, but there are " << ParamNames.size() << " parameters." << endl;
			exit(-1);
		}
		ThetaList.push_back(Theta);
	}
	//Likelihood_Current=0;

	if(CREATE_TRACE){
		Visualizer = new VizHandler(this);
	}

	if(APPEND_TRACE){
		string addon = "";
		bool Done = false;
		int filecount = 0;
		while(!Done){
			struct stat st;
			stringstream ss;
			string tempfile = tracedir + addon;
			if(stat(tempfile.c_str(), &st)==0){
				//directory exists.
				cout << tempfile << " exists, trying next option..." << endl;
				filecount++;
				ss << "_" << filecount;
				addon = ss.str();
				ss.str(string());
			}else{
				//doesn't exist
				Done = true;
				tracedir = tempfile;
			}
		}
	}else{
		cout << "Deleting prior trace data." << endl;
		string cmd = "rm " + tracedir + "/output*.dat " + tracedir + "/trace.dat";
		system(cmd.c_str());
	}
	
	string command = "mkdir -p "+ tracedir;
	
	system(command.c_str());
}

/** This runs MAXITERATIONS samplings */
void MCMC::Run(){
	double Likelihood_Current, Likelihood_New;
	double Prior_Current, Prior_New;
	double Proposal_Current, Proposal_New;
	float Scale_Current, Scale_New;

	double LOGBF,alpha;
	Proposed_Theta = Theta;
	
	Likelihood_Current = Likelihood->Evaluate(Theta);
	Scale_Current = (rand() / double(RAND_MAX));
	Prior_Current = Prior->Evaluate(Theta);

	Accept_Count = 0;
	for(int i =1; i<=MAXITERATIONS; i++){
		if(LOGLIKE){
			LOGBF = 0;
		}else{
			LOGBF = 1;
		}

		Scale_New = (rand() / double(RAND_MAX));
		Proposed_Theta = Proposal->Iterate(Theta,Scale_New);
		//JFN 12/5/12: right from Scott's code
		/*for(int k=0; k<ParamNames.size();k++){
			Proposed_Theta[k] = Theta[k]+0.1*randnum->gauss();
			while((Proposed_Theta[k]>Max_Ranges[k])||(Proposed_Theta[k]<Min_Ranges[k])){
				if(Proposed_Theta[k]>Max_Ranges[k]) Proposed_Theta[k]=2*Max_Ranges[k]-Proposed_Theta[k];
				if(Proposed_Theta[k]<Min_Ranges[k]) Proposed_Theta[k]=2*Min_Ranges[k]-Proposed_Theta[k];
			}
		}*/

		Likelihood_New = Likelihood->Evaluate(Proposed_Theta);
		if(i==1){
			bestlikelihood = Likelihood_New;
			BestParameterSet = Theta;
			Scaled_Theta = Theta;
		}
		Prior_New = Prior->Evaluate(Proposed_Theta);
		Proposal_New = Proposal->Evaluate(Theta,Proposed_Theta,Scale_Current);
		Proposal_Current = Proposal->Evaluate(Proposed_Theta,Theta,Scale_New);
		
		// cout << "Likelihood of proposed set: " << Likelihood_New << endl;
		// cout << "Likelihood of current set: " << Likelihood_Current << endl;
		// cout << "Prior of proposed set: " << Prior_New << endl;
		// cout << "Prior of current set: " << Prior_Current << endl;
		// cout << "Proposal of proposed set: " << Proposal_New << endl;
		// cout << " Proposal of current set: " << Proposal_Current << endl;
		

		// Likelihood
		if(LOGLIKE){
			if(!QUIET) printf(" ll_new=%g, ll_current=%g\n",Likelihood_New,Likelihood_Current);
			LOGBF += Likelihood_New-Likelihood_Current;
		} else {
			if(!QUIET) printf(" l_new=%g, l_current=%g\n",exp(Likelihood_New),exp(Likelihood_Current));
			LOGBF *= exp(Likelihood_New)/exp(Likelihood_Current);
		}

		// Prior
		/*if(LOGPRIOR){
			LOGBF += (log(Prior_New)-log(Prior_Current));
		} else {
			LOGBF *= (Prior_New/Prior_Current);
		}

		if(!QUIET) printf(" Prior_New=%g, Prior_Current=%g\n",Prior_New,Prior_Current);*/

		// Proposal		
		if(LOGLIKE){
			LOGBF += (log(Proposal_Current)-log(Proposal_New));
		} else {
			LOGBF *= Proposal_Current/Proposal_New;
		}

		if(!QUIET) printf(" Proposal_New=%g, Proposal_Current=%g\n",Proposal_New,Proposal_Current);
		
		if(LOGLIKE){
			alpha = min(1.0,exp(LOGBF));
		} else {
			alpha = min(1.0,LOGBF);
		}

		//cout << LOGBF << endl;

		if(!QUIET) printf("%5d\talpha=%6.5f\t",i,alpha);
		if(alpha > (randnum->ran())) { //Accept the proposed set.
		//double ll = Likelihood_New;
		//double oldll = Likelihood_Current;
		//if((ll>oldll || randnum->ran()<exp(ll-oldll)) && (ll-oldll>-40)){
		//if(ll>oldll || randnum->ran()<exp(ll-oldll)){
			if(!QUIET) printf("Accept\n");
			if(i > BURN_IN){
				Accept_Count++;
			}
			Likelihood_Current = Likelihood_New;
			Prior_Current = Prior_New;
			//Proposal_Current = Proposal_New;
			Theta = Proposed_Theta;
			Scale_Current = Scale_New;
			if(Likelihood_Current>bestlikelihood && i>1){
				bestlikelihood = Likelihood_New;
				BestParameterSet = Theta;
				if(!QUIET){
					if(LOGLIKE){
						printf("XXXXXXXXX YIPPEE!! Best parameters so far, loglikelihood=%g\n",bestlikelihood);
					}
					else{
						printf("XXXXXXXXX YIPPEE!! Best parameters so far, likelihood=%g\n",bestlikelihood);
					}
				}
			}
		}else{
			if(!QUIET) printf("Reject\n");
		}

		for(int k = 0; k < ParamNames.size(); k++){ //These ParamValus are used for the density plots
			Scaled_Theta[k] = (Theta[k] - Min_Ranges[k])/(Max_Ranges[k]-Min_Ranges[k]);
		}

		if(i >= BURN_IN){ // We are just tossing everything in the burn in period.
			ThetaList.push_back(Theta);
			Scaled_ThetaList.push_back(Scaled_Theta);
			VizThetaList.push_back(Theta);
			VizScaled_ThetaList.push_back(Scaled_Theta);
			if(CREATE_TRACE && ((i+1) % VIZOUT == 0)){
				Visualizer->UpdateTraceFig();
				VizThetaList.clear();
				VizScaled_ThetaList.clear();
			}
			if((i+1) % WRITEOUT == 0){
				if(QUIET) cout << "."; cout.flush();
				PrintDataToFile();
				ThetaList.clear();
				Scaled_ThetaList.clear();
			}
		}
	}
	
	PrintDataToFile();
	MakeFinalTrace();

	if(CREATE_TRACE){
		Visualizer->FinalTrace();
	}
	double ratio = (double)Accept_Count/((double)MAXITERATIONS-(double)BURN_IN);
	cout << "Accepts: " << Accept_Count << endl;
	cout << "Iterations-Burn in: " << MAXITERATIONS-BURN_IN << endl;
	cout << "Acceptance ratio: " << ratio << endl;
	printf("-------- Best Parameter Set, likelihood=%g -------------\n",bestlikelihood);
	for(int i=0;i<ParamNames.size();i++){
		cout << ParamNames[i] << " " << BestParameterSet[i] << endl;
	}
}

void MCMC::GetRanges(){
	string filename = dir_name + "/ranges.dat";
	fstream ranges;
	ranges.open(filename.c_str(),fstream::in);
	if(ranges){
		string temps,name,line,type;
		int index = 0;
		while(ranges >> type){
			if(strcmp(type.c_str(), "double") == 0){
				ranges >> name;
				ranges >> Min_Ranges[index]; //minimum
				ranges >> Max_Ranges[index]; //maximum
				if(Min_Ranges[index] > Max_Ranges[index]){ //Flip them
					double temp2 = Min_Ranges[index];
					Min_Ranges[index] = Max_Ranges[index];
					Max_Ranges[index] = temp2;
				}
			}else{
				if(strncmp(type.c_str(), "#", 1) == 0){
					string temp;
					getline(ranges, temp, '\n');
				}
			}
			ParamNames.push_back(name);
			EmulatorParams = EmulatorParams + name + " ";
		}

		ranges.close();
	}
	else{
		cout << "Ranges.dat wont open" << endl;
		exit(1);
	}

	cout << "Ranges and parameters loaded" << endl;

	for(int i=0;i<ParamNames.size();i++) {
		printf("%s ",ParamNames[i].c_str());
	}
	printf("\n");
}

void MCMC::LoadObservables(){
	string observables_filename = dir_name + "/pcanames.dat";
	fstream PcaNames;
	PcaNames.open(observables_filename.c_str(),fstream::in);

	if(PcaNames){
		string line,name;
		getline(PcaNames,line,'\n');
		while(!PcaNames.eof()){
			while(line.compare(0,1,"#") == 0 || line.empty() ) { // keep reading until it's not a comment
				getline(PcaNames,line,'\n'); 
			} 
	
			stringstream Observableslist;
			//line.erase(line.end()-1,line.end());
			Observableslist << line;
			
			while(!Observableslist.eof()){
				getline(Observableslist,name,' ');
				if(!name.empty() || name!=" ")
				ObservablesNames.push_back(name);
				//cout << "Observable: " << ObservablesNames.back() << endl;
			}
			getline(PcaNames,line,'\n');
		}
		if(ObservablesNames.back().compare(0,1," ")==0 || ObservablesNames.back().empty()){ // if for some reason the last element is empty, drop it
			ObservablesNames.pop_back();
		}
	}
	else{
		cout << "Could not open " << observables_filename.c_str() << endl;
		exit(1);
	}
	PcaNames.close();
	
	cout << "The emulated observables are: " << endl;
	for(int i=0;i<ObservablesNames.size();i++) {
		cout << ObservablesNames[i] << endl;
	}
}

void MCMC::PrintDataToFile(){
	stringstream ss;
	ss << WriteOutCounter+1;
	string filename = tracedir + "/output"+ ss.str() +".dat";
	
	ofstream outputfile;
	outputfile.open(filename.c_str());

	if(!QUIET){
		cout << "Printing data to file." << endl;
		cout << "Writing out to: " << filename << endl;
	}

	if(outputfile){
		outputfile << "#ITERATION,";
		for(int i = 0; i<ParamNames.size(); i++){
			outputfile << ParamNames[i] << ",";
		}
		// outputfile << "LIKELIHOOD" << endl;

		for(int i =0; i < WRITEOUT; i++){
			outputfile << i+WriteOutCounter*WRITEOUT << ",";
			for(int j=0; j < ParamNames.size(); j++){
				if(RESCALED_TRACE){
					outputfile << (Scaled_ThetaList[i][j]) << ",";
				}else{
					outputfile << ThetaList[i][j] << ",";
				}
			}
			outputfile << endl;
			// outputfile << "," << LikelihoodArray[i] << endl;
		}
		outputfile.close();
	}else{
		cout << "Error in writing output file." << endl;
		exit(1);
	}
	WriteOutCounter++;

	if(!QUIET) cout << "Done printing to file." << endl;
}

void MCMC::MakeFinalTrace(){
	stringstream ss;
	ss << "cat ";
	
	for(int i = 1; i <=ceil((double)(MAXITERATIONS)/(double)(WRITEOUT)); i++){
		cout << "Parsing output" << i << ".dat" << endl;
		ss << tracedir << "/output" << i << ".dat ";
	}
	ss << "> " << tracedir << "/trace.dat" << endl;
	
	string command = ss.str();
	system(command.c_str());
	
	ss.str(string());
}

#endif
