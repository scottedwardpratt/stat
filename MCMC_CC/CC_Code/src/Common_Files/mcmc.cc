#ifndef __MCMC_CC__
#define __MCMC_CC__

#include "mcmc.h"

using namespace std;

MCMC::MCMC(string info_dir, string configuration){
	configname = configuration;
	dir_name = info_dir;
	parameterfile = info_dir+"/defaultpars/";
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
	
	LoadParams();

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

	randnum = new CRandom(1234);
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

MCMC::FirstPass(){
	tracedir = dir_name + "/trace/" + configname;
	
	MAXITERATIONS = parameter::getI(parmap, "MAX_ITERATIONS", 500);
	WRITEOUT = parameter::getI(parmap, "WRITEOUT", 100);
	BURN_IN = parameter::getI(parmap, "BURN_IN", 0);
	RANDOM_THETA0 = parameter::getB(parmap, "RANDOM_THETA0", false);
	VIZTRACE = parameter::getB(parmap, "VISUALIZE_TRACE", true);
	QUIET = parameter::getB(parmap, "QUIET", false);
	
	if(RANDOM_THETA0){
		ThetaList = new ParameterSetList(this, time(NULL), PRESCALED_PARAMS);
	}
	else{
		ThetaList = new ParameterSetList(this);
	}
	//Likelihood_Current=0;

	if(CREATE_TRACE){
		Visualizer = new VizHandler(this);
		Viz_Count = parameter::getI(parmap, "VIZ_COUNT", floor(MAXITERATIONS/200));
		//Visualizer->UpdateTraceFig();
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
	/*if(!QUIET){
		printf("Iteration\tAlpha\tResult\n");
	}*/
}

/** This runs MAXITERATIONS samplings */
MCMC::Run(){
	srand(time(NULL));
	double Likelihood_Current, Likelihood_New;
	double Prior_Current, Prior_New;
	double Proposal_Current, Proposal_New;
	float Scale_Current, Scale_New;

	double LOGBF,alpha;
	ParameterSet * ThetaZeroPtr = ThetaList->Theta[0];
	ParameterSet CurrentParameters = *ThetaZeroPtr;
	
	if(!ThetaZeroPtr){
		cout << "Error getting zeroth iteration parameters." << endl;
	}
	
	Likelihood_Current = Likelihood->Evaluate(*ThetaZeroPtr);
	Scale_Current = (rand() / double(RAND_MAX));
	//Proposal_Current = Proposal->Evaluate(*ThetaZeroPtr);
	Prior_Current = Prior->Evaluate(*ThetaZeroPtr);

	for(int i = 0; i < ThetaList->ParamNames.size(); i++){
		ParamValues.push_back(0);
	}

	Accept_Count = 0;
	for(int i =1; i<=MAXITERATIONS; i++){
		if(LOGLIKE){
			LOGBF = 0;
		}else{
			LOGBF = 1;
		}
		Scale_New = (rand() / double(RAND_MAX));
		ParameterSet Temp_Theta = Proposal->Iterate(CurrentParameters,Scale_New);
		Likelihood_New = Likelihood->Evaluate(Temp_Theta);
		if(i==1){
			bestlikelihood=Likelihood_New;
			BestParameterSetPtr=&CurrentParameters;
		}
		Prior_New = Prior->Evaluate(Temp_Theta);
		Proposal_New = Proposal->Evaluate(CurrentParameters,Temp_Theta,Scale_Current);
		Proposal_Current = Proposal->Evaluate(Temp_Theta,CurrentParameters,Scale_New);
		
		// cout << "Likelihood of proposed set: " << Likelihood_New << endl;
		// cout << "Likelihood of current set: " << Likelihood_Current << endl;
		// cout << "Prior of proposed set: " << Prior_New << endl;
		// cout << "Prior of current set: " << Prior_Current << endl;
		// cout << "Proposal of proposed set: " << Proposal_New << endl;
		// cout << " Proposal of current set: " << Proposal_Current << endl;
		

		// Likelihood
		if(LOGLIKE){
			if(!QUIET){
				printf(" ll_new=%g, ll_current=%g\n",Likelihood_New,Likelihood_Current);
			}
			LOGBF += Likelihood_New-Likelihood_Current;
		} else {
			if(!QUIET){
				printf(" l_new=%g, l_current=%g\n",exp(Likelihood_New),exp(Likelihood_Current));
			}
			LOGBF *= exp(Likelihood_New)/exp(Likelihood_Current);
		}

		// Prior
		/*if(LOGPRIOR){
			LOGBF += (log(Prior_New)-log(Prior_Current));
		} else {
			LOGBF *= (Prior_New/Prior_Current);
		}

		if(!QUIET){
			printf(" Prior_New=%g, Prior_Current=%g\n",Prior_New,Prior_Current);
		}*/

		// Proposal		
		if(LOGLIKE){
			LOGBF += (log(Proposal_Current)-log(Proposal_New));
		} else {
			LOGBF *= Proposal_Current/Proposal_New;
		}

		if(!QUIET){
			printf(" Proposal_New=%g, Proposal_Current=%g\n",Proposal_New,Proposal_Current);
		}
		
		if(LOGLIKE){
			alpha = min(1.0,exp(LOGBF));
		} else {
			alpha = min(1.0,LOGBF);
		}

		//cout << LOGBF << endl;

		if(!QUIET){
			printf("%5d\talpha=%6.5f\t",i,alpha);
			//printf("LOGBF=%6.5f\t",alpha);
		}
		if(alpha > (randnum->ran())) { //Accept the proposed set.
			if(!QUIET){
				printf("Accept\n");
			}
			if(i > BURN_IN){
				Accept_Count++;
			}
			Likelihood_Current = Likelihood_New;
			Prior_Current = Prior_New;
			//Proposal_Current = Proposal_New;
			CurrentParameters = Temp_Theta;
			Scale_Current = Scale_New;
			if(Likelihood_Current>bestlikelihood && i>1){
				bestlikelihood=Likelihood_New;
				BestParameterSetPtr=&CurrentParameters;
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
			if(!QUIET){
				printf("Reject\n");
			}
		}

		if(i > BURN_IN){ // We are just tossing everything in the burn in period.
			ThetaList->Add(CurrentParameters);
		}

		for(int k = 0; k < ThetaList->ParamNames.size(); k++){ //These ParamValus are used for the density plots
			//cout << ThetaList->ParamNames[k] << " " << CurrentParameters.Values[k] << endl;
			//cout << "(" << CurrentParameters.Values[k] << " - " << Min_Ranges[k] << ") / (" << Max_Ranges[k] << " - " << Min_Ranges[k] << ")" << endl; cout.flush();
			//cout << Max_Ranges[k] << " " << Min_Ranges[k] << endl;
			ParamValues[k] = (CurrentParameters.Values[k] - Min_Ranges[k])/(Max_Ranges[k]-Min_Ranges[k]);
			//cout << ParamValues[k] << endl;
		}

		/*if((i > BURN_IN) && (CREATE_TRACE)){
			if((i+1) % Viz_Count == 0){
				Visualizer->UpdateTraceFig();
			}
		}*/

		if((i > BURN_IN) && ((i+1) % WRITEOUT == 0)){
			cout << "Writing out." << endl;
			if(CREATE_TRACE &&(i!=1)){
				Visualizer->UpdateTraceFig();
			}
			ThetaList->WriteOut();
		}
	}
	
	ThetaList->WriteOut();
	ThetaList->MakeTrace();
	if(CREATE_TRACE){
		Visualizer->FinalTrace();
	}
	double ratio = (double)Accept_Count/((double)MAXITERATIONS-(double)BURN_IN);
	cout << "Accepts: " << Accept_Count << endl;
	cout << "Iterations-Burn in: " << MAXITERATIONS-BURN_IN << endl;
	cout << "Acceptance ratio: " << ratio << endl;
	printf("-------- Best Parameter Set, likelihood=%g -------------\n",bestlikelihood);
	BestParameterSetPtr->Print();
}

MCMC::LoadParams(){
	if(EmulatorParams!=""){
		string line;
		stringstream Emparams;
		Emparams << EmulatorParams;
		getline(Emparams,line,' ');
		ParamNames.push_back(line);
		while(!Emparams.eof()){
			getline(Emparams,line,' ');
			if(line.compare(ParamNames.back())!=0) {ParamNames.push_back(line);}
		} 
		if(ParamNames.back().compare(0,1," ")==0 || ParamNames.back().empty()){ // if for some reason the last element is empty, drop it
			ParamNames.pop_back();
		}
	}
	for(int i=0;i<ParamNames.size();i++) {
		printf("%s ",ParamNames[i].c_str());
	}
	printf("\n");
}

MCMC::GetRanges(){
	string filename = info_dir + "/ranges.dat";
	fstream ranges;
	ranges.open(filename.c_str(),fstream::in);
	if(ranges){
		string temps,name,line,type;
		while(ranges >> type){
			if(strcmp(type.c_str(), "double") == 0){
				ranges >> name;
				for(int i=0;i<ParamNames.size();i++){
					if(strcmp(ParamNames[i].c_str(),name.c_str())==0){
						ranges >> Min_Ranges[i]; //minimum
						ranges >> Max_Ranges[i]; //maximum
						if(Min_Ranges[i] > Max_Ranges[i]){ //Flip them
							double temp2 = Min_Ranges[i];
							Min_Ranges[i] = Max_Ranges[i];
							Max_Ranges[i] = temp2;
						}
					}
				}
			}else{
				if(strncmp(type.c_str(), "#", 1) == 0){
					string temp;
					getline(ranges, temp, '\n');
				}
			}
			EmulatorParams = EmulatorParams + name + " ";
			index++;
		}

		ranges.close();
	}
	else{
		cout << "Ranges.dat wont open" << endl;
		exit(1);
	}

	cout << "Ranges loaded" << endl;
}

MCMC::LoadObservables(){
	string observables_filename = info_dir + "/pcanames.dat";
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
			line.erase(line.end()-1,line.end());
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

#endif
