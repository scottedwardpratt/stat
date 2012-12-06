#include "Distribution.h"
#include <time.h>

madai::LikelihoodDistribution_Cosmo::LikelihoodDistribution_Cosmo(madai::Model *in_Model){
  m_Model = in_Model;
	m_SepMap = parameter::getB(m_Model->m_ParameterMap, "LIKELIHOOD_PARAMETER_MAP", false);
	
	if(m_SepMap){
		std::string parmapfile = m_Model->m_DirectoryName + "defaultpars/likelihood.param";
		m_ParameterMap = new parameterMap;
		parameter::ReadParsFromFile(*m_ParameterMap, parmapfile);
	}else{
    m_ParameterMap = &(m_Model->m_ParameterMap);
	}
	
	// cout << "Like param map made." << endl;
	
	m_UseEmulator = parameter::getB(*m_ParameterMap, "USE_EMULATOR", false);
	m_Timing = parameter::getB(*m_ParameterMap, "TIMING", false) || parameter::getB(*m_ParameterMap, "TIME_LIKELIHOOD", false);
	m_Verbose = parameter::getB(*m_ParameterMap, "VERBOSE", false) || parameter::getB(*m_ParameterMap, "VERBOSE_LIKELIHOOD", false);
	
	// cout << "params declared." << endl;
	if(m_UseEmulator){
		std::cout << "Turn off USE_EMULATOR" << std::endl;
		exit(-1);
		//emulator = new EmulatorHandler(m_ParameterMap, mcmc_in);
	}

	m_Data = GetData();
	m_intData.resize(m_Data.size());
	for(int i = 1; i < m_intData.size(); i++){
		m_intData[i]=gsl_ran_poisson(m_RandNumGen,m_Data[i]);
	}
}

madai::LikelihoodDistribution_Cosmo::~LikelihoodDistribution_Cosmo(){
	//delete emulator;
}

double madai::LikelihoodDistribution_Cosmo::Evaluate(std::vector<double> Theta){
	clock_t begintime;
	vector<double> ModelMeans;
	vector<double> ModelErrors;
	double likelihood = 0.0,dll;
	std::stringstream ss;
	std::ifstream inputfile;
	
	if(m_Timing){
		begintime = clock();
	}
	
  ss << "cosmosurvey";
  std::vector<madai::Parameter> const * parameters = &(m_Model->GetParameters());
  for(int i=0; i<parameters->size();i++)
    ss << " -" << (*parameters)[i].m_Name << " " << Theta[i];

	ss << " -nz 10 -nf 10 -ob .0406 > output.dat" << std::endl;
	
	std::cout << ss.str() << std::endl;
	std::cout << "Waiting on cosmosurvey...";
	std::cout.flush();
	int result = system((ss.str()).c_str());
	std::cout << "Done." << std::endl;
	
	inputfile.open("output.dat");
	
	while(!inputfile.eof()){
		double temp;
		inputfile >> temp;
		ModelMeans.push_back(temp);
	}
	
	
	for(int i = 1; i < ModelMeans.size(); i++){
		dll=log(gsl_ran_poisson_pdf(m_intData[i],ModelMeans[i]));
			//printf("ModelMeans[%d]=%g, intDATA[%d]=%d, Dloglikelihood=%g\n",i,ModelMeans[i],i,intDATA[i],dll);
		//likelihood += log(gsl_ran_poisson_pdf(static_cast<unsigned int>(ModelMeans[i] + 0.5), DATA[i]));
		likelihood += dll;
	}
	printf("XXXXX LogLikelihood=%g XXXXXX\b",likelihood);
	if(likelihood>-0.0001){
		for(int i = 1; i < ModelMeans.size(); i++){
			dll=log(gsl_ran_poisson_pdf(m_intData[i],ModelMeans[i]));
			printf("ModelMeans[%d]=%g, intDATA[%d]=%d, Dloglikelihood=%g\n",i,ModelMeans[i],i,m_intData[i],dll);
				//likelihood += log(gsl_ran_poisson_pdf(static_cast<unsigned int>(ModelMeans[i] + 0.5), DATA[i]));
		}
		exit(1);
	}
	
	if(!(m_Model->m_LogLike)){
		likelihood = exp(likelihood);
	}
	
	if(m_Timing){
		std::cout << "Likelihood evaluation took " << (clock()-begintime)*1000/CLOCKS_PER_SEC << " ms." << std::endl;
	}
	
	return likelihood;
}

std::vector<double> madai::LikelihoodDistribution_Cosmo::GetData(){
	std::vector<double> datameans;
	std::stringstream ss;
	parameterMap actualparmap;
	std::ifstream inputfile;
	
	std::string actual_filename = m_Model->m_DirectoryName + "/defaultpars/actual.param";
	parameter::ReadParsFromFile(actualparmap, actual_filename);
	
	std::vector<std::string> temp_names = parameter::getVS(actualparmap, "NAMES", "");
	std::vector<double> temp_values = parameter::getV(actualparmap, "VALUES", "");
	
	std::vector<double> ActualParamValues;
  std::vector<std::string> ActualParamNames;
	ActualParamValues = temp_values;
  ActualParamNames = temp_names;
	
	ss << "cosmosurvey";
	for(int i = 0; i < temp_names.size(); i++){
		ss << " -" << temp_names[i] << " " << temp_values[i];
	}
	ss << " -nz 10 -nf 10 -ob .0406 > data_output.dat" << endl;
	
	std::cout << ss.str() << std::endl;
	std::cout << "Waiting on cosmosurvey...";
	std::cout.flush();
	int result = system((ss.str()).c_str());
	std::cout << "Done." << std::endl;
	
	inputfile.open("data_output.dat");
	
	int counter = 0;
	while(!inputfile.eof()){
		counter++;
		double temp;
		inputfile >> temp;
		std::cout << "Data point " << counter << ": " << temp << std::endl;
		datameans.push_back(temp);
	}
	return datameans;
}