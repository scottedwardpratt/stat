#ifndef __LIKELIHOOD_CC__
#define __LIKELIHOOD_CC__

#include "Model.h"
#include <time.h>

madai::LikelihoodDistribution_Test::LikelihoodDistribution_Test(madai::Model *in_Model){
  m_Model = in_Model;
	m_SepMap = parameter::getB(m_Model->m_ParameterMap, "LIKELIHOOD_PARAMETER_MAP", false);
	
	if(m_SepMap){
		std::string parmapfile = m_Model->m_ParameterFile + "/likelihood.param";
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
		std::cout << "This is a guassian test. The UseEmulator flag will be ignored." << std::endl;
	}
	else{
		exit(1);
	}

	m_Data = GetData();

	//testing the outputs of the emulator at various points	// 
	m_EmulatorTest.open("PCA0.dat");
	m_EmulatorTest.close();
	
}

madai::LikelihoodDistribution_Test::~LikelihoodDistribution_Test(){
}

double madai::LikelihoodDistribution_Test::Evaluate(std::vector<double> Theta){
	clock_t begintime;
	std::vector<double> ModelMeans;
	std::vector<double> ModelErrors;
	double likelihood;
	
	if(m_Timing){
		begintime = clock();
	}

	//Fill vectors:
	if(m_UseEmulator){
		// This is a gaussian test, so there shouldn't be an 'emulator' per say.
		//emulator->QueryEmulator(Theta, ModelMeans, ModelErrors); //fills vectors with emulator output
	}
	else{
		int errors[]={0,0,0,0}; //For now
		ModelErrors.assign (errors,errors+4);
		int means[4];
		for(int i=0;i<4;i++){
			means[i]=Theta[i];
		}
		ModelMeans.assign(means,means+4);
	}
	
	
	//Initialize GSL containers
	int N = ModelErrors.size();
	//gsl_matrix * sigma = gsl_matrix_calloc(N,N);
	gsl_vector * model = gsl_vector_alloc(N);
	gsl_vector * mu = gsl_vector_alloc(N);
	// cout << "Done allocating gsl containers." << endl;
	
	
	//Read in appropriate elements
	for(int i = 0; i<N; i++){
		// gsl_matrix_set(sigma, i,i,Theta.GetValue("SIGMA"));
		//gsl_matrix_set(sigma,i,i,ModelErrors[i]);
		gsl_vector_set(model, i,ModelMeans[i]);
		gsl_vector_set(mu, i, m_Data[i]);
	}
	
	//likelihood = Log_MVNormal(*model, *mu, *sigma);
	
	for(int i=0;i<4;i++){
		likelihood+=(ModelMeans[i]-m_Data[i])*(ModelMeans[i]-m_Data[i])/2; //This is a gaussian, but it would be "log(exp(stuff))", so I just left it "stuff"
	}

	if(!(m_Model->m_LogLike)){
		likelihood = exp(likelihood);
	}
	
	if(m_Verbose){
		/*double sum = 0.0;
		
		for(int i = 0; i< N; i++){
			sum += (gsl_vector_get(model, i) - gsl_vector_get(mu, i));
		}
		sum = sum/(double)N;
		cout << "Average difference between outputs:" << sum << endl;*/
	}
	
	//deallocate GSL containers.
	gsl_vector_free(model);
	gsl_vector_free(mu);
	//gsl_matrix_free(sigma);
	
	if(m_Timing){
		std::cout << "Likelihood evaluation took " << (clock()-begintime)*1000/CLOCKS_PER_SEC << " ms." << std::endl;
	}
	
	std::cout << "PCA 0: " << ModelMeans[0] << std::endl;
	
	m_EmulatorTest.open("PCA0.dat", ios_base::app);
	m_EmulatorTest << ModelMeans[0] << std::endl;
	m_EmulatorTest.close();
	// emulator_test << ModelMeans[0] << endl;
	
	return likelihood;
}

std::vector<double> madai::LikelihoodDistribution_Test::GetData(){
	//Four separate gaussians with means of 1, 2, 3, and 4 and standard deviations of 2. The means will be the input, 
	//and three values on the gaussian will be the output.

	std::vector<double> datameans;
	int mymeans[]={1,2,3,4};
	datameans.assign (mymeans,mymeans+4);
	//vector<double> dataerror;
	/*
	parameterMap actualparmap;
	
	string actual_filename = mcmc->parameterfile + "/actual.param";
	parameter::ReadParsFromFile(actualparmap, actual_filename);
	
	vector<string> temp_names = parameter::getVS(actualparmap, "NAMES", "");
	vector<double> temp_values = parameter::getV(actualparmap, "VALUES", "");
	
	ParameterSet ActualParams;
	ActualParams.Initialize(temp_names, temp_values);
	
	emulator->QueryEmulator(ActualParams, datameans, dataerror);*/

	return datameans;
}
#endif