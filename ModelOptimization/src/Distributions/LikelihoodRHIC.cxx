#ifndef __LIKELIHOOD_CC__
#define __LIKELIHOOD_CC__

#include "Model.h"
#include <time.h>
#include "EmuPlusPlus/EmuPlusPlus.h"

madai::LikelihoodDistribution_RHIC::LikelihoodDistribution_RHIC(madai::Model *in_Model){
  m_Model = in_Model;
  m_SuppressErrors = parameter::getB(m_Model->m_ParameterMap, "SUPPRESS_ERRORS", false);
  m_SepMap = parameter::getB(m_Model->m_ParameterMap, "LIKELIHOOD_PARAMETER_MAP", false);
  if(m_SepMap){
    std::string parmapfile = m_Model->m_ParameterFile + "/likelihood.param";
    m_ParameterMap = new parameterMap;
    parameter::ReadParsFromFile(*m_ParameterMap,parmapfile);
  }else{
    m_ParameterMap = &(m_Model->m_ParameterMap);
  }

  m_UseEmulator = parameter::getB(*m_ParameterMap, "USE_EMULATOR", false);
  m_Timing = parameter::getB(*m_ParameterMap, "TIMING", false) || parameter::getB(*m_ParameterMap, "TIME_LIKELIHOOD", false);
  m_Verbose = parameter::getB(*m_ParameterMap, "VERBOSE", false) || parameter::getB(*m_ParameterMap, "VERBOSE_LIKELIHOOD", false);
  m_FakeData = parameter::getB(*m_ParameterMap, "FAKE_DATA", false);

  if(m_UseEmulator && !(m_Model->m_ProcessPipe)){
    std::cout << "Emulator is being loaded from: " << m_Model->m_DirectoryName + "/Emulator.statefile" << std::endl;
    m_Emulator = new ::emulator(m_Model->m_DirectoryName + "/Emulator.statefile");
    //emulator m_Emulator(mcmc->dir_name + "/Emulator.statefile");
    std::cout << "Emulator loaded. Test (number of params): _" << m_Emulator->number_params << "_" << std::endl;
  }else if(m_Model->m_ProcessPipe){
    std::cout << "Emulator loaded in MultiModel" << std::endl;
  }else{
    std::cout << "The UseEmulator flag is set to false (or not set). We can't do anything without an emulator" << std::endl;
    exit(1);
  }

  if(m_FakeData){
    m_Data = GetFakeData();
  }else{
    m_Data = GetRealData();
  }

  //ERROR = GetRealError();

  //testing the outputs of the emulator at various points	// 
  //emulator_test.open("PCA0.dat");
  //emulator_test.close();

}

madai::LikelihoodDistribution_RHIC::~LikelihoodDistribution_RHIC(){
  if(m_UseEmulator && !(m_Model->m_ProcessPipe)){
    delete m_Emulator;
  }
}

double madai::LikelihoodDistribution_RHIC::Evaluate(std::vector<double> Theta){
  clock_t begintime;
  std::vector<double> ModelMeans;
  std::vector<double> ModelErrors;
  double likelihood;

  if(m_Timing){
    begintime = clock();
  }

  if(m_UseEmulator && !(m_Model->m_ProcessPipe)){
    m_Emulator->QueryEmulator(Theta, ModelMeans, ModelErrors); //fills vectors with emulator output
  }else if(m_Model->m_ProcessPipe){
    std::vector<double> temp_outs;
    m_Model->GetScalarOutputs(Theta, temp_outs); //fills temp_outs with emulator output
    for(int i = 0; i < int(temp_outs.size()/2); i++){ //have to sort the means and variances
      ModelMeans.push_back(temp_outs[2*i]);
      ModelErrors.push_back(temp_outs[2*i+1]);
    }
  }else{
    //determine another way to fill the vectors
  }

  //Initialize GSL containers
  int N = ModelErrors.size();
  gsl_matrix * sigma = gsl_matrix_calloc(N,N);
  //gsl_matrix * sigma_data = gsl_matrix_calloc(N,N);
  gsl_vector * model = gsl_vector_alloc(N);
  gsl_vector * mu = gsl_vector_alloc(N);
  // cout << "Done allocating gsl containers." << endl;

  if(m_Verbose){
    std::cout << "Theta: ";
    for(int i = 0; i < Theta.size(); i++){
      std::cout << Theta[i] << " ";
    }
    std::cout << std::endl << "Observable, Model means, Model errors, Data" << std::endl;
    for(int i = 0; i<N; i++){
      std::cout << m_Model->GetScalarOutputNames()[i] << " " << ModelMeans[i] << " " << ModelErrors[i] << " " << m_Data[i] << std::endl;
    }
  }
  //Read in appropriate elements
  for(int i = 0; i<N; i++){
    //cout << " Data: " << DATA[i] << " Emu: " << ModelMeans[i] << " +/-: " << ModelErrors[i] << endl;
    if(m_SuppressErrors){
      //ModelErrors[i]=ModelMeans[i]*0.1; //What a reasonable error is depends on the observable
      ModelErrors[i]=1;
    }
    gsl_matrix_set(sigma, i, i,ModelErrors[i]);
    gsl_vector_set(model, i, ModelMeans[i]);
    gsl_vector_set(mu, i, m_Data[i]);
    //cout << "i: " << i << " Data: " << DATA[i] <<  " Mean: " << ModelMeans[i] << " Error: " << ModelErrors[i] << endl;
  }

  likelihood = Log_MVNormal(*model, *mu, *sigma);
  //likelihood = Gaussian(*model, *mu, *sigma);
  //likelihood = Gaussian(*model, *mu, *sigma, *sigma_data); //This is the integrated likelihood

  if(m_Verbose){
    double sum = 0.0;

    for(int i = 0; i< N; i++){
      sum += (gsl_vector_get(model, i) - gsl_vector_get(mu, i));
    }
    sum = sum/(double)N;
    std::cout << "Average difference between outputs:" << sum << std::endl;
  }

  //deallocate GSL containers.
  gsl_vector_free(model);
  gsl_vector_free(mu);
  gsl_matrix_free(sigma);

  if(m_Timing){
    std::cout << "Likelihood evaluation took " << (clock()-begintime)*1000/CLOCKS_PER_SEC << " ms." << std::endl;
  }

  //cout << "PCA 0: " << ModelMeans[0] << endl;

  //emulator_test.open("PCA0.dat", ios_base::app);
  //emulator_test << ModelMeans[0] << endl;
  //emulator_test.close();
  //emulator_test << ModelMeans[0] << endl;

  return likelihood;
}

std::vector<double> madai::LikelihoodDistribution_RHIC::GetFakeData(){
	//This makes some fake results by querying the emulator with the parameters from actual.param
  std::vector<double> datameans;
  std::vector<double> dataerror;
	
	parameterMap actualparmap;
	
  std::string actual_filename = m_Model->m_ParameterFile + "/actual.param";
	parameter::ReadParsFromFile(actualparmap, actual_filename);
	
  std::vector<std::string> temp_names = parameter::getVS(actualparmap, "NAMES", "");
  std::vector<double> temp_values = parameter::getV(actualparmap, "VALUES", "");
	
  if(m_UseEmulator && !(m_Model->m_ProcessPipe)){
    m_Emulator->QueryEmulator(temp_values, datameans, dataerror);
  }else if(m_Model->m_ProcessPipe){
    std::vector<double> temp_outs;
    m_Model->GetScalarOutputs(temp_values, temp_outs);
    for(int i = 0; i < int(temp_outs.size()/2); i++){
      datameans.push_back(temp_outs[2*i]);
      dataerror.push_back(temp_outs[2*i+1]);
    }
	}

  std::cout << "We are using FAKE DATA!!!!!!!!" << std::endl;
  std::cout << "The parameter values read in from actual.param are:" << std::endl;
	for(int i = 0; i < temp_values.size(); i++){
    std::cout << temp_values[i] << " ";
	}
  std::cout << std::endl << "Thses have given us parameter values of:" << std::endl;
	for(int i = 0; i<datameans.size(); i++){
    std::cout << m_Model->GetScalarOutputNames()[i] << " " << datameans[i] << std::endl;
	}

	return datameans;
}

std::vector<double> madai::LikelihoodDistribution_RHIC::GetRealData(){
  std::vector<double> datameans;

  std::string data_filename = m_Model->m_DirectoryName + "/exp_data/results.dat";
  std::fstream data;
  std::string type, obsv_name;
  int count=0;
  std::vector<std::string> PNames;
  PNames = m_Model->GetScalarOutputNames();

  int numparams = PNames.size();
  std::cout << "There are " << numparams  << " observables used in the emulator." << std::endl;
  m_DataMean=new double[numparams];
  m_DataError=new double[numparams];
  std::vector<double> temp (numparams, .01);
  double dump;

  data.open(data_filename.c_str(), std::fstream::in);
  if(data){
    while(data >> type){
      if(std::strcmp(type.c_str(), "double") == 0){
        data >> obsv_name;
        int index = FindParam(obsv_name, PNames);
        if(index != -1){ //returns -1 if not found
          data >> m_DataMean[index]; //Mean
          data >> m_DataError[index]; //Error
          count++;
          std::cout << obsv_name << " index: " << index << " " << m_DataMean[index] << std::endl;
        }else{
          std::cout << "Not using observable: " << obsv_name << std::endl;
          data >> dump; //we aren't using the observable, so we need to get the data out of the stream
          data >> dump;
          //exit(1);
        }
      }else{
        if(std::strncmp(type.c_str(), "#", 1) == 0){
          std::string temp;
          std::getline(data, temp, '\n');
        }
        else{
          std::cout << "Unrecognized variable type " << type << std::endl;
          exit(1);
        }
      }
    }
    data.close();
  }else{
    std::cout << "Warning: Unable to open data file in model directory." << std::endl;
    exit(1);
  }

  if(count!=numparams){
    std::cout << "Not all emulated observables found! count=" << count << " numparams= " << numparams << std::endl;
    exit(1);
  }

  datameans.assign(m_DataMean,m_DataMean+numparams);

  if(m_Verbose){
    std::cout << "Data array: " << std::endl;
    for(int i = 0; i<numparams; i++){
      cout << m_Model->GetScalarOutputNames()[i] << " " << datameans[i] << endl;
    }
  }
  return datameans;
}

std::vector<double> madai::LikelihoodDistribution_RHIC::GetRealError(){
  std::vector<double> dataerrors;
  std::string data_filename = m_Model->m_DirectoryName + "/exp_data/results.dat";

  std::fstream data;
  std::string type, obsv_name;
  int count = 0;

  std::vector<std::string> PNames;
  PNames = m_Model->GetScalarOutputNames();

  int numparams = PNames.size();
  m_DataMean=new double[numparams];
  m_DataError=new double[numparams];
  std::vector<double> temp (numparams, .01);
  double dump;

  data.open(data_filename.c_str(), std::fstream::in);
  if(data){
    while(data >> type){
      if(std::strcmp(type.c_str(), "double") == 0){
        data >> obsv_name;
        int index = FindParam(obsv_name, PNames);
        if(index != -1){ //returns -1 if not found
          data >> m_DataMean[index]; //Mean
          data >> m_DataError[index]; //Error
          count++;
          std::cout << obsv_name << " index: " << index << " " << m_DataMean[index] << std::endl;
        }else{
          data >> dump; //we aren't using the observable, so we need to get the data out of the stream
          data >> dump;
        }
      }else{
        if(std::strncmp(type.c_str(), "#", 1) == 0){
          std::string temp;
          std::getline(data, temp, '\n');
        }
        else{
          std::cout << "Unrecognized variable type " << type << std::endl;
          exit(1);
        }
      }
    }
    data.close();
  }else{
    std::cout << "Warning: Unable to open data file in model directory." << std::endl;
    exit(1);
  }

  if(count!=numparams){
    std::cout << "Not all emulated observables found!" << std::endl;
    exit(1);
  }

  dataerrors.assign(m_DataError,m_DataError+numparams);

  if(m_Verbose){
    for(int i = 0; i<numparams; i++){
      std::cout << "Error array: " << std::endl;
      std::cout << "i: " << i << " Error: " << dataerrors[i] << std::endl;
    }
  }
  return dataerrors;
}

int madai::LikelihoodDistribution_RHIC::FindParam(std::string name, std::vector<std::string> PNames){ 
	int out = -1;
	int i = 0;
	bool Found = false;
	
	while(i < PNames.size()){
		// cout << "FindParam: Comparing " << name << " to " << PNames[i] << endl;
		if(std::strcmp(PNames[i].c_str(), name.c_str()) == 0){
			if(!Found){
				out = i;
				Found = true;
			}else{ //A matching parameter has already been found, multiple parameters with the same name.
                std::cout << PNames[out] << std::endl;
                std::cout << PNames[i] << std::endl;
                std::cout << "In ProposalName::FindParam; Duplicate parameter names found. Please change parameter names." << std::endl;
				exit(1);
			}
		}
		i++;
	}
	return out;
}
#endif
