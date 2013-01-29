#include "Distribution.h"
#include <time.h>

madai::LikelihoodDistribution_RHIC_PCA::LikelihoodDistribution_RHIC_PCA(madai::Model *in_Model){
  m_Model = in_Model;
  m_SepMap = parameter::getB(m_Model->m_ParameterMap, "LIKELIHOOD_PARAMETER_MAP", false);

  if(m_SepMap){
  std::string parmapfile = m_Model->m_DirectoryName + "/defaultpars/likelihood.param";
  m_ParameterMap = new parameterMap;
  parameter::ReadParsFromFile(*m_ParameterMap, parmapfile);
  }else{
  m_ParameterMap = &(m_Model->m_ParameterMap);
  }

  m_Timing = parameter::getB(*m_ParameterMap, "TIMING", false) || parameter::getB(*m_ParameterMap, "TIME_LIKELIHOOD", false);
  m_Verbose = parameter::getB(*m_ParameterMap, "VERBOSE", false) || parameter::getB(*m_ParameterMap, "VERBOSE_LIKELIHOOD", false);

  m_Data = GetRealData();
}

madai::LikelihoodDistribution_RHIC_PCA::~LikelihoodDistribution_RHIC_PCA(){
  //delete m_Quad;
}

double madai::LikelihoodDistribution_RHIC_PCA::Evaluate(std::vector<double> ModelMeans,
                                                        std::vector<double> ModelErrors){
  clock_t begintime;
  double likelihood;

  if(m_Timing){
  begintime = clock();
  }

  //Initialize GSL containers
  int N = ModelErrors.size();
  gsl_matrix * sigma = gsl_matrix_calloc(N,N);
  gsl_vector * model = gsl_vector_alloc(N);
  gsl_vector * mu = gsl_vector_alloc(N);


  //Read in appropriate elements
  for(int i = 0; i<N; i++){
  gsl_matrix_set(sigma,i,i,ModelErrors[i]);
  gsl_vector_set(model, i,ModelMeans[i]);
  gsl_vector_set(mu, i, m_Data[i]);
  }

  likelihood = Log_MVNormal(*model, *mu, *sigma);

  if(!(m_Model->m_LogLike)){ //If you don't want the loglikelihood, and we've used Log_MVN, we have to exponentiate.
  likelihood = exp(likelihood);
  }

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

  return likelihood;
}

std::vector<double> madai::LikelihoodDistribution_RHIC_PCA::GetRealData(){
  std::vector<double> datameans;
  float myints[] = {403.723,467.153,743.79,1042.58,5.28,4.81,5.47,180.682,460.801,739.43,1025.55,0.0891618,0.0271948,0.0498033,4.27,3.99,4.53};
  datameans.assign (myints,myints+17);
  //float myints[] = {0.658581,0.545024,1.4456,-2.03515,1.22693,0.731763};
  //datameans.assign (myints,myints+6);
  return datameans;
}
