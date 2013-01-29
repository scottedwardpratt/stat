#include "RHIC_PCA_Model.h"

madai::RHIC_PCA_Model::RHIC_PCA_Model()
{
  this->stateFlag = UNINITIALIZED;
  this->m_Quad = NULL;
}

madai::RHIC_PCA_Model::RHIC_PCA_Model( const std::string info_dir )
{
  this->stateFlag = UNINITIALIZED;
  this->m_Quad = NULL;
  this->LoadConfiguration(info_dir);
}

madai::RHIC_PCA_Model::~RHIC_PCA_Model()
{
  if( m_Likelihood != NULL )
    delete m_Likelihood;
  if( m_Prior != NULL )
    delete m_Prior;
  if( m_Quad != NULL )
    delete m_Quad;
}

madai::RHIC_PCA_Model::ErrorType
madai::RHIC_PCA_Model::GetScalarOutputs( const std::vector< double > & parameters,
                                         std::vector< double > & scalars ) const
{
  std::vector< double > Means;
  std::vector< double > Errors;
  if(m_UseEmulator){
  m_Quad->QueryQuad( parameters, Means, Errors );
  } else {
  // determine another way to fill the vectors
  }
  if(Means.size() != Errors.size()){
  std::cerr << "Means and Errors from PCA Model arn't the same size" << std::endl;
  return OTHER_ERROR;
  } else {
  for(unsigned int i = 0; i < Means.size(); i++){
  scalars.push_back(Means[i]);
  scalars.push_back(Errors[i]);
  }
  }
  return NO_ERROR;
}

// Not implemented yet.  Should we do this numerically?
/** Get both scalar values and the gradient of the parameters. */
madai::RHIC_PCA_Model::ErrorType
madai::RHIC_PCA_Model::GetScalarAndGradientOutputs(const std::vector< double > & parameters,
                                                   const std::vector< bool > & activeParameters,
                                                   std::vector< double > & scalars,
                                                   unsigned int outputIndex, std::vector< double > & gradient) const
{
  return OTHER_ERROR;
}

// For interaction with the mcmc
madai::RHIC_PCA_Model::ErrorType
madai::RHIC_PCA_Model::GetLikeAndPrior( const std::vector< double > & parameters,
                                        double & Like,
                                        double & Prior) const
{
  std::vector< double > outputs;
  std::vector< double > ModelMeans;
  std::vector< double > ModelErrors;
  this->GetScalarOutputs(parameters, outputs);
  for(unsigned int i = 0 ; i < outputs.size(); i++ ){
  if( i%2 == 0 ){
  ModelMeans.push_back(outputs[i]);
  } else {
  ModelErrors.push_back(outputs[i]);
  }
  }

  Like = m_Likelihood->Evaluate( ModelMeans, ModelErrors );
  Prior = m_Prior->Evaluate( parameters );

  return NO_ERROR;
}

madai::RHIC_PCA_Model::ErrorType
madai::RHIC_PCA_Model::LoadDistributions()
{
  if(m_UseEmulator){
  bool sep_map = parameter::getB(m_ParameterMap, "LIKELIHOOD_PARAMETER_MAP", false);
  if(sep_map){
  std::string pmf = m_DirectoryName + "/defaultpars/likelihood.param";
  parameterMap* ParMap = new parameterMap;
  parameter::ReadParsFromFile(*ParMap, pmf);
  m_Quad = new QuadHandler(ParMap, this);
  } else {
  m_Quad = new QuadHandler(&m_ParameterMap, this);
  }
  } else {
  std::cerr << "RHIC_PCA requires the use of the Quad emulator" << std::endl;
  this->stateFlag = ERROR;
  return OTHER_ERROR;
  }

  m_Likelihood = new LikelihoodDistribution_RHIC_PCA(this);
  m_Prior = new PriorDistribution_RHIC_PCA(this);
}
