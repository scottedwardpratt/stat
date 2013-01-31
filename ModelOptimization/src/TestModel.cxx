#include "TestModel.h"

madai::TestModel::TestModel()
{
  this->stateFlag = UNINITIALIZED;
}

madai::TestModel::TestModel( const std::string info_dir )
{
  this->stateFlag = UNINITIALIZED;
  this->LoadConfiguration(info_dir);
}

madai::TestModel::~TestModel()
{
  if( m_Likelihood != NULL )
    delete m_Likelihood;
  if( m_Prior != NULL )
    delete m_Prior;
}

madai::TestModel::ErrorType
madai::TestModel::GetScalarOutputs( const std::vector< double > & parameters,
                                    std::vector< double > & scalars ) const
{
  std::vector< double > Means;
  std::vector< double > Errors;
  if(m_UseEmulator){
  std::cerr << "TestModel should not use an emulator. Set the UseEmulator flag to false" << std::endl;
  return OTHER_ERROR;
  } else {
  double errors[]={0,0,0,0};
  Errors.assign (errors,errors+4);
  double means[4];
  for(int i = 0; i < 4; i++){
  means[i] = parameters[i];
  }
  Means.assign(means,means+4);
  }
  for( int j = 0; j < 4; j++ ) {
  scalars.push_back(Means[j]);
  scalars.push_back(Errors[j]);
  }
  return NO_ERROR;
}

// Not implemented yet.  Should we do this numerically?
/** Get both scalar values and the gradient of the parameters. */
madai::TestModel::ErrorType
madai::TestModel::GetScalarAndGradientOutputs(const std::vector< double > & parameters,
                                              const std::vector< bool > & activeParameters,
                                              std::vector< double > & scalars,
                                              unsigned int outputIndex, std::vector< double > & gradient) const
{
  return OTHER_ERROR;
}

madai::TestModel::ErrorType
madai::TestModel::GetLikeAndPrior( const std::vector< double > & parameters,
                                   double & Like,
                                   double & Prior) const
{
  std::vector< double > outputs;
  std::vector< double > ModelMeans;
  std::vector< double > ModelErrors;
  this->GetScalarOutputs(parameters, outputs);
  for(unsigned int i = 0; i < outputs.size(); i++){
  if( i%2 == 0 ){
  ModelMeans.push_back( outputs[i] );
  } else {
  ModelErrors.push_back( outputs[i] );
  }
  }

  Like = m_Likelihood->Evaluate( ModelMeans, ModelErrors );
  Prior = m_Prior->Evaluate( parameters );

  return NO_ERROR;
}

madai::TestModel::ErrorType
madai::TestModel::LoadDistributions()
{
  if( m_UseEmulator || m_ProcessPipe ){
  std::cerr << "Emulator is not used for the test model. Turn off USE_EMULATOR and PROCESS_PIPE" << std::endl;
  this->stateFlag = ERROR;
  return OTHER_ERROR;
  }
  m_Likelihood = new LikelihoodDistribution_Test(this);
  m_Prior = new PriorDistribution_Test(this);
}
