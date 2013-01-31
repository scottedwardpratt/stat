#include "CosmoModel.h"

namespace madai {

CosmoModel::CosmoModel()
{
  this->m_StateFlag = UNINITIALIZED;
}

CosmoModel::CosmoModel( const std::string info_dir )
{
  this->m_StateFlag = UNINITIALIZED;
  this->LoadConfiguration(info_dir);
}

CosmoModel::~CosmoModel()
{
  if( m_Likelihood != NULL )
    delete m_Likelihood;
  if( m_Prior != NULL )
    delete m_Prior;
}

CosmoModel::ErrorType
CosmoModel::GetScalarOutputs( const std::vector< double > & parameters,
                              std::vector< double > & scalars ) const
{
  std::stringstream ss;
  std::ifstream inputfile;
  
  ss << "cosmosurvey";
    
  for(int i = 0; i< m_Parameters.size(); i++)
    ss << " -" << m_Parameters[i].m_Name << " " << parameters[i];
  
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
    scalars.push_back(temp);
  }
  return NO_ERROR;
}

// Not implemented yet.  Should we do this numerically?
/** Get both scalar values and the gradient of the parameters. */
CosmoModel::ErrorType 
CosmoModel::GetScalarAndGradientOutputs(const std::vector< double > & parameters,
                                        const std::vector< bool > & activeParameters,
                                        std::vector< double > & scalars,
                                        unsigned int outputIndex, std::vector< double > & gradient) const 
{
  return OTHER_ERROR;
}

// For interaction with the mcmc
CosmoModel::ErrorType
CosmoModel::GetLikeAndPrior( const std::vector< double > & parameters,
                             double & Like,
                             double & Prior) const
{
  std::vector< double > ModelMeans;
  std::vector< double > ModelErrors;
  this->GetScalarOutputs(parameters, ModelMeans);
  
  Like = m_Likelihood->Evaluate( ModelMeans, ModelErrors );
  Prior = m_Prior->Evaluate( parameters );
  
  return NO_ERROR;
}

CosmoModel::ErrorType
CosmoModel::LoadDistributions()
{
  if(m_UseEmulator || m_ProcessPipe){
    std::cerr << "Emulator is not needed for this model. Turn off USE_EMULATOR and PROCESS_PIPE" << std::endl;
    this->m_StateFlag = ERROR;
    return OTHER_ERROR;
  }
  m_Likelihood = new LikelihoodDistribution_Cosmo(this);
  m_Prior = new PriorDistribution_Cosmo(this);
}

} // end namespace madai
