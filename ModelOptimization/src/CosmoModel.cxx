#include "CosmoModel.h"

madai::CosmoModel::CosmoModel()
{
  this->stateFlag = UNINITIALIZED;
}

madai::CosmoModel::CosmoModel( const std::string info_dir )
{
  this->stateFlag = UNINITIALIZED;
  this->LoadConfiguration(info_dir);
}

madai::CosmoModel::~CosmoModel()
{
  if( m_Likelihood != NULL )
    delete m_Likelihood;
  if( m_Prior != NULL )
    delete m_Prior;
}

madai::CosmoModel::ErrorType
madai::CosmoModel::GetScalarOutputs( const std::vector< double > & parameters,
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
madai::CosmoModel::ErrorType 
madai::CosmoModel::GetScalarAndGradientOutputs(const std::vector< double > & parameters,
                                              const std::vector< bool > & activeParameters,
                                              std::vector< double > & scalars,
                                              unsigned int outputIndex, std::vector< double > & gradient) const 
{
	return OTHER_ERROR;
}

// For interaction with the mcmc
madai::CosmoModel::ErrorType
madai::CosmoModel::GetLikeAndPrior( const std::vector< double > & parameters,
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

madai::CosmoModel::ErrorType
madai::CosmoModel::LoadDistributions()
{
  if(m_UseEmulator || m_ProcessPipe){
    std::cerr << "Emulator is not needed for this model. Turn off USE_EMULATOR and PROCESS_PIPE" << std::endl;
    this->stateFlag = ERROR;
    return OTHER_ERROR;
  }
  m_Likelihood = new LikelihoodDistribution_Cosmo(this);
  m_Prior = new PriorDistribution_Cosmo(this);
}
