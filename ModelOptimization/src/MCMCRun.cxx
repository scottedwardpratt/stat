#include "MCMCRun.h"

madai::MCMCRun::MCMCRun(const madai::Model *in_model) :
  Optimizer( in_model ) {
	m_LocalParameterMap = m_Model->m_ParameterMap;
  
	m_BurnIn = parameter::getI(m_LocalParameterMap, "BURN_IN", 0);
	m_RandomTheta0 = parameter::getB(m_LocalParameterMap, "RANDOM_THETA0", false);
	m_VizTrace = parameter::getB(m_LocalParameterMap, "VISUALIZE_TRACE", true);
	m_Quiet = parameter::getB(m_LocalParameterMap, "QUIET", false);
  m_AppendTrace = parameter::getB(m_LocalParameterMap, "APPEND_TRACE", false);
  m_RescaledTrace = parameter::getB(m_LocalParameterMap, "RESCALED_TRACE", false);
	m_LogPrior = parameter::getB(m_LocalParameterMap, "LOGPRIOR", true);
	m_LogProposal = parameter::getB(m_LocalParameterMap, "LOGPROPOSAL", true);
	m_CreateTrace = parameter::getB(m_LocalParameterMap, "CREATE_TRACE", true);
    
  m_RandomNumber = new CRandom(1234);
    
  if(m_RandomTheta0){
    m_InitialTheta = this->GetRandomTheta0(time(NULL));
  }
  else{
    m_InitialTheta = this->GetTheta0FromFile();
  }
  m_CurrentParameters = m_InitialTheta;
}

madai::MCMCRun::~MCMCRun(){
    
}

void madai::MCMCRun::NextIteration(madai::Trace *ThetaOutsList){
  double LOGBF, alpha;
  if(m_Model->m_LogLike){
    LOGBF = 0;
  }else{
    LOGBF = 1;
  }
  m_ScaleNew = (rand() / double(RAND_MAX));
  
  std::vector<double> Temp_Theta;
  Temp_Theta = m_Model->m_Proposal->Iterate(m_CurrentParameters,m_ScaleNew,m_ActiveParameters);
  m_LikelihoodNew = m_Model->m_Likelihood->Evaluate(Temp_Theta);
  if(m_IterationNumber == 1){
    m_BestLikelihood=m_LikelihoodNew;
    m_BestParameterSet = m_InitialTheta;
    if(m_CreateTrace){
      m_Visualizer = new VizHandler(this);
      m_VizCount = parameter::getI(m_LocalParameterMap, "VIZ_COUNT", floor(ThetaOutsList->m_MaxIterations/200));
    }
  }
  m_PriorNew = m_Model->m_Prior->Evaluate(Temp_Theta);
  m_ProposalNew = m_Model->m_Proposal->Evaluate(m_CurrentParameters,Temp_Theta,m_ScaleCurrent);
  m_ProposalCurrent = m_Model->m_Proposal->Evaluate(Temp_Theta,m_CurrentParameters,m_ScaleNew);
  
  // Likelihood
  if(m_Model->m_LogLike){
    if(!m_Quiet){
      printf(" ll_new=%g, ll_current=%g\n",m_LikelihoodNew,m_LikelihoodCurrent);
    }
    LOGBF += m_LikelihoodNew-m_LikelihoodCurrent;
  } else {
    if(!m_Quiet){
      printf(" l_new=%g, l_current=%g\n",exp(m_LikelihoodNew),exp(m_LikelihoodCurrent));
    }
    LOGBF *= exp(m_LikelihoodNew)/exp(m_LikelihoodCurrent);
  }
  
  // Prior
  if(m_LogPrior){
    LOGBF += (log(m_PriorNew)-log(m_PriorCurrent));
  } else {
    LOGBF *= (m_PriorNew/m_PriorCurrent);
  }
  
  if(!m_Quiet){
    printf(" Prior_New=%g, Prior_Current=%g\n",m_PriorNew,m_PriorCurrent);
  }
  
  // Proposal		
  if(m_Model->m_LogLike){
    LOGBF += (log(m_ProposalCurrent)-log(m_ProposalNew));
  } else {
    LOGBF *= m_ProposalCurrent/m_ProposalNew;
  }
  
  if(!m_Quiet){
    printf(" Proposal_New=%g, Proposal_Current=%g\n",m_ProposalNew,m_ProposalCurrent);
  }
  
  if(m_Model->m_LogLike){
    alpha = min(1.0,exp(LOGBF));
  } else {
    alpha = min(1.0,LOGBF);
  }
  
  if(!m_Quiet){
    printf("%5d\talpha=%6.5f\t",m_IterationNumber,alpha);
  }
  
  if(alpha > (this->m_RandomNumber->ran())) { //Accept the proposed set.
    if(!m_Quiet){
      printf("Accept\n");
    }
    if(m_IterationNumber > m_BurnIn){
      m_AcceptCount++;
    }
    m_LikelihoodCurrent = m_LikelihoodNew;
    m_PriorCurrent = m_PriorNew;
    m_CurrentParameters = Temp_Theta;
    m_ScaleCurrent = m_ScaleNew;
    if(m_LikelihoodCurrent > m_BestLikelihood && m_IterationNumber > 1){
      m_BestLikelihood=m_LikelihoodNew;
      m_BestParameterSet=m_CurrentParameters;
      if(!m_Quiet){
        if(m_Model->m_LogLike){
          printf("XXXXXXXXX YIPPEE!! Best parameters so far, loglikelihood=%g\n",m_BestLikelihood);
        }
        else{
          printf("XXXXXXXXX YIPPEE!! Best parameters so far, likelihood=%g\n",m_BestLikelihood);
        }
      }
    }
  }else{
    if(!m_Quiet){
      printf("Reject\n");
    }
  }
  
  if(m_IterationNumber > m_BurnIn){ // We are just tossing everything in the burn in period.
    if(m_RescaledTrace){
      std::vector<double> RescaledParameters;
      double* range = new double[2]();
      for( int k = 0; k < m_CurrentParameters.size(); k++){
        m_Model->GetRange(k, range);
        RescaledParameters.push_back((m_CurrentParameters[k]-range[0])/(range[1]-range[0]));
      }
      ThetaOutsList->add(RescaledParameters);
    }else{
      ThetaOutsList->add(m_CurrentParameters);
    }
  }
  
  double range[2];
  for(int k = 0; k < m_Model->GetNumberOfParameters(); k++){ //These ParamValus are used for the density plots
    m_Model->GetRange(k,range);
    m_ParameterValues[k] = (m_CurrentParameters[k] - range[0])/(range[1]-range[0]);
  }
  
  if((m_IterationNumber > m_BurnIn) && ((m_IterationNumber+1) % (ThetaOutsList->m_Writeout) == 0)){
    std::cout << "Writing out." << std::endl;
    if(m_CreateTrace && (m_IterationNumber!=1)){
      m_Visualizer->UpdateTraceFig(ThetaOutsList);
    }
    ThetaOutsList->WriteOut(m_Model->GetParameters());
  }
}

std::vector<double> 
madai::MCMCRun::GetRandomTheta0(int seed){ //This creates random theta 0s
	srand(seed);
  std::cout << "We are using random theta0 values. They are:" << std::endl;
	double *range = new double[2];
  std::vector<double> temp_values(this->m_Model->GetNumberOfParameters(), 0.0);
    
  for(unsigned int i = 0; i<this->m_Model->GetNumberOfParameters();i++){
    if(m_Model->m_PrescaledParams){
      range[0]=0;
      range[1]=1;
    }else{
      this->m_Model->GetRange(i,range);
    }
    temp_values[i] = double(rand() % int((range[1] - range[0])*1000))/1000+range[0];
  }
    
  return temp_values;
}

std::vector<double> 
madai::MCMCRun::GetTheta0FromFile(){
  parameterMap parmap;
  std::string theta0_filename = this->m_Model->m_ParameterFile + "/theta0.param";
  parameter::ReadParsFromFile(parmap, theta0_filename);
  std::vector<std::string> temp_names = parameter::getVS(parmap,"NAMES","");
  std::vector<double> temp_values = parameter::getV(parmap, "VALUES","");
  std::vector<double> read_theta;
  unsigned int j, i=0;
    
  if(!((m_Model->GetParameters()).empty())){
    double *range = new double[2];
    std::vector<madai::Parameter>::const_iterator itr = (m_Model->GetParameters()).begin();
    std::vector<std::string>::const_iterator titr;
    for( itr; itr < (m_Model->GetParameters()).end(); itr++){
      j=0;
      for( titr = temp_names.begin(); titr < temp_names.end(); titr++){
        if((itr->m_Name) == *titr){
          read_theta.push_back(temp_values[j]);
          break;
        }else{
          j++;
        }
      }
      m_Model->GetRange(i,range);
      if(read_theta[i] < range[0]){
        std::cout << "Parameter below min value.\n\n";
        exit(1);
      } else if(read_theta[i] > range[1]){
        std::cout << "Parameter above max value.\n\n";
        exit(1);
      }
      i++;
    }
  }
  return read_theta;
}

