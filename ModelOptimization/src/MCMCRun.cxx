#include "MCMCRun.h"

madai::MCMCRun::MCMCRun(const madai::Model *in_model) :
Optimizer( in_model ) {
	m_LocalParameterMap = m_Model->m_ParameterMap;
	m_TraceDirectory = m_Model->m_DirectoryName + "/trace/" + m_Model->m_ConfigurationName;
	
	m_MaxIterations = parameter::getI(m_LocalParameterMap, "MAX_ITERATIONS", 500);
	m_Writeout = parameter::getI(m_LocalParameterMap, "WRITEOUT", 100);
	m_BurnIn = parameter::getI(m_LocalParameterMap, "BURN_IN", 0);
	m_RandomTheta0 = parameter::getB(m_LocalParameterMap, "RANDOM_THETA0", false);
	m_VizTrace = parameter::getB(m_LocalParameterMap, "VISUALIZE_TRACE", true);
	m_Quiet = parameter::getB(m_LocalParameterMap, "QUIET", false);
  m_RescaledTrace = parameter::getB(m_LocalParameterMap, "RESCALED_TRACE", false);
	m_AppendTrace = parameter::getB(m_LocalParameterMap, "APPEND_TRACE", false);
  //m_LogLike = parameter::getB(m_LocalParameterMap, "LOGLIKE", true);
	m_LogPrior = parameter::getB(m_LocalParameterMap, "LOGPRIOR", true);
	m_LogProposal = parameter::getB(m_LocalParameterMap, "LOGPROPOSAL", true);
  //SUPPRESS_ERRORS = parameter::getB(m_LocalParameterMap, "SUPPRESS_ERRORS", false);
	m_CreateTrace = parameter::getB(m_LocalParameterMap, "CREATE_TRACE", true);
  //MODEL = parameter::getS(m_LocalParameterMap,"MODEL","NOMODEL");
    
  m_RandomNumber = new CRandom(1234);
    
  if(m_RandomTheta0){
    m_InitialTheta = this->GetRandomTheta0(time(NULL));
  }
  else{
    m_InitialTheta = this->GetTheta0FromFile();
  }
  m_ThetaOutsList = new madai::Trace(this);
  m_ThetaOutsList->add(m_InitialTheta);
  //CurrentEl = madai::TraceElement(ThetaOutsList);
    
  //Likelihood_Current=0;
    
	if(m_CreateTrace){
		m_Visualizer = new VizHandler(this);
		m_VizCount = parameter::getI(m_LocalParameterMap, "VIZ_COUNT", floor(m_MaxIterations/200));
		//Visualizer->UpdateTraceFig();
	}
    
	if(m_AppendTrace){
    std::string addon = "";
		bool Done = false;
		int filecount = 0;
		while(!Done){
			struct stat st;
      std::stringstream ss;
      std::string tempfile = m_TraceDirectory + addon;
			if(stat(tempfile.c_str(), &st)==0){
				//directory exists.
        std::cout << tempfile << " exists, trying next option..." << std::endl;
				filecount++;
				ss << "_" << filecount;
				addon = ss.str();
				ss.str(string());
			}else{
				//doesn't exist
				Done = true;
				m_TraceDirectory = tempfile;
			}
		}
	}else{
    std::cout << "Deleting prior trace data." << std::endl;
    std::string cmd = "rm " + m_TraceDirectory + "/output*.dat " + m_TraceDirectory + "/trace.dat";
    std::system(cmd.c_str());
	}
	
  std::string command = "mkdir -p "+ m_TraceDirectory;
	
  std::system(command.c_str());
	/*if(!m_Quiet){
     printf("Iteration\tAlpha\tResult\n");
     }*/
}


madai::MCMCRun::MCMCRun(const madai::Model *in_model, const std::vector<double> Theta0) :
Optimizer( in_model ) {
  m_LocalParameterMap = m_Model->m_ParameterMap;
	m_TraceDirectory = m_Model->m_DirectoryName + "/trace/" + m_Model->m_ConfigurationName;
	
	m_MaxIterations = parameter::getI(m_LocalParameterMap, "MAX_ITERATIONS", 500);
	m_Writeout = parameter::getI(m_LocalParameterMap, "WRITEOUT", 100);
	m_BurnIn = parameter::getI(m_LocalParameterMap, "BURN_IN", 0);
	m_RandomTheta0 = parameter::getB(m_LocalParameterMap, "RANDOM_THETA0", false);
	m_VizTrace = parameter::getB(m_LocalParameterMap, "VISUALIZE_TRACE", true);
	m_Quiet = parameter::getB(m_LocalParameterMap, "QUIET", false);
  m_RescaledTrace = parameter::getB(m_LocalParameterMap, "RESCALED_TRACE", false);
	m_AppendTrace = parameter::getB(m_LocalParameterMap, "APPEND_TRACE", false);
  //m_LogLike = parameter::getB(m_LocalParameterMap, "LOGLIKE", true);
	m_LogPrior = parameter::getB(m_LocalParameterMap, "LOGPRIOR", true);
	m_LogProposal = parameter::getB(m_LocalParameterMap, "LOGPROPOSAL", true);
  //SUPPRESS_ERRORS = parameter::getB(m_LocalParameterMap, "SUPPRESS_ERRORS", false);
	m_CreateTrace = parameter::getB(m_LocalParameterMap, "CREATE_TRACE", true);
  //MODEL = parameter::getS(m_LocalParameterMap,"MODEL","NOMODEL");
    
  m_RandomNumber = new CRandom(1234);
  /*if(std::strcmp(MODEL.c_str(),"RHIC")==0){
    Likelihood = new LikelihoodDistribution_RHIC(this);
    Prior = new PriorDistribution_RHIC(this);
  }
	else{
		printf("Must define parameter MODEL in parameter file, or yours is unrecognized\n");
    std::cout << "Note that the only model supported at the moment is RHIC" << std::endl;
    exit(1);
	} 
  Proposal = new ProposalDistribution(this);*/
    
	if(m_RandomTheta0){
    m_InitialTheta = this->GetRandomTheta0(time(NULL));
  }
	else{
    m_InitialTheta = Theta0;
	}
  m_ThetaOutsList = new madai::Trace(this);
  m_ThetaOutsList->add(m_InitialTheta);
  //CurrentEl = madai::TraceElement(ThetaOutsList);
    
	//Likelihood_Current=0;
    
	if(m_CreateTrace){
		m_Visualizer = new VizHandler(this);
		m_VizCount = parameter::getI(m_LocalParameterMap, "VIZ_COUNT", floor(m_MaxIterations/200));
		//Visualizer->UpdateTraceFig();
	}
	
	if(m_AppendTrace){
    std::string addon = "";
		bool Done = false;
		int filecount = 0;
		while(!Done){
			struct stat st;
        std::stringstream ss;
        std::string tempfile = m_TraceDirectory + addon;
			if(stat(tempfile.c_str(), &st)==0){
				//directory exists.
        std::cout << tempfile << " exists, trying next option..." << std::endl;
				filecount++;
				ss << "_" << filecount;
				addon = ss.str();
				ss.str(string());
			}else{
				//doesn't exist
				Done = true;
				m_TraceDirectory = tempfile;
			}
		}
	}else{
    std::cout << "Deleting prior trace data." << std::endl;
    std::string cmd = "rm " + m_TraceDirectory + "/output*.dat " + m_TraceDirectory + "/trace.dat";
    std::system(cmd.c_str());
	}
	
  std::string command = "mkdir -p "+ m_TraceDirectory;
	
  std::system(command.c_str());
	/*if(!m_Quiet){
     printf("Iteration\tAlpha\tResult\n");
     }*/
}

madai::MCMCRun::~MCMCRun(){
    
}

// Not implemented yet. May change over NextIteration() funct to this at some point
// don't know if structure is conducive to the change though.
void madai::MCMCRun::NextIteration(madai::Trace *ThetaOutsList){
  std::cerr << "MCMCRun::NextIteration(Trace*) not defined" << std::endl;
  exit(1);
}

void madai::MCMCRun::NextIteration(){
  double LOGBF, alpha;
  if(m_Model->m_LogLike){
    LOGBF = 0;
  }else{
    LOGBF = 1;
  }
  m_Scale_New = (rand() / double(RAND_MAX));
    
  std::vector<double> Temp_Theta;
  Temp_Theta = m_Model->m_Proposal->Iterate(m_CurrentParameters,m_Scale_New,m_ActiveParameters);
    
  m_Likelihood_New = m_Model->m_Likelihood->Evaluate(Temp_Theta);
  if(m_IterationNumber==1){
    m_BestLikelihood=m_Likelihood_New;
    m_BestParameterSet = m_ThetaZeroPtr->m_ParameterValues;
    //BestParameterSetPtr=ThetaZeroPtr;
  }
  m_Prior_New = m_Model->m_Prior->Evaluate(Temp_Theta);
  m_Proposal_New = m_Model->m_Proposal->Evaluate(m_CurrentParameters,Temp_Theta,m_Scale_Current);
  m_Proposal_Current = m_Model->m_Proposal->Evaluate(Temp_Theta,m_CurrentParameters,m_Scale_New);
    
    //std::cerr << "Likelihood of proposed set: " << m_Likelihood_New << std::endl;
    //std::cerr << "Likelihood of current set: " << m_Likelihood_Current << std::endl;
    //std::cerr << "Prior of proposed set: " << m_Prior_New << std::endl;
    //std::cerr << "Prior of current set: " << m_Prior_Current << std::endl;
    //std::cerr << "Proposal of proposed set: " << m_Proposal_New << std::endl;
    //std::cerr << " Proposal of current set: " << m_Proposal_Current << std::endl;
    
    
    // Likelihood
  if(m_Model->m_LogLike){
    if(!m_Quiet){
      printf(" ll_new=%g, ll_current=%g\n",m_Likelihood_New,m_Likelihood_Current);
    }
    LOGBF += m_Likelihood_New-m_Likelihood_Current;
  } else {
    if(!m_Quiet){
      printf(" l_new=%g, l_current=%g\n",exp(m_Likelihood_New),exp(m_Likelihood_Current));//200
    }
    LOGBF *= exp(m_Likelihood_New)/exp(m_Likelihood_Current);
  }
    
  // Prior
  if(m_LogPrior){
    LOGBF += (log(m_Prior_New)-log(m_Prior_Current));
  } else {
    LOGBF *= (m_Prior_New/m_Prior_Current);
  }
    
  if(!m_Quiet){
    printf(" Prior_New=%g, Prior_Current=%g\n",m_Prior_New,m_Prior_Current);
  }
    
  // Proposal		
  if(m_Model->m_LogLike){
    LOGBF += (log(m_Proposal_Current)-log(m_Proposal_New));
  } else {
    LOGBF *= m_Proposal_Current/m_Proposal_New;
  }
    
  if(!m_Quiet){
    printf(" Proposal_New=%g, Proposal_Current=%g\n",m_Proposal_New,m_Proposal_Current);
  }
    
  if(m_Model->m_LogLike){
    alpha = min(1.0,exp(LOGBF));
  } else {
    alpha = min(1.0,LOGBF);
  }
    
  //std::cout << LOGBF << std::endl;
    
  if(!m_Quiet){
    printf("%5d\talpha=%6.5f\t",m_IterationNumber,alpha);
    //printf("LOGBF=%6.5f\t",alpha);
  }
  if(alpha > (this->m_RandomNumber->ran())) { //Accept the proposed set.
    if(!m_Quiet){
      printf("Accept\n");
    }
    if(m_IterationNumber > m_BurnIn){
      m_AcceptCount++;
    }
    m_Likelihood_Current = m_Likelihood_New;
    m_Prior_Current = m_Prior_New;
    //m_Proposal_Current = m_Proposal_New;
    m_CurrentParameters = Temp_Theta;
    m_Scale_Current = m_Scale_New;
    if(m_Likelihood_Current>m_BestLikelihood && m_IterationNumber>1){
      m_BestLikelihood=m_Likelihood_New;
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
    m_ThetaOutsList->add(m_CurrentParameters);
  }
    
  double range[2];
  for(int k = 0; k < m_Model->GetNumberOfParameters(); k++){ //These ParamValus are used for the density plots
    m_Model->GetRange(k,range);
    //cout << ThetaList->ParamNames[k] << " " << CurrentParameters.Values[k] << endl;
    //cout << "(" << CurrentParameters.Values[k] << " - " << mcmcconfig->Min_Ranges[k] << ") / (" << mcmcconfig->Max_Ranges[k] << " - " << mcmcconfig->Min_Ranges[k] << ")" << endl; cout.flush();
    //cout << mcmcconfig->Max_Ranges[k] << " " << mcmcconfig->Min_Ranges[k] << endl;
    m_ParameterValues[k] = (m_CurrentParameters[k] - range[0])/(range[1]-range[0]);
    //cout << ParamValues[k] << endl;
  }
    
    /*if((niter > BURN_IN) && CREATE_TRACE){
        if((niter+1) % 5 == 0){
            Visualizer->UpdateTraceFig();
        }
    }*/
    
  if((m_IterationNumber > m_BurnIn) && ((m_IterationNumber+1) % m_Writeout == 0)){
    std::cout << "Writing out." << std::endl;
    if(m_CreateTrace &&(m_IterationNumber!=1)){
      m_Visualizer->UpdateTraceFig();
    }
    m_ThetaOutsList->WriteOut();
    //BestParameterSetPtr->Print();
  }
}

std::vector<double> madai::MCMCRun::GetRandomTheta0(int seed){ //This creates random theta 0s
	srand(seed);
  std::cout << "We are using random theta0 values. They are:" << std::endl;
	double *range = new double[2];
  std::vector<double> temp_values(this->m_Model->GetNumberOfParameters(), 0.0);
    
  for(unsigned int i = 0; i<this->m_Model->GetNumberOfParameters();i++){
    //This might be redundant
    /*if(m_Model->PRESCALED_PARAMS){
      range[0]=0;
      range[1]=1;
    }else{
      this->m_Model->GetRange(i,range);
    }*/
    this->m_Model->GetRange(i,range);
    temp_values[i] = double(rand() % int((range[1] - range[0])*1000))/1000+range[0];
  }
    
  return temp_values;
}

std::vector<double> madai::MCMCRun::GetTheta0FromFile(){
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
    for(itr;itr<(m_Model->GetParameters()).end();itr++){
      j=0;
      for(titr=temp_names.begin();titr<temp_names.end();titr++){
        if((itr->m_Name)==*titr){
          read_theta.push_back(temp_values[j]);
          break;
        }else{
          j++;
        }
      }
      m_Model->GetRange(i,range);
      if(read_theta[i]<range[0]){
        std::cout << "Parameter below min value.\n\n";
        exit(1);
      } else if(read_theta[i]>range[1]){
        std::cout << "Parameter above max value.\n\n";
        exit(1);
      }
      i++;
    }
  }
  return read_theta;
}

double madai::MCMCRun::RescaleTheta(int itern, int parn){
  double range[2];
  m_Model->GetRange(parn,range);
  return (((*m_ThetaOutsList)[itern].m_ParameterValues[parn]-range[0])/(range[1]-range[0]));
}

