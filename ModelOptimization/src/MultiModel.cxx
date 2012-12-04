#include "MultiModel.h"

madai::MultiModel::MultiModel(){
  this->stateFlag = UNINITIALIZED;
}

madai::MultiModel::MultiModel(const std::string info_dir){
  this->stateFlag = UNINITIALIZED;
  this->LoadConfiguration(info_dir);
}

madai::MultiModel::~MultiModel()
{
  delete m_Likelihood;
  delete m_Proposal;
  delete m_Prior;
}

void discard_line(std::FILE * fp) {
	static int buffersize = 1024;
	char buffer[buffersize];
	std::fgets(buffer, buffersize, fp);
}


bool discard_comments(std::FILE * fp, char comment_character) {
	int c = std::getc(fp);
	if ((c == EOF) || std::ferror(fp)) {
		std::cerr << "premature end of file:(\n";
		return false;
	}
	while (c == comment_character) {
		discard_line(fp);
		c = std::getc(fp);
	}
	if (EOF == std::ungetc(c, fp)) {
		std::cerr << "ungetc error :(\n";
		return false;
	}
	return true;	
}

#include <cctype>

void eat_whitespace(std::istream & i) {
  while (true) {
    if (! std::isspace(i.peek()))
      return;
    if (i.get() == '\n')
      return;
  }
}
void eat_whitespace(std::FILE * fp) {
  while (true) {
  int c = std::fgetc(fp);
  if (! std::isspace(c)) {
    std::ungetc(c,fp);
    return;
  }
  if (c == '\n')
    return;
  }
}

bool discard_comments(std::istream & i, char comment_character) {
	int c = i.peek();
	while (i.good() && ((c == comment_character) || (c == '\n'))) {
		std::string s;
		std::getline(i, s);
		c = i.peek();
	}
}

madai::MultiModel::ErrorType
madai::MultiModel::LoadConfiguration(std::string info_dir){
	m_DirectoryName = info_dir;
	m_ParameterFile = info_dir+"/defaultpars/";
  std::cout << "There should be something here ->" << m_ParameterFile << "<-" << std::endl;
  std::cout << "In config: " << m_ParameterFile << std::endl;
	m_ParameterFileName = m_ParameterFile + "/mcmc.param";
  std::cout << "Reading in " << m_ParameterFileName << std::endl;
    
	parameter::ReadParsFromFile(m_ParameterMap, m_ParameterFileName.c_str());
	m_LogLike = parameter::getB(m_ParameterMap, "LOGLIKE", true);
  m_ModelType = parameter::getS(m_ParameterMap,"MODEL","NOMODEL");
   
    //============================================
	// Reading the parameters out of ranges.dat
  m_PrescaledParams = parameter::getB(m_ParameterMap, "PRESCALED_PARAMS", false);
  std::string filename = info_dir + "/ranges.dat";
  this->LoadConfigurationFile(filename);
  std::cout << "Ranges loaded" << std::endl;
    
  if(!(this->m_Parameters.empty())){
    std::vector<madai::Parameter>::const_iterator itr = this->m_Parameters.begin();
    for( itr++; itr < this->m_Parameters.end(); itr++)
      std::cout << itr->m_Name << " ";
  }else{
    std::cout << "Parameters were not read in!" << std::endl;
    this->stateFlag = ERROR;
    return OTHER_ERROR;
  }
  std::cout << std::endl;
    
  std::string observables_filename = info_dir + "/pcanames.dat";
  this->LoadConfigurationFile(observables_filename);
    
  std::cout << "Emulated Observables Are:" << std::endl;
  if(!(this->m_ScalarOutputNames.empty())){
    std::vector<std::string>::const_iterator itr = this->m_ScalarOutputNames.begin();
    for( itr; itr < this->m_ScalarOutputNames.end(); itr++)
      std::cout << *itr << " ";
    std::cout << std::endl;
  }else{
    std::cout << "Observables were not read in!" << std::endl;
    this->stateFlag=ERROR;
    return OTHER_ERROR;
  }
    
  std::vector<std::string> temp_logparam = parameter::getVS(m_ParameterMap, "LOG_PARAMETERS", "");
	
	for(int i = 0; i < temp_logparam.size(); i++){
		if( std::strcmp(temp_logparam[i].c_str(), "true") == 0 || std::strcmp(temp_logparam[i].c_str(), "True") == 0){
			m_LogParam.push_back(true);
		}else if(std::strcmp(temp_logparam[i].c_str(), "false") == 0 || std::strcmp(temp_logparam[i].c_str(), "False") == 0){
			m_LogParam.push_back(false);
		}else{
      std::cout << "Unrecognized LogParam value " << temp_logparam[i] << std::endl;
			this->stateFlag=ERROR;
      return OTHER_ERROR;
		}
	}
	
  if(std::strcmp(m_ModelType.c_str(),"RHIC")==0){
    m_Likelihood = new LikelihoodDistribution_RHIC(this);
    m_Prior = new PriorDistribution_RHIC(this);
  }else if(std::strcmp(m_ModelType.c_str(), "Cosmosurvey")==0){
    m_Likelihood = new LikelihoodDistribution_Cosmo(this);
    m_Prior = new PriorDistribution_Cosmo(this);
  }else if(std::strcmp(m_ModelType.c_str(), "RHIC_PCA")==0){
    m_Likelihood = new LikelihoodDistribution_RHIC_PCA(this); 
    m_Prior = new PriorDistribution_RHIC_PCA(this);
  }else if(std::strcmp(m_ModelType.c_str(), "TEST")==0){
    m_Likelihood = new LikelihoodDistribution_Test(this);
    m_Prior = new PriorDistribution_Test(this);
  }else{
    printf("Must define parameter MODEL in parameter file, or yours is unrecognized\n");
    this->stateFlag = ERROR;
    return OTHER_ERROR;
  } 
  m_Proposal = new ProposalDistribution(this);

  this->stateFlag=READY;
  return NO_ERROR;
}

/** 
 * Loads a configuration from a file. Reads pcanames.dat or ranges.dat based on fileName
 */
madai::MultiModel::ErrorType
madai::MultiModel::LoadConfigurationFile( const std::string fileName )
{
  std::fstream config_file;
  config_file.open(fileName.c_str(),std::fstream::in);
  if(config_file){
    int num_inputs;
    config_file >> num_inputs;
    if(num_inputs < 1){
      std::cerr << "Number of inputs is < 1" << std::endl;
      this->stateFlag = ERROR;
      return OTHER_ERROR;
    }
    bool fs, ro;
    double MinR, MaxR;
    std::string name, type;
        
    // Check for type of file
    if(fileName == (m_DirectoryName + "/pcanames.dat")){
      fs = false; //File switch to determine how the data is read in
      ro = true; //Tells the program whether to keep reading a line or not
      this->number_of_outputs = num_inputs;
      this->m_ScalarOutputNames.reserve(num_inputs);
    }else if(fileName == (m_DirectoryName+ "/ranges.dat")){
      fs = true; ro = false;
      this->number_of_parameters = num_inputs;
      this->m_Parameters.reserve(num_inputs);
    }else{
      std::cerr << "That fileName doesn't correspond to an mcmc input file" << std::endl;
      this->stateFlag = ERROR;
      return OTHER_ERROR;
    }
    eat_whitespace(config_file);
    int index = 0;
    while(!config_file.eof() && index<num_inputs){
      discard_comments(config_file,'#');
      if(fs){//check if reading ranges
        config_file >> type;
        if(std::strcmp(type.c_str(), "double") == 0){
          ro = true;
        }else{
          ro = false;
        }
      }
      if(ro){//check to see if keep reading
        config_file >> name;
        if(index!=-1 && fs){//check for reading ranges
          if(!m_PrescaledParams){
            config_file >> MinR >> MaxR;
            if(MinR>MaxR){//Flip them
              double temp2 = MinR;
              MinR = MaxR;
              MaxR = temp2;
            }
          }else{
            MinR=0;
            MaxR=1;
          }
        }
        if(fs){
          if(index==0){
            this->m_Parameters.push_back(Parameter(name,MinR,MaxR));
          }else{
            if(name.compare(this->m_Parameters.back().m_Name)!=0)
              this->m_Parameters.push_back(Parameter(name,MinR,MaxR));
          }
        }else{
          if(index==0){
            this->m_ScalarOutputNames.push_back(name);
          }else{
            if(name.compare(this->m_ScalarOutputNames.back())!=0)
              this->m_ScalarOutputNames.push_back(name);
          }
        }
      }
      eat_whitespace(config_file);
      index++;
    }
    if(fs){
      if(this->m_Parameters.back().m_Name.compare(0,1," ")==0 || this->m_Parameters.back().m_Name.empty())
      this->m_Parameters.pop_back();
    }else{
      if(this->m_ScalarOutputNames.back().compare(0,1," ")==0 || this->m_ScalarOutputNames.back().empty())
        this->m_ScalarOutputNames.pop_back();
    }
  }else{
    std::cout << "Could not open " << fileName.c_str() << std::endl;
    this->stateFlag=ERROR;
    return OTHER_ERROR;
  }
	config_file.close();
  
  return NO_ERROR;
}

/** 
 * Take a step in parameter space and calculate LiklihoodNew, PriorNew, ProposalNew,
 * and ProposalCurrent.
 */
madai::MultiModel::ErrorType 
madai::MultiModel::GetScalarOutputs(const std::vector< double > & parameters,
                                    std::vector< double > & scalars ) const
{
  double ScaleC = parameters[parameters.size()-2];
  double ScaleN = parameters[parameters.size()-1];
  std::vector<double> Temp_Theta; 
  std::vector<double> CurrentTheta;
  for(std::vector<double>::const_iterator par_it = parameters.begin(); par_it < parameters.end()-2; par_it++){
    CurrentTheta.push_back(*par_it);
  }
  Temp_Theta = m_Proposal->Iterate(CurrentTheta, ScaleN);
  unsigned int i;
  for( i = 0; i < Temp_Theta.size(); i++){
    scalars.push_back(Temp_Theta[i]);
  }
  scalars.push_back( m_Likelihood->Evaluate(Temp_Theta) ); // m_LikelihoodNew
  scalars.push_back( m_Prior->Evaluate(Temp_Theta) ); // m_PriorNew
  scalars.push_back( m_Proposal->Evaluate(CurrentTheta, Temp_Theta, ScaleC) ); // m_ProposalNew
  scalars.push_back( m_Proposal->Evaluate(Temp_Theta, CurrentTheta, ScaleN) ); // m_ProposalCurrent
	return NO_ERROR;
}

// Not implemented yet.  Should we do this numerically?
/** Get both scalar values and the gradient of the parameters. */
madai::MultiModel::ErrorType 
madai::MultiModel::GetScalarAndGradientOutputs(const std::vector< double > & parameters,
                                              const std::vector< bool > & activeParameters,
                                              std::vector< double > & scalars,
                                              unsigned int outputIndex, std::vector< double > & gradient) const 
{
	return OTHER_ERROR;
}
