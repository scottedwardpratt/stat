#include "MultiModel.h"

madai::MultiModel::MultiModel(){
  this->stateFlag = UNINITIALIZED;
}

madai::MultiModel::MultiModel(const std::string info_dir){
  this->stateFlag = UNINITIALIZED;
  this->LoadConfiguration(info_dir);
}

madai::MultiModel::MultiModel(const std::string info_dir, const std::string configuration){
  this->stateFlag = UNINITIALIZED;
  this->LoadConfiguration(info_dir, configuration);
}

madai::MultiModel::~MultiModel(){}

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

/*
madai::RHICModel::ErrorType
madai::RHICModel::LoadConfiguration(std::string info_dir){
    ErrorType r = this->LoadConfiguration(info_dir, "default");
    return r;
}
 */

madai::MultiModel::ErrorType
madai::MultiModel::LoadConfiguration(std::string info_dir){
  m_ConfigurationName = "default";
	m_DirectoryName = info_dir;
	m_ParameterFile = info_dir+"/defaultpars/";
  std::cout << "There should be something here ->" << m_ParameterFile << "<-" << std::endl;
  std::cout << "In config: " << m_ParameterFile << std::endl;
	m_ParameterFileName = m_ParameterFile + "/mcmc.param";
  std::cout << "Reading in " << m_ParameterFileName << std::endl;
    
	parameter::ReadParsFromFile(m_ParameterMap, m_ParameterFileName.c_str());
	m_LogLike = parameter::getB(m_ParameterMap, "LOGLIKE", true);
	/*LOGPRIOR = parameter::getB(parmap, "LOGPRIOR", true);
	LOGPROPOSAL = parameter::getB(parmap, "LOGPROPOSAL", true);
	CREATE_TRACE = parameter::getB(parmap, "CREATE_TRACE", true);*/
	//m_SuppressErrors = parameter::getB(m_ParameterMap, "SUPPRESS_ERRORS", false);
  m_ProcessPipe = parameter::getB(m_ParameterMap, "PROCESS_PIPE", false);
  m_ModelType = parameter::getS(m_ParameterMap,"MODEL","NOMODEL");
   
    //============================================
	// Reading the parameters out of ranges.dat
  m_PrescaledParams = parameter::getB(m_ParameterMap, "PRESCALED_PARAMS", false);
  std::string filename = info_dir + "/ranges.dat";
  this->LoadConfigurationFile(filename);
  std::cout << "Ranges loaded" << std::endl;
    
  if(!(this->m_Parameters.empty())){
    std::vector<madai::Parameter>::const_iterator itr = this->m_Parameters.begin();
    for(itr++;itr<this->m_Parameters.end();itr++)
      std::cout << itr->m_Name << " " << itr->m_MinimumPossibleValue << " " << itr->m_MaximumPossibleValue << std::endl;
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
    for(itr;itr<this->m_ScalarOutputNames.end();itr++)
      std::cout << *itr << " ";
    std::cout << std::endl;
  }else{
    std::cout << "Observables were not read in!" << std::endl;
    this->stateFlag=ERROR;
    return OTHER_ERROR;
  }
    
  std::vector<std::string> temp_logparam = parameter::getVS(m_ParameterMap, "LOG_PARAMETERS", "");
	
	for(int i =0; i<temp_logparam.size(); i++){
		if(std::strcmp(temp_logparam[i].c_str(), "true") == 0 || std::strcmp(temp_logparam[i].c_str(), "True") == 0){
			m_LogParam.push_back(true);
		}else if(std::strcmp(temp_logparam[i].c_str(), "false") == 0 || std::strcmp(temp_logparam[i].c_str(), "False") == 0){
			m_LogParam.push_back(false);
		}else{
      std::cout << "Unrecognized LogParam value " << temp_logparam[i] << std::endl;
			this->stateFlag=ERROR;
      return OTHER_ERROR;
		}
	}
    
  if(m_ProcessPipe)
    this->LoadProcess();
    
    /*cout << "_____________________" << endl;
     cout << "parameterfile: " << parameterfile << endl;
     cout << "MODEL: " << MODEL << endl;
     cout << "dir_name: " << dir_name << endl;
     cout << "parameterfile: " << parameterfile << endl;
     cout << "configname: " << configname << endl;
     cout << "parameter_file_name: " << parameter_file_name << endl;
     cout << "EmulatorParams: " << EmulatorParams << endl;
     cout << "---------------------" << endl;*/
    
	// cout << "stuff done." << endl;
	
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

madai::MultiModel::ErrorType
madai::MultiModel::LoadConfiguration(std::string info_dir, std::string configuration){
  m_ConfigurationName = configuration;
	m_DirectoryName = info_dir;
	m_ParameterFile = info_dir+"/defaultpars/";
  std::cout << "There should be something here ->" << m_ParameterFile << "<-" << std::endl;
  std::cout << "In config: " << m_ParameterFile << std::endl;
	m_ParameterFileName = m_ParameterFile + "/mcmc.param";
  std::cout << "Reading in " << m_ParameterFileName << std::endl;
	
	parameter::ReadParsFromFile(m_ParameterMap, m_ParameterFileName.c_str());
	m_LogLike = parameter::getB(m_ParameterMap, "LOGLIKE", true);
	/*LOGPRIOR = parameter::getB(parmap, "LOGPRIOR", true);
	LOGPROPOSAL = parameter::getB(parmap, "LOGPROPOSAL", true);
	CREATE_TRACE = parameter::getB(parmap, "CREATE_TRACE", true);*/
	//m_SuppressErrors = parameter::getB(m_ParameterMap, "SUPPRESS_ERRORS", false);
  m_ProcessPipe = parameter::getB(m_ParameterMap, "PROCESS_PIPE", false);
  m_ModelType = parameter::getS(m_ParameterMap,"MODEL","NOMODEL");
	
	//============================================
	// Reading the parameters out of ranges.dat
	m_PrescaledParams = parameter::getB(m_ParameterMap, "PRESCALED_PARAMS", false);
  std::string filename = info_dir + "/ranges.dat";
  this->LoadConfigurationFile(filename);
  std::cout << "Ranges loaded" << std::endl;
    
  if(!(this->m_Parameters.empty())){
    std::vector<madai::Parameter>::const_iterator itr = this->m_Parameters.begin();
    for(itr++;itr<this->m_Parameters.end();itr++)
    std::cout << itr->m_Name << " ";
  }else{
    std::cout << "The parameters were not read in." << std::endl;
    this->stateFlag=ERROR;
    return OTHER_ERROR;
  }
  std::cout << std::endl;
    
  std::string observables_filename = info_dir + "/pcanames.dat";
  this->LoadConfigurationFile(observables_filename);
    
  std::cout << "The emulated observables are: " << std::endl;
	if(!(this->m_ScalarOutputNames.empty())){
    std::vector<std::string>::const_iterator itr = this->m_ScalarOutputNames.begin();
    for(itr;itr<this->m_ScalarOutputNames.end();itr++)
      std::cout << *itr << " ";
    std::cout << std::endl;
  }else{
    std::cout << "The observables were not read in." << std::endl;
    this->stateFlag=ERROR;
    return OTHER_ERROR;
  }
    
	//====================================================
	
  std::vector<std::string> temp_logparam = parameter::getVS(m_ParameterMap, "LOG_PARAMETERS", "");
	
	for(int i =0; i<temp_logparam.size(); i++){
		if(std::strcmp(temp_logparam[i].c_str(), "true") == 0 || std::strcmp(temp_logparam[i].c_str(), "True") == 0){
			m_LogParam.push_back(true);
		}else if(std::strcmp(temp_logparam[i].c_str(), "false") == 0 || std::strcmp(temp_logparam[i].c_str(), "False") == 0){
			m_LogParam.push_back(false);
		}else{
      std::cout << "Unrecognized LogParam value " << temp_logparam[i] << std::endl;
			this->stateFlag=ERROR;
      return OTHER_ERROR;
		}
	}
    
  if(m_ProcessPipe)
    this->LoadProcess();
    
    /*cout << "_____________________" << endl;
     cout << "parameterfile: " << parameterfile << endl;
     cout << "MODEL: " << MODEL << endl;
     cout << "dir_name: " << dir_name << endl;
     cout << "parameterfile: " << parameterfile << endl;
     cout << "configname: " << configname << endl;
     cout << "parameter_file_name: " << parameter_file_name << endl;
     cout << "EmulatorParams: " << EmulatorParams << endl;
     cout << "---------------------" << endl;*/
    
	// cout << "stuff done." << endl;
    
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
    if(num_inputs<1){
      std::cerr << "Number of inputs is < 1" << std::endl;
      this->stateFlag = ERROR;
      return OTHER_ERROR;
    }
    bool fs, ro;
    double MinR, MaxR;
    std::string name, type;
        
    // Check for type of file
    if(fileName==(m_DirectoryName + "/pcanames.dat")){
      fs = false; //File switch to determine how the data is read in
      ro = true; //Tells the program whether to keep reading a line or not
      this->number_of_outputs = num_inputs;
      this->m_ScalarOutputNames.reserve(num_inputs);
    }else if(fileName==(m_DirectoryName+ "/ranges.dat")){
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

madai::MultiModel::ErrorType
madai::MultiModel::LoadProcess(){
  std::string EmuSnapFile_Name = m_DirectoryName + "/Emulator.statefile";
  std::cerr << "Loading emulator from " << EmuSnapFile_Name << std::endl;
  std::ofstream shell_script;
  std::string ssname = "begin_int_emu.sh";
  shell_script.open(ssname.c_str());
  std::string cmd = "~/local/bin/interactive_emulator interactive_mode " + EmuSnapFile_Name;
  shell_script << cmd;
  shell_script.close();
  /*std::string cmd2 = "chmod +x " + ssname;
  std::system(cmd2.c_str());*/
  unsigned int command_line_length=1;
  char ** argv = new char* [command_line_length + 1];
  argv[command_line_length] = NULL;
  unsigned int stringsize = ssname.size();
  argv[0] = new char[stringsize+1];
  ssname.copy(argv[0], stringsize);
  argv[0][stringsize] = '\0';
    
	/*
     Command for loading the emulator is set. We will noew open a pipe
     to the emulator and leave it going
     */
    
	/** function returns EXIT_FAILURE on error, EXIT_SUCCESS otherwise */
  if(EXIT_FAILURE == create_process_pipe(&(this->process), argv)){
    std::cerr << "create_process_pipe returned failure.\n";
    this->stateFlag=ERROR;
    return OTHER_ERROR;
  }
  for (unsigned int i = 0; i < command_line_length; i++) {
		delete argv[i];
	}
	delete[] argv;
    
	if (this->process.answer == NULL || this->process.question == NULL) {
		std::cerr << "create_process_pipe returned NULL fileptrs.\n";
		this->stateFlag = ERROR;
		return OTHER_ERROR;
	}
    
	discard_comments(this->process.answer, '#');
	// allow comment lines to BEGIN the interactive process
    
 	unsigned int n;
	if (1 != std::fscanf(this->process.answer, "%u", &n)) {
		std::cerr << "fscanf failure reading from the external process [1]\n";
		this->stateFlag = ERROR;
		return OTHER_ERROR;
	}
	if (n != this->number_of_parameters) {
		std::cerr << "number_of_parameters mismatch\n";
		this->stateFlag = ERROR;
		return OTHER_ERROR;
	}
	eat_whitespace(this->process.answer);
	for (unsigned int i = 0; i < this->number_of_parameters; i++) {
		discard_line(this->process.answer);
	}
	if (1 != std::fscanf(this->process.answer,"%d", &n)) {
		std::cerr << "fscanf failure reading from the external process [2]\n";
		this->stateFlag = ERROR;
		return OTHER_ERROR;
	}
	if (n != (2*(this->number_of_outputs))) {
		std::cerr << "number_of_outputs mismatch";
		this->stateFlag = ERROR;
		return OTHER_ERROR;
	}
	eat_whitespace(this->process.answer);
	for (unsigned int i = 0; i < (2*(this->number_of_outputs)); i++) {
		discard_line(this->process.answer);
	}
	return NO_ERROR;

}

/** 
 * Get the scalar outputs from the model evaluated at x.  If an
 * error happens, the scalar output array will be left incomplete.
 */
madai::MultiModel::ErrorType 
madai::MultiModel::GetScalarOutputs(const std::vector< double > & parameters,
                                    std::vector< double > & scalars ) const
{
  unsigned int index=0;
  double * range = new double[2]();
  for(std::vector<double>::const_iterator par_it = parameters.begin(); par_it < parameters.end(); par_it++){
    this->GetRange(index,range);
    if((*par_it)<range[0] || (*par_it)>range[1]){
      std::cerr << "Parameter out of bounds" << std::endl;
      return OTHER_ERROR;
    }
    std::fprintf(this->process.question, "%.17lf\n", *par_it);
    index++;
  }
  std::fflush(this->process.question);
  double dtemp;
  for(unsigned int i = 0; i<(2*(this->number_of_outputs)); i++){
    if(1!=fscanf(this->process.answer, "%lf%*c", &dtemp)){
      scalars.push_back(dtemp);
      std::cerr << "interprocess communication error [cj83A]n";
      return OTHER_ERROR;
    }
    scalars.push_back(dtemp);
  }
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
