#include "MultiModel.h"

madai::MultiModel::MultiModel()
{
  this->m_Process.question = NULL;
  this->m_Process.answer = NULL;
  this->m_Emulator = NULL;
  this->stateFlag = UNINITIALIZED;
}

madai::MultiModel::MultiModel(const std::string info_dir)
{
  this->m_Process.question = NULL;
  this->m_Process.answer = NULL;
  this->m_Emulator = NULL;
  this->stateFlag = UNINITIALIZED;
  this->LoadConfiguration(info_dir);
}

madai::MultiModel::~MultiModel()
{
  if( m_Likelihood != NULL )
    delete m_Likelihood;
  if( m_Prior != NULL )
    delete m_Prior;
  if( this->m_Process.question != NULL )
    std::fclose(this->m_Process.question);
  if( this->m_Process.answer != NULL )
    std::fclose(this->m_Process.answer);
  if(m_Emulator != NULL)
    delete m_Emulator;
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
  m_ProcessPipe = parameter::getB(m_ParameterMap, "PROCESS_PIPE", false);
  m_UseEmulator = parameter::getB(m_ParameterMap, "USE_EMULATOR", false);
  m_Optimizer = parameter::getS(m_ParameterMap, "OPTIMIZER", "NOOPTIMIZER" );
   
    //============================================
	// Reading the parameters out of ranges.dat
  m_PrescaledParams = parameter::getB(m_ParameterMap, "PRESCALED_PARAMS", false);
  std::string filename = info_dir + "/ranges.dat";
  this->LoadConfigurationFile(filename);
  std::cout << "Ranges loaded" << std::endl;
    
  if(!(this->m_Parameters.empty())){
    std::vector<madai::Parameter>::const_iterator itr = this->m_Parameters.begin();
    for( itr; itr < this->m_Parameters.end(); itr++)
      std::cout << itr->m_Name << std::endl;
  }else{
    std::cout << "Parameters were not read in!" << std::endl;
    this->stateFlag = ERROR;
    return OTHER_ERROR;
  }
  
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
  
  if(m_UseEmulator && !m_ProcessPipe){
    std::cerr << "Emulator is being loaded from: " << m_DirectoryName << + "/Emulator.statfile" << "using the EmuPlusPlus.h emulator handler" << std::endl;
    m_Emulator= new ::emulator(m_DirectoryName + "/Emulator.statefile");
    std::cerr << "Emulator loaded. Test (number of params): _" << m_Emulator->number_params << "_" << std::endl;
  } else if (m_ProcessPipe){
    this->LoadProcess();
  } else {
    std::cerr << "Neither the process pipe nor use emulator flag is set to true. We can't do anything without an emulator" << std::endl;
    exit(1);
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
	
  if(std::strcmp(m_Optimizer.c_str(), "MCMC")==0){
    if(std::strcmp(m_ModelType.c_str(),"RHIC")==0){
      m_Likelihood = new LikelihoodDistribution_RHIC(this);
      m_Prior = new PriorDistribution_RHIC(this);
    }else{
      printf("Must define parameter MODEL in parameter file, or yours is unrecognized\n");
      this->stateFlag = ERROR;
      return OTHER_ERROR;
    } 
  } else {
    std::cerr << "Likelihood and Prior not loaded since the MCMC isn't being used" << std::endl;
  }
  
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
              double temp = MinR;
              MinR = MaxR;
              MaxR = temp;
            }
          }else{
            std::string line;
            std::getline(config_file, line);
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
 * Get the scalar outputs for the model
 **/
madai::MultiModel::ErrorType
madai::MultiModel::GetScalarOutputs( const std::vector< double > & parameters,
                                     std::vector< double > & scalars ) const
{
  if(m_ProcessPipe){
    std::vector< double > temp_outs;
    unsigned int index = 0;
    double * range = new double[2]();
    for( std::vector< double >::const_iterator par_it = parameters.begin(); par_it < parameters.end(); par_it++ ){
      this->GetRange(index,range);
      if(*par_it < range[0] || *par_it > range[1]){
        std::cerr << "Parameter out of bounds" << std::endl;
        return OTHER_ERROR;
      }
      std::fprintf(this->m_Process.question, "%17lf\n", *par_it);
      index++;
    }
    std::fflush(this->m_Process.question);
    double dtemp;
    for(unsigned int i = 0; i < 2*number_of_outputs; i++){
      if(1 != fscanf(this->m_Process.answer, "%lf%*c", &dtemp) ){
        std::cerr << "Interprocess communication error [cj83A]\n";
        return OTHER_ERROR;
      }
      scalars.push_back(dtemp);
    }
    return NO_ERROR;
  } else if( m_UseEmulator && !m_ProcessPipe ) {
    std::vector< double > Means;
    std::vector< double > Errors;
    m_Emulator->QueryEmulator( parameters, Means, Errors );
    if(Means.size() != number_of_outputs){
      std::cerr << "Number of observables from emulator not the number of observables read in\n";
      return OTHER_ERROR; 
    }
    for(unsigned int j = 0; j < Means.size(); j++){
      scalars.push_back(Means[j]);
      scalars.push_back(Errors[j]);
    }
    return NO_ERROR;
  } else {
    // Find some other way to fill the scalars vector
    std::cerr << "Scalars not filled in GetScalarOutputs" << std::endl;
    return OTHER_ERROR;
  }
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

// For interaction with the mcmc
madai::MultiModel::ErrorType
madai::MultiModel::GetLikeAndPrior( const std::vector< double > & parameters,
                                       double & Like,
                                       double & Prior) const
{
  std::vector< double > ModelMeans;
  std::vector< double > ModelErrors;
  if(m_UseEmulator && !m_ProcessPipe){
    m_Emulator->QueryEmulator(parameters, ModelMeans, ModelErrors);
  } else if ( m_ProcessPipe ) {
    std::vector< double > outputs;
    this->GetScalarOutputs( parameters, outputs );
    for(unsigned int i = 0; i < outputs.size(); i++){
      if(i%2 == 0){
        ModelMeans.push_back( outputs[i] );
      } else {
        ModelErrors.push_back( outputs[i] );
      }
    }
  } else {
    std::cerr << "No emulator is being used. Can't do anything without an emulator" << std::endl;
    exit(1);
  }
  
  Like = m_Likelihood->Evaluate( ModelMeans, ModelErrors );
  Prior = m_Prior->Evaluate( parameters );

  return NO_ERROR;
}

madai::MultiModel::ErrorType
madai::MultiModel::LoadProcess()
{
  std::string EmuSnapFile_Name = m_DirectoryName + "/Emulator.statefile";
  std::cerr << "Loading emulator from " << EmuSnapFile_Name << " as a running process_pipe" << std::endl;
  std::ofstream shell_script;
  std::string ssname = "begin_int_emu.sh";
  shell_script.open(ssname.c_str());
  std::string cmd = "~/local/bin/interactive_emulator interactive_mode " + EmuSnapFile_Name;
  shell_script << cmd;
  shell_script.close();
  std::string mfe = "chmod +x begin_int_emu.sh";
  std::system(mfe.c_str());
  
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
  if(EXIT_FAILURE == create_process_pipe(&(this->m_Process), argv)){
    std::cerr << "create_process_pipe returned failure.\n";
    exit(1);
  }
  for (unsigned int i = 0; i < command_line_length; i++) {
		delete argv[i];
	}
	delete[] argv;
  
	if (this->m_Process.answer == NULL || this->m_Process.question == NULL) {
		std::cerr << "create_process_pipe returned NULL fileptrs.\n";
    this->stateFlag = ERROR;
    return OTHER_ERROR;
	}
  
	discard_comments(this->m_Process.answer, '#');
	// allow comment lines to BEGIN the interactive process
  
 	unsigned int n;
	if (1 != std::fscanf(this->m_Process.answer, "%u", &n)) {
		std::cerr << "fscanf failure reading from the external process [1]\n";
    this->stateFlag = ERROR;
    return OTHER_ERROR;
	}
	if (n != number_of_parameters) {
		std::cerr << "number_of_parameters mismatch\n";
    this->stateFlag = ERROR;
    return OTHER_ERROR;
	}
	eat_whitespace(this->m_Process.answer);
	for (unsigned int i = 0; i < number_of_parameters; i++) {
		discard_line(this->m_Process.answer);
	}
	if (1 != std::fscanf(this->m_Process.answer,"%d", &n)) {
		std::cerr << "fscanf failure reading from the external process [2]\n";
    this->stateFlag = ERROR;
    return OTHER_ERROR;
	}
	if (n != (2*number_of_outputs)) {
		std::cerr << "number_of_outputs mismatch\n";
    this->stateFlag = ERROR;
    return OTHER_ERROR;
	}
	eat_whitespace(this->m_Process.answer);
	for (unsigned int i = 0; i < 2*number_of_outputs; i++) {
		discard_line(this->m_Process.answer);
	}
  return NO_ERROR;
}
