#include "RHICModel.h"

madai::RHICModel::RHICModel()
{
  this->m_Process.question = NULL;
  this->m_Process.answer = NULL;
  this->m_Emulator = NULL;
  this->stateFlag = UNINITIALIZED;
}

madai::RHICModel::RHICModel(const std::string info_dir)
{
  this->m_Process.question = NULL;
  this->m_Process.answer = NULL;
  this->m_Emulator = NULL;
  this->stateFlag = UNINITIALIZED;
  this->LoadConfiguration(info_dir);
}

madai::RHICModel::~RHICModel()
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

/**
 * Get the scalar outputs for the model
 **/
madai::RHICModel::ErrorType
madai::RHICModel::GetScalarOutputs( const std::vector< double > & parameters,
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
madai::RHICModel::ErrorType 
madai::RHICModel::GetScalarAndGradientOutputs(const std::vector< double > & parameters,
                                              const std::vector< bool > & activeParameters,
                                              std::vector< double > & scalars,
                                              unsigned int outputIndex, std::vector< double > & gradient) const 
{
	return OTHER_ERROR;
}

// For interaction with the mcmc
madai::RHICModel::ErrorType
madai::RHICModel::GetLikeAndPrior( const std::vector< double > & parameters,
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

madai::RHICModel::ErrorType
madai::RHICModel::LoadProcess()
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

madai::RHICModel::ErrorType
madai::RHICModel::LoadDistributions()
{
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
  m_Likelihood = new LikelihoodDistribution_RHIC(this);
  m_Prior = new PriorDistribution_RHIC(this);
}
