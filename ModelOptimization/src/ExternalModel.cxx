/*********************************************************************
MADAI Model Statistical Tools
Copyright 2011-2012, The University of North Carolina at Chapel Hill.

This software was written in 2011-2012 by
	Hal Canary <hal AT cs.unc.edu>
while working for the MADAI project <http://madai.us/>.

See copyright.txt for more information.
*********************************************************************/
#include <string> // std::string
#include <fstream>
#include <cstdio> // std::fgets, std::fscanf, et cetera.
#include <cstdlib> // EXIT_FAILURE
#include <cassert>

#include "ExternalModel.h"

madai::ExternalModel::ExternalModel()
{
  this->process.question = NULL; // construcotr must do this
  // so that we know whether the destructor must act.
  this->process.answer = NULL; // construcotr must do this.
  this->stateFlag = UNINITIALIZED;
}

madai::ExternalModel::ExternalModel(const std::string & configFileName)
{
  this->process.question = NULL; // construcotr must do this
  // so that we know whether the destructor must act.
  this->process.answer = NULL; // construcotr must do this.
  this->stateFlag = UNINITIALIZED;
  this->LoadConfigurationFile(configFileName);
}

madai::ExternalModel::~ExternalModel()
{
  if (this->process.question != NULL)
    std::fclose(this->process.question);
  if (this->process.answer != NULL)
    std::fclose(this->process.answer);
  // Hopefully the external process will clean itself up after its
  // stdin is closed for reading.
}


/**
 * Loads a configuration from a file.  The format of the file is
 * defined by this function.  We'll lock it down later.
 */
madai::ExternalModel::ErrorType
madai::ExternalModel::LoadConfigurationFile( const std::string fileName )
{
  this->configFileName = fileName; // keep a copy of the file name
  // just in case we need it.
  std::ifstream configFile(fileName.c_str());
  ErrorType r = this->LoadConfigurationFile(configFile);
  configFile.close();
  return r;
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


madai::ExternalModel::ErrorType
madai::ExternalModel::LoadConfigurationFile( std::istream & configFile )
{
  discard_comments(configFile, '#');
  configFile >> this->number_of_parameters;
  if (this->number_of_parameters < 1) {
  std::cerr << "bad value [YRDR]\n"; // 1g6o8duPJOAzMJgo
  this->stateFlag = ERROR;
  return OTHER_ERROR;
  }
  eat_whitespace(configFile);
  //std::cerr << this->number_of_parameters << '\n';
  this->m_Parameters.reserve(this->number_of_parameters);
  for (unsigned int i = 0; i < this->number_of_parameters; i++) {
  std::string s;
  double minimum = 0.0, maximum = 0.0;
  std::getline(configFile, s);
  configFile >> minimum >> maximum;
  if (minimum >= maximum) {
  std::cerr << "bad range [Tw4X]: " << i << ' '<<minimum << ':' << maximum<<"\n";
  this->stateFlag = ERROR;
  return OTHER_ERROR;
  }
  eat_whitespace(configFile);
  this->m_Parameters.push_back(Parameter(s, minimum, maximum));
  //std::cerr << s << '\n';
  }
  configFile >> this->number_of_outputs;
  if (this->number_of_outputs < 1) {
  std::cerr << "bad value [YRDR]\n"; // 1g6o8duPJOAzMJgo
  this->stateFlag = ERROR;
  return OTHER_ERROR;
  }
  eat_whitespace(configFile);
  //std::cerr << this->number_of_outputs << '\n';
  this->m_ScalarOutputNames.reserve(this->number_of_parameters);
  for (unsigned int i = 0; i < this->number_of_outputs; i++) {
  std::string s;
  std::getline(configFile, s);
  this->m_ScalarOutputNames.push_back(s);
  //std::cerr << s << '\n';
  }
  unsigned int command_line_length;
  configFile >> command_line_length;
  eat_whitespace(configFile);
  //	std::cerr << command_line_length << '\n';
  assert(command_line_length > 0);
  char ** argv = new char* [command_line_length + 1];
  argv[command_line_length] = NULL; // termination signal for exec*();
  for (unsigned int i = 0; i < command_line_length; i++) {
  std::string s;
  std::getline(configFile, s);
  unsigned int stringsize = s.size();
  argv[i] = new char[stringsize + 1];
  s.copy (argv[i], stringsize);
  argv[i][stringsize] = '\0'; //nil terminated c-string.
  }

  /*
				We are now done reading configuration settings.  We will now
				open a pipe to the emulator and leave it going
	*/

  /** function returns EXIT_FAILURE on error, EXIT_SUCCESS otherwise */
  if (EXIT_FAILURE == create_process_pipe(&(this->process), argv)) {
  std::cerr << "create_process_pipe returned failure.\n";
  this->stateFlag = ERROR;
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
  if (1 != std::fscanf(this->process.answer,"%u",&n)) {
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
  // read but don't do anything with parameter names.
  // THINK ABOUT where is best to define these?
  }
  if (1 != std::fscanf(this->process.answer,"%d", &n)) {
  std::cerr << "fscanf failure reading from the external process [2]\n";
  this->stateFlag = ERROR;
  return OTHER_ERROR;
  }
  if (n != this->number_of_outputs) {
  std::cerr << "number_of_outputs mismatch";
  this->stateFlag = ERROR;
  return OTHER_ERROR;
  }
  eat_whitespace(this->process.answer);
  for (unsigned int i = 0; i < this->number_of_outputs; i++) {
  discard_line(this->process.answer);
  }
  /* We are now ready to go! */
  this->stateFlag = READY;
  return NO_ERROR;
}

/** Get the valid range for the parameter at parameterIndex. */
void madai::ExternalModel::GetRange(unsigned int parameterIndex, double range[2]) const
{
  range[0] = this->m_Parameters.at(parameterIndex).m_MinimumPossibleValue;
  range[1] = this->m_Parameters.at(parameterIndex).m_MaximumPossibleValue;
}

/**
 * Get the scalar outputs from the model evaluated at x.  If an
 * error happens, the scalar output array will be left incomplete.
 */
madai::ExternalModel::ErrorType
madai::ExternalModel::GetScalarOutputs(
  const std::vector< double > & parameters,
  std::vector< double > & scalars ) const
{
  for (std::vector<double>::const_iterator par_it = parameters.begin();
       par_it < parameters.end(); par_it++ ) {
  //FIXME to do: check against parameter range.
  std::fprintf(this->process.question,"%.17lf\n", *par_it);
  }
  std::fflush(this->process.question);
  for (std::vector<double>::iterator ret_it = scalars.begin();
       ret_it < scalars.end(); ret_it++ )
    if (1 != fscanf(this->process.answer, "%lf%*c", &(*ret_it))) {
    std::cerr << "interprocess communication error [cJ83A]\n";
    return OTHER_ERROR;
    }
  return NO_ERROR;
}

// Not implemented yet.  Should we do this numerically?
/** Get both scalar values and the gradient of the parameters. */
madai::ExternalModel::ErrorType
madai::ExternalModel::GetScalarAndGradientOutputs(
  const std::vector< double > & parameters,
  const std::vector< bool > & activeParameters,
  std::vector< double > & scalars,
  unsigned int outputIndex, std::vector< double > & gradient) const
{
  return OTHER_ERROR;
}

// Not implemented yet.
// Get the likelihood and prior at the point theta
madai::ExternalModel::ErrorType
madai::ExternalModel::GetLikeAndPrior( const std::vector< double > & parameters,
                                       double & Like,
                                       double & Prior ) const
{
  return OTHER_ERROR;
}
