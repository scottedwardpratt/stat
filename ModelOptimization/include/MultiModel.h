#ifndef __MultiModel_h__
#define __MultiModel_h__

#include <vector>
#include <string>
#include <map>
#include <math.h>
#include <iostream>
#include <fstream>
#include "parametermap.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <sys/stat.h>
#include "EmuPlusPlus/EmuPlusPlus.h"
#include "process_pipe.h"
#include "Model.h"

namespace madai {
    
class MultiModel : public Model {
private:
  unsigned int number_of_parameters, number_of_outputs;
  process_pipe process;
  typedef enum {
    UNINITIALIZED,
    READY,
    ERROR
  } internal_state;
  internal_state stateFlag;
        
public:
  std::vector<bool> m_LogParam;
  bool              m_PrescaledParams;
  std::string       m_ModelType;
        
  bool good() { return (this->stateFlag == READY);}
  
  MultiModel();
  MultiModel(std::string info_dir);
  MultiModel(std::string info_dir, std::string configuration);
  virtual ~MultiModel();

  /**
   * Load a configuration from multiple files. The format of these
   * files is defined by the following functions
   **/
  virtual ErrorType LoadConfigurationFile( const std::string fileName );
  virtual ErrorType LoadConfiguration(const std::string info_dir);
  virtual ErrorType LoadConfiguration(const std::string info_dir, const std::string configuration);
  /** 
   * Load a process pipe for emulation handling
   **/
  virtual ErrorType LoadProcess();
  /** 
   * Get the scalar outputs from the model evaluated at x.  If an
   * error happens, the scalar output array will be left incomplete.
   **/
  virtual ErrorType GetScalarOutputs( const std::vector< double > & parameters,
                                     std::vector< double > & scalars ) const;
        
  // Not implemented yet.  Should we do this numerically?
  /** Get both scalar values and the gradient of the parameters. */
  virtual ErrorType GetScalarAndGradientOutputs(const std::vector< double > & parameters,
                                                const std::vector< bool > & activeParameters,
                                                std::vector< double > & scalars,
                                                unsigned int outputIndex, std::vector< double > & gradient) const;
};

} // end namespace madai

#endif // end __RHICModel_h__