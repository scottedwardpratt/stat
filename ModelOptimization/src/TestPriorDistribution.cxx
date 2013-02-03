#include "TestPriorDistribution.h"

#include "Model.h"

namespace madai {

TestPriorDistribution
::TestPriorDistribution( Model * in_Model )
{
  m_Model = in_Model;
  m_SepMap = parameter::getB( m_Model->m_ParameterMap, "PRIOR_PARAMETER_MAP", false );

  if ( m_SepMap ) {
    std::string parmapfile = m_Model->m_DirectoryName + "/parameters/prior.param";
    m_ParameterMap = new parameterMap;
    parameter::ReadParsFromFile( *m_ParameterMap, parmapfile );
    //parameter::ReadParsFromFile(parmap, parameter_file_name);
  } else {
    m_ParameterMap = &(m_Model->m_ParameterMap);
  }
}


double
TestPriorDistribution
::Evaluate( std::vector< double > Theta ) {
  double mean  = parameter::getD( *m_ParameterMap, "PRIOR_MEAN", -3.7372 );
  double sigma = parameter::getD( *m_ParameterMap, "PRIOR_SIGMA", 1.6845 );
  double temp = 0.0;
  bool Found = false;

  std::vector< Parameter > const * parameters = &(m_Model->GetParameters());
  for ( int i = 0; i < parameters->size(); i++ ) {
    if ( (*parameters)[i].m_Name.compare( 0, 1, "SIGMA" ) == 0 ) {
      if ( !Found ) {
        temp = Normal( log( Theta[i] ), mean, sigma );
        Found = true;
      } else {
        std::cerr << "In RHIC_PCA_PRIOR::Evaluate; Duplicate parameter names found." << std::endl;
        exit( 1 );
      }
    }
  }
  if ( Found ) {
    return temp;
  } else {
    std::cerr << "SIGMA parameter not found!" << std::endl;
    std::cerr << "Will return the prior as 1.0" << std::endl;

    return 1.0;
  }
  //return Normal(log(Theta.GetValue("SIGMA")), mean, sigma);
}

} // end namespace madai
