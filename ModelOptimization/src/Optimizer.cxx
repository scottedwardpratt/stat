#include "Optimizer.h"

#include "Parameter.h"

namespace madai {

Optimizer
::Optimizer( const Model *model )
{
  m_Model = model;

  // Initialize the vector for current parameters.
  m_CurrentParameters.resize( m_Model->GetNumberOfParameters() );

  // Activate all parameters by default.
  const std::vector< Parameter > parameterDescriptions =
    m_Model->GetParameters();
  for ( unsigned int i = 0; i < parameterDescriptions.size(); ++i )
  {
    this->ActivateParameter( parameterDescriptions[i].m_Name );
  }
}


Optimizer
::~Optimizer()
{
}


void
Optimizer
::ActivateParameter( const std::string & parameterName )
{
  m_ActiveParameters.insert( parameterName );
}


void
Optimizer
::DeactivateParameter( const std::string & parameterName )
{
  m_ActiveParameters.erase( parameterName );
}


Optimizer::ErrorType
Optimizer
::SetParameterValue( const std::string & parameterName, double value )
{
  unsigned int parameterIndex = this->GetParameterIndex( parameterName );

  if ( parameterIndex == static_cast< unsigned int >( -1 ) )
    {
    // Error. Parameter not found.
    return INVALID_PARAMETER_INDEX_ERROR;
    }

  m_CurrentParameters[parameterIndex] = value;

  // TODO - set dirty flag

  return NO_ERROR;
}


Optimizer::ErrorType
Optimizer
::SetOutputScalarToOptimize( const std::string & scalarName )
{
  m_OutputScalarToOptimize = scalarName;
}


const std::vector< double > &
Optimizer
::GetCurrentParameters() const
{
  return m_CurrentParameters;
}


unsigned int
Optimizer
::GetParameterIndex( const std::string & parameterName )
{
  std::vector< Parameter > parameters = m_Model->GetParameters();
  for ( unsigned int i = 0; i < m_Model->GetNumberOfParameters(); i++ )
  {
    if ( parameters[i].m_Name == parameterName )
    {
      return i;
    }
  }

  return static_cast< unsigned int >(-1); // Intentional underflow
}

} // end namespace madai
