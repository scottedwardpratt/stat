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
  unsigned int idx = this->GetOutputScalarIndex(scalarName);
  if (idx == static_cast< unsigned int >(-1))
    return INVALID_PARAMETER_INDEX_ERROR;
  this->m_OutputScalarToOptimize = scalarName;
  this->m_OutputScalarToOptimizeIndex = idx;
  return NO_ERROR;
}


Optimizer::ErrorType
Optimizer
::SetOutputScalarToOptimizeIndex(unsigned int idx)
{
  if (idx >= m_Model->GetNumberOfScalarOutputs())
    return INVALID_PARAMETER_INDEX_ERROR;
  this->m_OutputScalarToOptimizeIndex = idx;
  this->m_OutputScalarToOptimize = this->m_Model->GetScalarOutputNames()[idx];
  return NO_ERROR;
}

std::string Optimizer::GetOutputScalarToOptimize()
{
  return this->m_OutputScalarToOptimize;
} 
unsigned int Optimizer::GetOutputScalarToOptimizeIndex() const
{
  return this->m_OutputScalarToOptimizeIndex;
}




const std::vector< double > &
Optimizer
::GetCurrentParameters() const
{
  return m_CurrentParameters;
}


unsigned int
Optimizer
::GetOutputScalarIndex( const std::string & scalarName ) const
{
  std::vector< std::string > const & outputs = this->m_Model->GetScalarOutputNames();
  for ( unsigned int i = 0; i < outputs.size(); i++ )
  {
    if ( outputs[i] == scalarName )
      return i;
  }

  return static_cast< unsigned int >(-1); // Intentional underflow
}


unsigned int
Optimizer
::GetParameterIndex( const std::string & parameterName ) const
{
  const std::vector< Parameter > & parameters = this->m_Model->GetParameters();
  for ( unsigned int i = 0; i < m_Model->GetNumberOfParameters(); i++ )
  {
    if ( parameters[i].m_Name == parameterName )
      return i;
  }

  return static_cast< unsigned int >(-1); // Intentional underflow
}

} // end namespace madai
