#include "TraceElement.h"


namespace madai {


TraceElement
::TraceElement( const std::vector< double > & parameter_values,
                const std::vector< double > & output_values ) :
  m_ParameterValues( parameter_values ),
  m_OutputValues( output_values )
{
}


TraceElement
::TraceElement( const std::vector< double > & parameter_values) :
  m_ParameterValues( parameter_values ),
  m_Used( true )
{
}


void
TraceElement
::Reset()
{
  m_ParameterValues.clear();
  m_OutputValues.clear();
  m_Comments.clear();
  m_Used = false;
  m_InTrace = false;
}


TraceElement
::TraceElement()
{
  m_Used=false;
}


void
TraceElement
::VizTrace()
{
  if ( !m_InTrace ) {
    m_InTrace=true;
  }
}

} // end namespace madai
