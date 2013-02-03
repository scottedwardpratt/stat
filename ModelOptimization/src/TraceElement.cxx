/*=========================================================================
 *
 *  Copyright (c) 2010-2012 The University of North Carolina at Chapel Hill
 *  All rights reserved.
 *
 *  Licensed under the MADAI Software License. You may obtain a copy of
 *  this license at
 *
 *         https://madai-public.cs.unc.edu/software/license/
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

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
