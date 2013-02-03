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

#include "RegularStepGradientDescentOptimizer.h"


namespace madai {


RegularStepGradientDescentOptimizer
::RegularStepGradientDescentOptimizer( const Model *model ) :
  Optimizer( model )
{
  m_StepSize = 1.0e-3;
}


RegularStepGradientDescentOptimizer
::~RegularStepGradientDescentOptimizer()
{
}


void
RegularStepGradientDescentOptimizer
::NextIteration( Trace *trace )
{
  std::vector< bool > activeParameters( m_CurrentParameters.size() );

  // TODO - set these from the set of active parameters
  for ( unsigned int i = 0; i < activeParameters.size(); i++ ) {
    activeParameters[i] = true;
    }

  std::vector< double > scalars;
  std::vector< double > gradient;

  // TODO - set the index from the active scalar

  Model::ErrorType error =
    m_Model->GetScalarAndGradientOutputs( m_CurrentParameters, activeParameters,
                                          scalars, 0, gradient );

  trace->add( m_CurrentParameters, scalars );

  // Update the current parameters to the new position
  unsigned int activeParameter = 0;
  for ( unsigned int i = 0; i < m_CurrentParameters.size(); ++i ) {
    if ( m_ActiveParameters.find( m_Model->GetParameters()[i].m_Name ) !=
         m_ActiveParameters.end() ) {
      m_CurrentParameters[i] = -m_StepSize * gradient[activeParameter++] +
        m_CurrentParameters[i];
    }
  }
}


void
RegularStepGradientDescentOptimizer
::SetStepSize( double stepSize )
{
  m_StepSize = stepSize;
}

} // end namespace madai
