/*=========================================================================
 *
 *  Copyright The University of North Carolina at Chapel Hill
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

#include "Model.h"

namespace madai {

Model
::Model()
{
  
}


Model
::~Model()
{
}


bool
Model
::IsReady() const
{
  return ( m_StateFlag == READY );
}


unsigned int
Model
::GetNumberOfParameters() const
{
  return static_cast<unsigned int>(m_Parameters.size());
}


const std::vector< Parameter > &
Model
::GetParameters() const
{
  return m_Parameters;
}


unsigned int
Model
::GetNumberOfScalarOutputs() const
{
  return static_cast<unsigned int>(m_ScalarOutputNames.size());
}


const std::vector< std::string > &
Model
::GetScalarOutputNames() const
{
  return m_ScalarOutputNames;
}


void
Model
::GetRange( unsigned int parameterIndex, double range[2] ) const
{
  range[0] = this->m_Parameters.at(parameterIndex).m_MinimumPossibleValue;
  range[1] = this->m_Parameters.at(parameterIndex).m_MaximumPossibleValue;
}


void
Model
::AddParameter( const std::string & name,
                double minimumPossibleValue,
                double maximumPossibleValue )
{
  m_Parameters.push_back( Parameter(name, minimumPossibleValue, maximumPossibleValue) );
}


void
Model
::AddScalarOutputName( const std::string & name )
{
  m_ScalarOutputNames.push_back( name );
}


} // end namespace madai
