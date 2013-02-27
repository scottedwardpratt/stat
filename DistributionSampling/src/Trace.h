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

#ifndef __Trace_h__
#define __Trace_h__

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

#include <cmath>
#include <map>

#include "Parameter.h"
#include "parametermap.h"

#include "TraceElement.h"

namespace madai {

/** \class Trace
 *
 * Traces contain all the state of the distribution sampling required to
 * replay it. */
class Trace {
public:
  Trace();
  virtual ~Trace();

  /** Add an entry from parameter, output values, and log-likelihood. */
  void Add( const std::vector< double > & parameterValues,
            const std::vector< double > & outputValues,
            double logLikelihood );

  /** Add an entry from parameter and output values.
   *
   * Sets log-likelihood to 0.0. */
  void Add( const std::vector< double > & parameterValues,
            const std::vector< double > & outputValues );

  /** Add an entry from parameter values alone.
   *
   * \todo It seems like you should always have to record an output
   * from a model. */
  void Add( const std::vector< double > & parameterValues );

  /** Get the number of entries in the Trace. */
  unsigned int GetSize() const;

  TraceElement & operator[]( unsigned int index );
  const TraceElement & operator[]( unsigned int index ) const;

  void Write( std::ostream & o ) const;  

  /*
    Assert:
      FOR ALL i < this->m_TraceElements.size():
        this->m_TraceElements[i].m_ParameterValues.size() == params.size()
        this->m_TraceElements[i].m_OutputValues.size() == outputs.size()
  */
  void WriteHead( std::ostream & o, 
                  const std::vector< Parameter > & params ) const;

  void WriteHead( std::ostream & o,
                  const std::vector< Parameter > & params,
                  const std::vector< std::string > & outputs) const;
  void PrintDataToFile( const std::vector< Parameter > & params );
  void WriteOut( const std::vector< Parameter > & params );
  void MakeTrace();

protected:
  std::vector< TraceElement > m_TraceElements;

}; // class Trace

} // namespace madai

#endif // __Trace_h__
