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

#ifndef __TraceElement_h__
#define __TraceElement_h__

#include <string>
#include <vector>


namespace madai {


/** \class TraceElement
 *
 * Individual entry in a Trace.
 */
class TraceElement {
public:
  //@{
  /**
     If LogLikelihood is not set, it defaults to 0.0.  If
     output_values is not specified, it defaults to an empty vector.
  */
  TraceElement(const std::vector< double > & parameterValues,
               const std::vector< double > & OutputValues,
               double LogLikelihood);

  TraceElement(const std::vector< double > & parameterValues,
               const std::vector< double > & OutputValues );

  TraceElement(const std::vector< double > & parameterValues);

  TraceElement();

  /** Clear all the parameter values, output values, and comments, and
   * set the log likelihood to 0.0. */
  void Reset();

  /** Returns true if there are any parameter values, output values,
   * or comments. */
  bool IsValid() const;

  /** a set of Parameters values from a model. */
  std::vector< double > m_ParameterValues;
  /** the model output that correspond to m_ParameterValues in some model. */
  std::vector< double > m_OutputValues;
  /**
     Given some set of field observations, the likelihood is the
     Likelihood that m_ParameterValues is the ground truth.
     this->m_LogLikelihood is (ln(C * Likelihood)) for some unknown
     constant C.
  */
  double                m_LogLikelihood;

  /** Comments may be used to store human-readable comments *or*
  record changes to state, such as changing an optimizer type,
  which parameters are activated, etc.. */
  std::vector< std::string > m_Comments;

}; // class TraceElement

} // end namespace madai

#endif // __TraceElement_h__
