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

#ifndef __PriorDistribution_h__
#define __PriorDistribution_h__


#include "Distribution.h"


namespace madai {

class PriorDistribution : public Distribution {
public:
  PriorDistribution();
  virtual ~PriorDistribution();
  virtual double Evaluate(std::vector<double> Theta);
};

} // end namespace madai


#endif // __PriorDistribution_h__

