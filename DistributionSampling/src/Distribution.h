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

#ifndef __Distribution_h__
#define __Distribution_h__

#include "Random.h"

namespace madai {

/** \class Distribution
 *
 * Base class for distributions. */
class Distribution {
public:
  Distribution();
  virtual ~Distribution();

  virtual double GetLogProbabilityDensity(double value) const = 0;
  virtual double GetProbabilityDensity(double value) const = 0;
  virtual double GetPercentile(double percentile) const = 0;
  virtual double GetSample(madai::Random & r) const = 0;

protected:

};

} // namespace madai

#endif
