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

#include "DiceModel.h"


namespace madai {


DiceModel
::DiceModel()
{
  this->m_StateFlag = UNINITIALIZED;
}


DiceModel
::DiceModel( const std::string info_dir )
{
  this->m_StateFlag = UNINITIALIZED;
  this->LoadConfigurationFile(info_dir);
}


Model
::ErrorType
DiceModel
::LoadConfigurationFile( const std::string info_dir )
{
  m_DirectoryName = info_dir;
  m_ConfigFile    = m_DirectoryName+"/DiceConfig.param";
  parameter::ReadParsFromFile( m_ParameterMap, m_ConfigFile.c_str() );
  m_Sum                = parameter::getB( m_ParameterMap, "SUM_DICE", true );
  m_Distinguishable    = parameter::getB( m_ParameterMap, "DISTINGUISHABLE", false );
  m_NumberOfParameters = parameter::getI( m_ParameterMap, "DICE", 5 );
  
  // Factor for calculating the likelihood and add parameters
  int m_Denom = 1;
  std::string par_name;
  for ( unsigned int i = 0; i < m_NumberOfParameters; i++ ) {
    std::stringstream ss;
    ss << i;
    std::string addon = ss.str();
    m_Denom *= 6;
    this->AddParameter( par_name+addon, 0.5, 6.5 );
  }
  
  this->m_StateFlag = READY;
  
  return Model::NO_ERROR;
}


Model
::ErrorType
DiceModel
::GetScalarOutputs( const std::vector< double > & parameters,
                    std::vector< double > & scalars ) const
{
  scalars.clear();
  for ( unsigned int i = 0; i < parameters.size(); i++ ) {
    double temp;
    if ( parameters[i] < 0.5 || parameters[i] > 6.5 ) {
      std::cerr << i << "th parameter out of range" << std::endl;
      return OTHER_ERROR;
    }
    if ( parameters[i] <= 1.5 ) {
      temp = 1.0;
    } else if ( parameters[i] <= 2.5 ) {
      temp = 2.0;
    } else if ( parameters[i] <= 3.5 ) {
      temp = 3.0;
    } else if ( parameters[i] <= 4.5 ) {
      temp = 4.0;
    } else if ( parameters[i] <= 5.5 ) {
      temp = 5.0;
    } else if ( parameters[i] <= 6.5 ) {
      temp = 6.0;
    }
    scalars.push_back( temp );
  }
  
  return NO_ERROR;
}


Model
::ErrorType
DiceModel
::GetScalarAndGradientOutputs( const std::vector< double > & parameters,
                               const std::vector< bool > & activeParameters,
                               std::vector< double > & scalars,
                               unsigned int outputIndex,
                               std::vector< double > & gradient) const
{
  return Model::OTHER_ERROR;
}


Model
::ErrorType
DiceModel
::GetLikeAndPrior( const std::vector< double > & parameters,
                   double & Like,
                   double & Prior) const
{
  if ( m_Sum ) {
    int sum = 0, s2 = 0;
    for ( unsigned int l = 0; l < parameters.size(); l++ ) {
      sum += parameters[l];
    }
    // Now calculate probability of getting that sum
    for ( unsigned int iter = 0; iter < floor((sum-parameters.size())/6.0); iter++ ) {
      double temp = 0;
      int iter2;
      temp += parameters.size();
      if ( iter % 2 == 1 ) {
        temp *= -1;
      }
      // (k-6*iter-1)!
      for ( iter2 = 0; iter2 < (sum-6*iter-1); iter2++ ) {
        temp *= (iter2 + 1);
      }
      // 1/(k-6*iter-n)!
      for ( iter2 = 0; iter2 < (sum-6*iter-parameters.size()); iter2++ ) {
        temp /= double( iter2 + 1 );
      }
      // 1/(n-iter)!
      for ( iter2 = 0; iter2 < (parameters.size()-iter); iter2++ ) {
        temp /= double( iter2 + 1 );
      }
      // 1/iter!
      for ( iter2 = 0; iter2 < iter; iter2++ ) {
        temp /= double( iter2 + 1 );
      }
      Like += temp;
    }
    Like /= double( m_Denom );
    Prior = 1;
  } else {
    if ( m_Distinguishable ) {
      Like = 1.0 / double( m_Denom );
      Prior = 1;
    } else {
      // Get outputs
      std::vector< double > outputs;
      int* bins = new int[6]();
      this->GetScalarOutputs( parameters, outputs );
      // Find total number of possible permutations of the dice
      int TotalPerms;
      for ( unsigned int m = 0; m < outputs.size(); m++ ) {
        TotalPerms *= (m + 1);
      }
      // Bin the dice values
      for ( unsigned int j = 0; j < outputs.size(); j++ ) {
        for ( unsigned int k = 0; k < 6; k++ ) {
          if ( outputs[j] == double(k + 1) ) {
            bins[k]++;
            break;
          }
        }
      }
      // Eliminate repeated permutations
      int factor=1;
      for ( unsigned int n = 0; n < 6; n++ ) {
        if ( bins[n] > 1 ) {
          for ( unsigned int index = 0; index < bins[n]; index++ ) {
            factor *= (index + 1);
          }
        }
      }
      TotalPerms /= factor;
      Like = double( TotalPerms ) / double( m_Denom );
      Prior = 1;
    }
  }
  Like = log( Like );
  
  return Model::NO_ERROR;
}

} // end namespace madai
