#ifndef __Parameter_h_
#define __Parameter_h_

#include <string>

namespace madai {

class Parameter {
public:
  Parameter() {};
  virtual ~Parameter() {};

  std::string m_Name;
  double      m_MinimumPossibleValue;
  double      m_MaximumPossibleValue;

}; // end class Parameter

} // end namespace madai

#endif // __Parameter_h_
