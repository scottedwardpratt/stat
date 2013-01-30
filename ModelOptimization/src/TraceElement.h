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
  TraceElement(const std::vector< double > & parameterValues,
               const std::vector< double > & OutputValues );
  TraceElement(const std::vector< double > & parameterValues);
  TraceElement();
  void Reset();
  void Print();
  void VizTrace();
    
  std::vector< double > m_ParameterValues;
  std::vector< double > m_OutputValues;
  bool                  m_InTrace;
  bool                  m_Used;
    
  /** Comments may be used to store human-readable comments *or*
  record changes to state, such as changing an optimizer type,
  which parameters are activated, etc.. */
  std::vector< std::string > m_Comments;

}; // class TraceElement

} // end namespace madai

#endif // __TraceElement_h__
