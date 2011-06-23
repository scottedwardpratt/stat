#ifndef __Trace_h_
#define __Trace_h_

#include <vector>

namespace madai {

class TraceElement {
public:
  std::vector< double >      m_ParameterValues;
  std::vector< double >      m_OutputValues;

  /** Comments may be used to store human-readable comments *or*
  record changes to state, such as changing an optimizer type,
  which parameters are activated, etc.. */
  std::vector< std::string > m_Comments;

}; // class TraceElement



/** Traces contain all the state of the optimization required to
 * replay it. */
class Trace {
public:
  Trace() {};
  virtual ~Trace();

  // TODO - add accessor methods

protected:
  std::vector< TraceElement > m_TraceElements;

  std::vector< std::string > m_ParameterNames;

}; // class Trace

} // namespace madai

#endif // __Trace_h_
