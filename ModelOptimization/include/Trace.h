/*********************************************************************
MADAI Model Statistical Tools
Copyright 2011-2012, The University of North Carolina at Chapel Hill.

This software was written in 2011-2012 by 
	Cory Quammen <cquammen AT cs.unc.edu>
	Russell Taylor <taylorr AT cs.unc.edu>
	Scott Pratt <pratt AT nscl.msu.edu>
	Kevin Novak <novakkev AT msu.edu>
	Hal Canary <hal AT cs.unc.edu>
while working for the MADAI project <http://madai.us/>.

See copyright.txt for more information.
*********************************************************************/
#ifndef __Trace_h_
#define __Trace_h_

#include <vector>
#include <string>
#include <iostream>

#include "Parameter.h"

namespace madai {

class TraceElement {
public:
	TraceElement(const std::vector< double > & parameterValues,
							 const std::vector< double > & OutputValues );
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
  virtual ~Trace() {}
	CParameterList *parlist;

	unsigned int length() const;
	void add(const std::vector< double > & parameterValues,
					 const std::vector< double > & OutputValues );
	TraceElement & operator[](unsigned int idx);
	const TraceElement & operator[](unsigned int idx) const;

	void write(std::ostream & o) const;  

	/*
    Assert:
	    FOR ALL i < this->m_TraceElements.size():
	      this->m_TraceElements[i].m_ParameterValues.size() == params.size()
	      this->m_TraceElements[i].m_OutputValues.size() == outputs.size()
	*/
	void writeHead(std::ostream & o, 
	               const std::vector<madai::Parameter> & params,
	               const std::vector<std::string> & outputs) const;  
protected:
  std::vector< TraceElement > m_TraceElements;
  std::vector< std::string > m_ParameterNames;

}; // class Trace

} // namespace madai

#endif // __Trace_h_
