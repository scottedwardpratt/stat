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
#ifndef __Trace_h__
#define __Trace_h__

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "math.h"
#include <map>

#include "Parameter.h"
#include "parametermap.h"

namespace madai {
    
class Trace;
class TraceElement;

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


/** \class Trace
 *
 * Traces contain all the state of the optimization required to
 * replay it. */
class Trace {
public:
  Trace() {};
  Trace(const std::string info_dir, 
        const std::string configuration);
  virtual ~Trace() {}

  unsigned int length() const;
  void add(const std::vector< double > & parameterValues,
           const std::vector< double > & OutputValues );
  void add(const std::vector< double > & parameterValues);
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
                 const std::vector< Parameter > & params) const;
  void writeHead(std::ostream & o,
                 const std::vector< Parameter > & params,
                 const std::vector< std::string > & outputs) const;
  void PrintDataToFile(const std::vector< Parameter > & params);
  void WriteOut(const std::vector< Parameter > & params);
  void MakeTrace();
  std::vector< std::string > GetParNames();
    
  std::string  m_TraceDirectory;
  int          m_Writeout;
  int          m_MaxIterations;
  int          m_WriteOutCounter;
  int          m_CurrentIteration;
  bool         m_AppendTrace;
  parameterMap m_TraceParameterMap;

protected:

  std::vector< TraceElement > m_TraceElements;
  std::vector< std::string >  m_ParameterNames;

}; // class Trace

} // namespace madai

#endif // __Trace_h__
