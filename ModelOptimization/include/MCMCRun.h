#ifndef __MCMCRun_h__
#define __MCMCRun_h__

#include <vector>
#include <string>
#include <map>

#include "parametermap.h"
#include "random.h"
#include <gsl/gsl_randist.h>

#include "Model.h"
#include "MultiModel.h"
#include "Trace.h"
#include "Visualization.h"
#include "Optimizer.h"

namespace madai {

class Model;
class VizHandler;
class Trace;
class TraceElement;
    
class MCMCRun : public Optimizer {
public:
  MCMCRun(const Model *in_model);
  ~MCMCRun();

  void NextIteration(Trace *trace);
        
  std::vector<double> GetRandomTheta0(int seed);
  std::vector<double> GetTheta0FromFile();

  parameterMap        m_LocalParameterMap;
  std::vector<double> m_BestParameterSet;
  std::vector<double> m_ParameterValues;
  std::vector<double> m_InitialTheta;
  int                 m_BurnIn;
  bool                m_RandomTheta0;
  bool                m_VizTrace;
  bool                m_Quiet;
  bool                m_RescaledTrace;
  bool                m_AppendTrace;
  bool                m_LogPrior;
  bool                m_LogProposal;
  bool                m_CreateTrace;
  
  VizHandler*         m_Visualizer;
  CRandom*            m_RandomNumber;
  
  std::vector<double> m_CurrentParameters;
  double              m_LikelihoodCurrent;
  double              m_LikelihoodNew;
  double              m_PriorCurrent;
  double              m_PriorNew;
  double              m_ProposalCurrent;
  double              m_ProposalNew;
  double              m_BestLikelihood;
  float               m_ScaleCurrent;
  float               m_ScaleNew;
  int                 m_AcceptCount;
  int                 m_VizCount;
  int                 m_IterationNumber;
protected:
};

} // end namespace madai

#endif // end __MCMCRun_h__