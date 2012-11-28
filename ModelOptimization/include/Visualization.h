#ifndef __VISUALIZATION_H__
#define __VISUALIZATION_H__

#include "MCMCRun.h"
#include <deque>

namespace madai {

class MCMCRun;
    
class VizHandler{
public:
	VizHandler(MCMCRun * mcmc_in);
	~VizHandler();

	void operator() (const std::string& command);
	void UpdateTraceFig();
	void FinalTrace();

private:

  std::string              m_GNUPlotTerm;
  std::string              m_GNUPlotStyle;
	int                      m_ThetaListSize;
	int                      m_HighestItnReadIn;
	bool                     m_MovingWindow;
	bool                     m_DensityPlot;
	MCMCRun*                 m_MCMC;
  std::string*             m_ParamValues;
  std::string              m_Header;
	deque<std::string>*      m_DequeParameterValues;
	FILE*                    m_GNUPlotPipe;
	FILE*                    m_GNUPlotMultipipe;
  std::vector<std::string> m_DensityPlotFileNames;
  std::vector<std::string> m_DensityPlotCommands;
  /**
   * The Densities array holds the density plot data. In principle the first two indices
   * should be NUM_PARAMS, and NUM_PARAMS-1, but I don't know how to do that here. The 
   * last two are the number of bins.
   **/
	int m_Densities[7][7][100][100];
	int m_Bins;
};

}

#endif