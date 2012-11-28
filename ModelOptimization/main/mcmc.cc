#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

#include "MCMCRun.h"
#include "MultiModel.h"

int main(int argc, char ** argv){
  srand(time(NULL));
  if(argc != 2) {
    std::cerr << 
    "Useage:\n\t mcmc info_dir_path\n\n"
    "where info_dir_path is the path to the "
    "directory containing all of the configuration "
    "files needed to run the mcmc.\n\n";
    return 0;
  }
  std::string info_dir(argv[1]);
  madai::MultiModel m_model;
  m_model.LoadConfiguration(info_dir,"default");
  if(!(m_model.good())){
    std::cerr << "Something is wrong with the model" << std::endl << std::endl;
    return 0;
  }
    
  madai::MCMCRun run(&m_model);
    
  std::vector<madai::Parameter> const * parameters = &(m_model.GetParameters());
  for(int i=0; i<parameters->size();i++)
    run.ActivateParameter((*parameters)[i].m_Name);
    
  run.m_ThetaZeroPtr = &(*run.m_ThetaOutsList)[0];
  run.m_CurrentParameters = run.m_ThetaZeroPtr->m_ParameterValues;
	
	if(!run.m_ThetaZeroPtr){
    std::cout << "Error getting zeroth iteration parameters." << std::endl;
	}
	
	run.m_Likelihood_Current = m_model.m_Likelihood->Evaluate(run.m_CurrentParameters);
	run.m_Scale_Current = (rand() / double(RAND_MAX));
	//run.Proposal_Current = m_model.m_Proposal->Evaluate(*ThetaZeroPtr);
	run.m_Prior_Current = m_model.m_Prior->Evaluate(run.m_CurrentParameters);
    
  for(int j=0; j<m_model.GetNumberOfParameters(); j++)
    run.m_ParameterValues.push_back(0);
    
	run.m_AcceptCount = 0;
	for(run.m_IterationNumber =1; run.m_IterationNumber<=run.m_MaxIterations; run.m_IterationNumber++){
		run.NextIteration();
  }
	run.m_ThetaOutsList->WriteOut();
	run.m_ThetaOutsList->MakeTrace();
    
	if(run.m_CreateTrace){
		run.m_Visualizer->FinalTrace();
	}
    
	double ratio = (double)run.m_AcceptCount/((double)run.m_MaxIterations-(double)run.m_BurnIn);
  std::cout << "Accepts: " << run.m_AcceptCount << std::endl;
  std::cout << "Iterations-Burn in: " << run.m_MaxIterations-run.m_BurnIn << std::endl;
  std::cout << "Acceptance ratio: " << ratio << std::endl;
  printf("-------- Best Parameter Set, likelihood=%g -------------\n",run.m_BestLikelihood);
  std::cout << "This parameter set contains " << parameters->size() << " parameters." << std::endl;
  for(int i=0;i<parameters->size();i++)
    std::cout << (*parameters)[i].m_Name << ":\t" << run.m_BestParameterSet[i] << std::endl;

    
	//run.BestParameterSetPtr->Print();
    
  std::cout << "Done Successfully." << std::endl;
  return 0;
}