#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

#include "MCMCRun.h"
#include "TestModel.h"

/***
 * Idea for a regression test:
 * Do a run with the mcmc with a fixed seed. Each run should be identicle with a
 * fixed seed since the sequence of random numbers from rand() will always be the
 * same for the same seed. If the trace changes between runs, then something is
 * wrong with the mcmc code.
 ***/

bool dis_comments(std::FILE * fp, char comment_character){
  int c = std::getc(fp);
  if((c == EOF) || std::ferror(fp)) {
    std::cerr << "premature endof file:(\n";
    return false;
  }
  while ( c == comment_character ) {
    static int buffersize = 1024;
    char buffer[buffersize];
    std::fgets(buffer, buffersize, fp);
    c = std::getc(fp);
  }
  if(EOF == std::ungetc(c, fp)) {
    std::cerr << "ungetc error :(\n";
    return false;
  }
  return true;
}

int main(int argc, char ** argv){
  srand(0);
  if(argc != 2) {
    "Useage:\n\t mcmc info_dir_path\n\n"
    "where info_dir_path is the path to the "
    "directory containing all of the configuration "
    "files needed to run the mcmc.\n\n";
    return 0;
  }
  std::string info_dir(argv[1]);
  madai::TestModel t_model;
  t_model.LoadConfiguration(info_dir);
  if(!(t_model.good())){
    std::cerr << "Something is wrong with the model\n\n";
    return 0;
  }
  
  madai::MCMCRun run(&t_model, info_dir);
  
  std::vector< madai::Parameter > const * parameters = &(t_model.GetParameters());
  for(int i = 0; i < parameters->size(); i++)
    run.ActivateParameter( (*parameters)[i].m_Name );
    
  madai::Trace trace(info_dir, "default" );
  if(run.m_BurnIn == 0){
    trace.add(run.m_InitialTheta);
  }
  
  for(int j = 0; j < t_model.GetNumberOfParameters(); j++)
    run.m_ParameterValues.push_back(0);
    
  run.m_AcceptCount = 0;
  for(run.m_IterationNumber = 1; run.m_IterationNumber < trace.m_MaxIterations; run.m_IterationNumber++){
    run.NextIteration(&trace);
  }
  if(trace.m_CurrentIteration!=0){
    trace.WriteOut(t_model.GetParameters());
  }
  trace.MakeTrace();
  
  // At this point the trace for the new run has been created.
  // Now want to compare it to a previous run.
  
  FILE* fp1 = fopen("test_trace.dat", "r");
  std::string trace_file_name;
  trace_file_name = trace.m_TraceDirectory.c_str();
  trace_file_name += "/trace.dat";
  FILE* fp2 = fopen(trace_file_name.c_str(), "r");
  if(fp1 == NULL){
    std::cerr << "Error opening test_trace.dat" << std::endl;
    return 0;
  }
  if(fp2 == NULL){
    std::cerr << "Error opening " << trace.m_TraceDirectory+"/trace.dat" << std::endl;
    return 0;
  }
  double* set1 = new double[parameters->size()+5]();
  double* set2 = new double[parameters->size()+5]();
  while(!feof(fp1) && !feof(fp2)){
    int k, iter_num1, iter_num2;
    dis_comments(fp1, '#');
    dis_comments(fp2, '#');
    fscanf(fp1, "%d,", &iter_num1);
    fscanf(fp2, "%d,", &iter_num2);
    if(iter_num1!=iter_num2){
      std::cerr << "Iteration number error :(" << std::endl;
      std::cerr << "Iter1: " << iter_num1;
      std::cerr << "Iter2: " << iter_num2;
      exit(1);
    }
    for(k = 0; k < parameters->size(); k++){
      fscanf(fp1, "%lf,", &set1[k]);
      fscanf(fp2, "%lf,", &set2[k]);
      if(set1[k]!=set2[k]){
        std::cerr << "Different values for parameter " << k << " on iteration " << iter_num1 << std::endl;
        if(iter_num1 == 0){
          std::cerr << "MCMC Core Test Failure: 0th iteration error. Initial theta needs to be the same" << std::endl;
          exit(1);
        } else {
          std::cerr << "MCMC Core Test Failure: Something must have gone wrong in taking a step in MCMCRun.cxx" << std::endl;
          exit(1);
        }
      }
    }
    // Scan in Likelihood, Prior, Proposal (Current and New), and alpha.
    for(int l = k; l < k+5; l++){
      fscanf(fp1, "%lf,", &set1[l]);
      fscanf(fp2, "%lf,", &set2[l]);
    }
    fscanf(fp1,"\n");
    fscanf(fp2,"\n");
    
    if(set1[k]!=set2[k]){
      std::cerr << "Different values for the Likelihood on iteration " << iter_num1 << "\n\n";
      std::cerr << "MCMC Core Test Failure: Model/Distribution Error." << std::endl;
      std::cerr << "Either TestModel.cxx's GetScalarOutputs or LikelihoodTest.cxx's Evaluate has changed" << std::endl;
      exit(1);
    }
    k++;
    if(set1[k]!=set2[k]){
      std::cerr << "Different values for the Prior on iteration " << iter_num1 << "\n\n";
      std::cerr << "MCMC Core Test Failure: Prior Distribution Error." << std::endl;
      std::cerr << "Check PriorTest.cxx's Evaluate" << std::endl;
      exit(1);
    }
    k++;
    if(set1[k]!=set2[k] || set1[k+1]!=set2[k+1]){
      std::cerr << "Different values for the Current Proposal on iteration " << iter_num1 << "\n\n";
      std::cerr << "MCMC Core Test Failure: MCMCRun Step Taking Error." << std::endl;
      std::cerr << "Check MCMCRun.cxx's EvaluateProposal and TakeStep" << std::endl;
      exit(1);
    }
    k+=2;
    if(set1[k]!=set2[k]){
      std::cerr << "Different values for alpha on iteration " << iter_num1 << "\n\n";
      std::cerr << "MCMC Core Test Failure: MCMC Error." << std::endl;
      std::cerr << "Areas that might be responsible:" << std::endl;
      std::cerr << "1) MCMCRun::TakeStep" << std::endl;
      std::cerr << "2) MCMCRun::EvaluateProposal" << std::endl;
      std::cerr << "3) MCMCRun::NextIteration" << std::endl;
    }
  }
  if(!feof(fp1)){
    std::cerr << "MCMC Core Test Failure: Size Error." << std::endl;
    std::cerr << "End of test_trace.dat has not been reached. Make sure run is of the right length." << std::endl;
  }
  std::cerr << "MCMC Core Test Passed:" << std::endl;
  std::cerr << "New run is identicle to test_trace.dat" << std::endl;
  
  return 0;
}