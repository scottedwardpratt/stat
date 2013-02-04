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

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "MCMCRun.h"
#include "Gaussian2DModel.h"

/***
 * Idea for a regression test:
 * Do a large run with the mcmc using a 2D gaussian as the model. Binning the
 * point in parameter space should approximate a gaussian shape.
 ***/
 
bool skip_comments(std::FILE * fp, char comment_character){
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
  srand(time(NULL));
  if(argc != 2){
    std::cerr <<
    "Useage:\n\t mcmc info_dir_path\n\n"
    "where info_dir_path is the path to the "
    "directory containing all of the configuration "
    "files needed to run the mcmc.\n\n";
    return 0;
  }
  std::string info_dir(argv[1]);
  madai::Gaussian2DModel g2d_model;
  madai::MCMCRun run(&g2d_model, info_dir);
  
  std::vector< madai::Parameter > const * parameters = &(g2d_model.GetParameters());
  for(int i = 0; i < parameters->size(); i++)
    run.ActivateParameter( (*parameters)[i].m_Name );
  
  madai::Trace trace(info_dir, "default" );
  if(run.m_BurnIn == 0){
    trace.add(run.m_InitialTheta);
  }
  
  for(int j = 0; j < g2d_model.GetNumberOfParameters(); j++)
    run.m_ParameterValues.push_back(0);
    
  run.m_AcceptCount = 0;
  for(run.m_IterationNumber = 1; run.m_IterationNumber < trace.m_MaxIterations; run.m_IterationNumber++){
    run.NextIteration(&trace);
  }
  if(trace.m_CurrentIteration!=0){
    trace.WriteOut(g2d_model.GetParameters());
  }
  trace.MakeTrace();
  
  // At this point the run has completed. Now I want to bin the points in the trace.
  std::string trace_file_name;
  trace_file_name = trace.m_TraceDirectory.c_str();
  trace_file_name +="/trace.dat";
  FILE* fp = fopen(trace_file_name.c_str(), "r");
  if(fp == NULL){
    std::cerr << "Error opening trace.dat [1]" << std::endl;
    return 0;
  }
  double range[2][2];
  double x, y;
  int iter;
  range[0][0] = range[1][0] = DBL_MAX;
  range[0][1] = range[1][1] = -DBL_MAX;
  while(!feof(fp)){
    skip_comments(fp, '#');
    fscanf(fp, "%d,%lf,%lf\n", &iter,&x,&y);
    if(x < range[0][0]){
      range[0][0]=x;
    } else if (x > range[0][1]){
      range[0][1]=x;
    }
    if(y < range[1][0]){
      range[1][0]=y;
    } else if (y > range[1][1]){
      range[1][1]=y;
    }
  }
  range[0][0]-=0.5;
  range[1][0]-=0.5;
  range[0][1]+=0.5;
  range[1][1]+=0.5;
  
  fclose(fp);
  FILE* tfile = fopen(trace_file_name.c_str(), "r");
  if(tfile == NULL){
    std::cerr << "Error opening trace.dat [2]" << std::endl;
    return 0;
  }
  int** densities = new int*[100]();
  for(unsigned int n = 0; n < 100; n++){
    densities[n] = new int[100]();
  }
  while(!feof(tfile)){
    skip_comments(tfile, '#');
    fscanf(tfile, "%d,%lf,%lf\n", &iter, &x, &y);
    // scale x and y to [0,100)
    x = 100*(x-range[0][0])/(range[0][1]-range[0][0]);
    y = 100*(y-range[1][0])/(range[1][1]-range[1][0]);
    densities[int(x)][int(y)]++;
  }
  // Calculate integrated densities (rho(x) and rho(y))
  int* rhox = new int[100]();
  int* rhoy = new int[100]();
  for(unsigned int i = 0; i < 100; i++){
    for(unsigned int j = 0; j < 100; j++){
      rhox[i]+=densities[i][j];
      rhoy[j]+=densities[i][j];
    }
  }
  
  // Print densities to files
  std::ofstream o1;
  o1.open("DensityPlot.txt");
  std::ofstream ox;
  ox.open("DensityPlotX.txt");
  std::ofstream oy;
  oy.open("DensityPlotY.txt");
  double binx, biny, temp, sum=0, absdev=0, MeanX, MeanY;
  std::vector< double > tempv;
  binx = (range[0][1] - range[0][0]) / 100;
  biny = (range[1][1] - range[1][0]) / 100;
  for(unsigned int k = 0; k < 100; k++){
    ox << range[0][0]+(double(k)+0.5)*binx << " " << double(rhox[k])/double(trace.m_MaxIterations) << std::endl;
    oy << range[1][0]+(double(k)+0.5)*biny << " " << double(rhoy[k])/double(trace.m_MaxIterations) << std::endl;
    for(unsigned int l = 0; l < 100; l++){
      tempv.clear();
      temp = range[0][0] + (double(k)+0.5)*binx;
      o1 << temp << " ";  // Print x value of center of bin
      tempv.push_back(temp);
      temp = range[1][0]+(double(l)+0.5)*biny;
      o1 << temp << " ";  // Print y value of center of bin
      tempv.push_back(temp);
      temp = double(densities[k][l])/double(trace.m_MaxIterations);
      o1 << temp << std::endl;  // Print density
      std::vector< double > Prob;
      g2d_model.GetScalarOutputs(tempv, Prob);
      g2d_model.GetMeans( MeanX, MeanY );
      Prob[0]/=(2*MeanX*MeanY);
      sum+=(Prob[0]-temp);
      absdev+=abs(Prob[0]-temp);
    }
  }
  sum/=10000;
  std::cerr << "The average deviation from the model is: " << sum << std::endl;
  std::cerr << "Average Absolute Deviation: " << (absdev/10000) << std::endl;
  o1.close();
  ox.close();
  oy.close();
  delete rhox;
  delete rhoy;
  for(unsigned int q = 0; q < 100; q++){
    delete densities[q];
  }
  delete densities;
  // Open pipe to create plots of rho(x) and rho(y)
  FILE* XPlot = popen("gnuplot -persist", "w");
  bool PlotDensities = true;
  if(!XPlot){
    std::cerr << "GNUPlot not found :(" << std::endl;
    PlotDensities = false;
    //exit(1);
  }
  FILE* YPlot = popen("gnuplot -persist", "w");
  if(!YPlot){
    std::cerr << "GNUPlot not found :(" << std::endl;
    PlotDensities = false;
    //exit(1);
  }
  if(PlotDensities){
    std::string command1;
    std::string command2;
    fprintf(XPlot, "%s\n", "set term x11");
    fprintf(XPlot, "%s\n", "set xlabel \"X\"");
    fprintf(XPlot, "%s\n", "set ylabel \"Density\"");
    fprintf(XPlot, "%s\n", "set title \"Density Over X\"");
    fprintf(XPlot, "%s\n", "set key out vert right top");
    command1 = "plot \"DensityPlotX.txt\" t \"rho(x)\" w lines";
    fprintf(XPlot, "%s\n", command1.c_str());
    fprintf(YPlot, "%s\n", "set term x11");
    fprintf(YPlot, "%s\n", "set ylabel \"Density\"");
    fprintf(YPlot, "%s\n", "set xlabel \"Y\"");
    fprintf(YPlot, "%s\n", "set title \"Density Over Y\"");
    fprintf(YPlot, "%s\n", "set key out vert right top");
    command2 = "plot \"DensityPlotY.txt\" t \"rho(y)\" w lines";
    fprintf(YPlot, "%s\n", command2.c_str());
  }
  
  return 0;
}
