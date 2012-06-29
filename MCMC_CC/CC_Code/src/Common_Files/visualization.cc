#ifndef __VISUALIZATION_CC__
#define __VISUALIZATION_CC__

#include "visualization.h"
#include <deque>

using namespace std;

VizHandler::VizHandler(MCMCRun *mcmc_in){
	mcmc = mcmc_in;
	ThetaListSize = mcmc->WRITEOUT;
	
	if(!mcmc){
		cout << "MCMC not loaded!" << endl;
		exit(1);
	}

	if(mcmc->VIZTRACE){
		// Setting up the pipe for the trace plot
		gnuplotpipe = popen("gnuplot -persist", "w");
		if(!gnuplotpipe){
			cout << "Gnuplot not found!" << endl;
			exit(1);
		}
		gnuplotterm = parameter::getS(mcmc->local_parmap, "GNUPLOT_TERMINAL", "x11");
		gnuplotstyle = parameter::getS(mcmc->local_parmap, "GNUPLOT_TRACESTYLE", "linespoints");
		
		fprintf(gnuplotpipe, "%s\n", ("set term " + gnuplotterm).c_str());
		fprintf(gnuplotpipe, "%s\n", "set xlabel 'iteration'");
		fprintf(gnuplotpipe, "%s\n", "set ylabel 'parameter value'");
		fprintf(gnuplotpipe, "%s\n", "set key out vert right top");
		fflush(gnuplotpipe);
	
		// Setting up the pipe for the density plots
		/*gnuplotmultipipe = popen("gnuplot -persist", "w");
		if(!gnuplotmultipipe){
			cout << "Gnuplot not found!" << endl;
			exit(1);
		}
	
		stringstream ss;
		ss << "set term " << gnuplotterm << "\n set multiplot\n";
	
		fprintf(gnuplotmultipipe, "%s\n", ss.str().c_str());
		fflush(gnuplotmultipipe);
		*/
	}
	// The header for the trace plot
	header = "plot ";
	for(int i = 0; i < mcmc->ThetaList->ParamNames.size(); i++){
		header = header + "'-' w " + gnuplotstyle + " t '" + mcmc->ThetaList->ParamNames[i]+ "'";
		if(i != mcmc->ThetaList->ParamNames.size()-1){
			header = header + ", ";
		}
	}
	header = header + "\n";

	//Check if DensityPlots file exists
	string tempfile;
	struct stat st;
	tempfile = mcmc->mcmcconfig->dir_name + "/DensityPlots";
	if(stat(tempfile.c_str(), &st)!=0){
		string command = "mkdir -p " + tempfile;
		system(command.c_str());
	}

	// The plot commands for the density plots
	for(int i = 0; i < mcmc->ThetaList->ParamNames.size(); i++){
		//cout << mcmc->ThetaList->ParamNames[i] << endl;
		stringstream command;
		command << "clear\nset size square " << 1./float(mcmc->ThetaList->ParamNames.size()-2) << "," << 1./float(mcmc->ThetaList->ParamNames.size()-2) << "\n";
		for(int j = i+1; j < mcmc->ThetaList->ParamNames.size(); j++){
			//cout << "     " << mcmc->ThetaList->ParamNames[j] << endl
			//Now we mke our 2d plots
			string filename = mcmc->mcmcconfig->dir_name + "/DensityPlots/" + mcmc->ThetaList->ParamNames[i] + "_" + mcmc->ThetaList->ParamNames[j];
			DensityPlotFileNames.push_back(filename);
			command << "set bmargin 0.001\n set tmargin 0.001\n set lmargin 0.001\n set rmargin 0.001\nunset key\nunset xtics\nunset ytics\n";
			command << "\nset origin " << (float(i)*(1./float(mcmc->ThetaList->ParamNames.size()))) << "," << (float(j-1)*(1./float(mcmc->ThetaList->ParamNames.size()))) << "\n";
			command << "set xlabel \'" << mcmc->ThetaList->ParamNames[i] << "\'\nset ylabel \'" << mcmc->ThetaList->ParamNames[j] << "\'\nset view map\nsplot \'" << DensityPlotFileNames.back() << ".txt\' matrix with image\n";
			DensityPlotCommands.push_back(command.str());
                        //cout << command.str();
			command.str("");
			if(mcmc->mcmcconfig->APPEND_TRACE){
				fstream densityfile;
				string line,val;
				stringstream ss2;
				filename = filename + ".txt";
				densityfile.open(filename.c_str());
				if(densityfile){
					for(int l = 0; l < 100; l++){
						getline(densityfile,line,'\n');
						ss2 << line;
						for(int k = 0; k < 100; k++){
							getline(ss2,val,' ');
							Densities[i][j][l][k]=atoi(val.c_str());
						}
					}
				}
			}
			else{
				for(int l = 0; l < 100; l++){
					for(int k = 0; k < 100; k++){
						Densities[i][j][l][k] = 0;
					}
				}
			}
		}
	}

	//cout << "The gnuplot header is: \n" << header << endl;

	paramvalues = new string[mcmc->ThetaList->ParamNames.size()];
	DequeParameterValues = new deque<string>[mcmc->ThetaList->ParamNames.size()];
	MovingWindow = true; //This should be made into an option which is set in mcmc.param
	DensityPlot = true; //Again, this should probably be set in a parameter file
	HighestItnReadIn = 0;
	BINS = 100;
}

VizHandler::~VizHandler(){
	if(mcmc->VIZTRACE){
		fprintf(gnuplotpipe, "exit\n");
		pclose(gnuplotpipe);
		//fprintf(gnuplotmultipipe, "exit\n");
		//pclose(gnuplotmultipipe);
	}
}

//not sure if I'm going to use this, it should allow for something like:
//VizHandler plotter(this);
//plotter("plot cos(x)");
void VizHandler::operator() (const string& command) {
	if(mcmc->VIZTRACE){
		fprintf(gnuplotpipe, "%s\n", command.c_str());
		fflush(gnuplotpipe);
	}
}

void VizHandler::UpdateTraceFig(){
	stringstream ss;
	string plotcommand = header;
	
	if(MovingWindow){
		int DequeSize = 500;
		for(int i = 0; i < ThetaListSize; i++){
			if(mcmc->ThetaList->Theta[i]->Used && !(mcmc->ThetaList->Theta[i]->InTrace)){
				for(int j = 0; j< mcmc->ThetaList->Theta[i]->Values.size(); j++){
					if(mcmc->mcmcconfig->RESCALED_TRACE){
						ss << mcmc->WRITEOUT*mcmc->ThetaList->WriteOutCounter + i + 1 << " " << (mcmc->ThetaList->Theta[i]->Values[j]-mcmc->mcmcconfig->Min_Ranges[j])/(mcmc->mcmcconfig->Max_Ranges[j]-mcmc->mcmcconfig->Min_Ranges[j])<< "\n";
					}
					else{
						ss << mcmc->WRITEOUT*mcmc->ThetaList->WriteOutCounter + i + 1 << " " << mcmc->ThetaList->Theta[i]->Values[j] << "\n";
					}
					if(DequeParameterValues[j].size() > DequeSize){
						DequeParameterValues[j].pop_front();
						DequeParameterValues[j].push_back(ss.str());
					}
					else{
						DequeParameterValues[j].push_back(ss.str());
					}
					ss.str(string()); //clears the stringstream.
				}
				mcmc->ThetaList->Theta[i]->VizTrace();
			}
		}
		for(int i = 0; i < mcmc->ThetaList->ParamNames.size(); i++){
			for(int j =0; j < DequeParameterValues[i].size(); j++){
				plotcommand = plotcommand + (DequeParameterValues[i])[j];
			}
			plotcommand = plotcommand + "e\n";
			//cout << plotcommand << endl;
		}
	}else{
		if(!mcmc){
			cout << "MCMC pointer not found." << endl;
		}
		for(int i = 0; i < ThetaListSize; i++){
			if(mcmc->ThetaList->Theta[i]->Used && !(mcmc->ThetaList->Theta[i]->InTrace)){
				// count++;
				for(int j =0; j< mcmc->ThetaList->Theta[i]->Values.size(); j++){
					if(paramvalues[j].empty()){
						paramvalues[j] = "";
					}
					if(mcmc->mcmcconfig->RESCALED_TRACE){
						ss << mcmc->WRITEOUT*mcmc->ThetaList->WriteOutCounter + i + 1 << " " << (mcmc->ThetaList->Theta[i]->Values[j]-mcmc->mcmcconfig->Min_Ranges[j])/(mcmc->mcmcconfig->Max_Ranges[j]-mcmc->mcmcconfig->Min_Ranges[j])<< "\n";
					}
					else{
						ss << mcmc->WRITEOUT*mcmc->ThetaList->WriteOutCounter + i + 1 << " " << mcmc->ThetaList->Theta[i]->Values[j] << "\n";
					}
					paramvalues[j] = ss.str();
					ss.str(string()); //clears the stringstream.
				}
				mcmc->ThetaList->Theta[i]->VizTrace();
			}
		}
		for(int i = 0; i < mcmc->ThetaList->ParamNames.size(); i++){
			plotcommand = plotcommand + paramvalues[i]+ "e\n";
			//cout << paramvalues[i] << endl;
		}
	}
	
	if(mcmc->VIZTRACE){
		//cout << "The plot command is: \n" << plotcommand << endl;
		fprintf(gnuplotpipe, "%s", plotcommand.c_str());
		fflush(gnuplotpipe);
	}

	if(DensityPlot){
		int index = 0;
		for(int i = 0; i < mcmc->ThetaList->ParamNames.size(); i++){
			for(int j = i+1; j < mcmc->ThetaList->ParamNames.size(); j++){
				Densities[i][j][int(mcmc->ParamValues[i]*100)][int(mcmc->ParamValues[j]*100)]++;
				ofstream outputfile;
				string filename = DensityPlotFileNames[index] + ".txt";
				outputfile.open(filename.c_str());
				for(int k = 0; k < BINS; k++){
					for(int l = 0; l < BINS; l++){
						outputfile << Densities[i][j][k][l] << " ";
					}
					outputfile << endl;
				}
				outputfile.close();
				//cout << DensityPlotCommands[index].c_str() << endl;
				//fprintf(gnuplotmultipipe, "%s", DensityPlotCommands[index].c_str());
				//fflush(gnuplotmultipipe);
				index++;
			}
		}
	}
}

void VizHandler::FinalTrace(){
	if(mcmc->VIZTRACE){
		stringstream ss;
		
		fprintf(gnuplotpipe,"%s\n","set datafile separator ','");
		fflush(gnuplotpipe);
		ss << "plot '"<< mcmc->tracedir << "/trace.dat' using 1:2 w " << gnuplotstyle << " t '" << mcmc->ThetaList->ParamNames[0] << "', ";
		for(int i = 1; i < mcmc->ThetaList->ParamNames.size(); i++){
			ss<< "'' u 1:"<< i+2 <<" w "<< gnuplotstyle << " t '"<< mcmc->ThetaList->ParamNames[i];
			if(i != mcmc->ThetaList->ParamNames.size()-1){
				ss << "', ";
			}
		}
		string gnuplotcmd = ss.str() + "\n";
		
		fprintf(gnuplotpipe, "%s", gnuplotcmd.c_str());
		fflush(gnuplotpipe);
	
		ss.str(string()); //clears the stringstream
	
		ss << "set term postscipt\n" << "set ouput " << mcmc->tracedir << "/trace.ps\n" << "replot\n" << "set term x11\n"; //This last may need to be changed to "set term win" on windows machines
		fprintf(gnuplotpipe, "%s", gnuplotcmd.c_str());
		fflush(gnuplotpipe);
	}
}

#endif
