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

	gnuplotpipe = popen("gnuplot -persist", "w");
	if(!gnuplotpipe){
		cout << "Gnuplot not found!" << endl;
		exit(1);
	}
	gnuplotterm = parameter::getS(mcmc->local_parmap, "GNUPLOT_TERMINAL", "x11");
	gnuplotstyle = parameter::getS(mcmc->local_parmap, "GNUPLOT_TRACESTYLE", "linespoints");
	
	fprintf(gnuplotpipe, "%s\n", ("set term " + gnuplotterm).c_str());
	// fprintf(gnuplotpipe, "%s\n", ("set title '" + mcmc->runnickname + "'").c_str());
	fprintf(gnuplotpipe, "%s\n", "set xlabel 'iteration'");
	fprintf(gnuplotpipe, "%s\n", "set ylabel 'parameter value'");
	fprintf(gnuplotpipe, "%s\n", "set key out vert right top");
	fflush(gnuplotpipe);
	
	header = "plot ";
	for(int i = 0; i < mcmc->ThetaList->ParamNames.size(); i++){
		header = header + "'-' w " + gnuplotstyle + " t '" + mcmc->ThetaList->ParamNames[i]+ "'";
		if(i != mcmc->ThetaList->ParamNames.size()-1){
			header = header + ", ";
		}
	}
	//header = header + ", '-' w "+ gnuplotstyle + " t 'Likelihood'"; //Adding the Likelihood to the plot
	header = header + "\n";
	
	//cout << "The gnuplot header is: \n" << header << endl;

	paramvalues = new string[mcmc->ThetaList->ParamNames.size()];
	DequeParameterValues = new deque<string>[mcmc->ThetaList->ParamNames.size()];
	MovingWindow = true; //This should be made into an option which is set in mcmc.param
	HighestItnReadIn = 0;
}

VizHandler::~VizHandler(){
	fprintf(gnuplotpipe, "exit\n");
	pclose(gnuplotpipe);
}

//not sure if I'm going to use this, it should allow for something like:
//VizHandler plotter(this);
//plotter("plot cos(x)");
void VizHandler::operator() (const string& command) {
	fprintf(gnuplotpipe, "%s\n", command.c_str());
	fflush(gnuplotpipe);
}

void VizHandler::UpdateTraceFig(){
	stringstream ss;
	string plotcommand = header;
	//string likelihood;
	
	if(MovingWindow){
		int DequeSize = 500;
		int k;
		for(int i = 0; i < ThetaListSize; i++){
			if(mcmc->ThetaList->Theta[i]->Used && !(mcmc->ThetaList->Theta[i]->InTrace)){
				for(int j = 0; j< mcmc->ThetaList->Theta[i]->Values.size(); j++){
					k = mcmc->WRITEOUT*mcmc->ThetaList->WriteOutCounter + i + 1;
					ss << mcmc->WRITEOUT*mcmc->ThetaList->WriteOutCounter + i + 1 << " " << mcmc->ThetaList->Theta[i]->Values[j] << "\n";
					//cout << mcmc->WRITEOUT*mcmc->ThetaList->WriteOutCounter + i + 1 << " " << mcmc->ThetaList->Theta[i]->Values[j] << "\n";
					if(DequeParameterValues[j].size() > DequeSize){
						DequeParameterValues[j].pop_front();
						DequeParameterValues[j].push_back(ss.str());
					}
					else{
						DequeParameterValues[j].push_back(ss.str());
					}
					ss.str(string()); //clears the stringstream.
					//cout << i << " " << j << " " << ThetaListSize << " " << mcmc->ThetaList->Theta[i]->Values.size() << endl;
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
					ss << paramvalues[j] << mcmc->WRITEOUT*mcmc->ThetaList->WriteOutCounter + i + 1 << " "\
					<< mcmc->ThetaList->Theta[i]->Values[j] << "\n";
					paramvalues[j] = ss.str();
					ss.str(string()); //clears the stringstream.
				}
				mcmc->ThetaList->Theta[i]->VizTrace();
			}
		}
		for(int i = 0; i < mcmc->ThetaList->ParamNames.size(); i++){
			plotcommand = plotcommand + paramvalues[i]+ "e\n";
		}
	}
	
	//cout << "The plot command is: \n" << plotcommand << endl;
	fprintf(gnuplotpipe, "%s", plotcommand.c_str());
	fflush(gnuplotpipe);
}

void VizHandler::FinalTrace(){
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
}

#endif