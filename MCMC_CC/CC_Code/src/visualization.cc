#ifndef __VISUALIZATION_CC__
#define __VISUALIZATION_CC__

#include "visualization.h"

using namespace std;

VizHandler::VizHandler(MCMC *mcmc_in){
	mcmc = mcmc_in;
	gnuplotpipe = popen("gnuplot -persist", "w");
	if(!gnuplotpipe){
		cout << "Gnuplot not found!" << endl;
		exit(1);
	}
	string gnuplotterm = parameter::getS(mcmc->parmap, "GNUPLOT_TERMINAL", "x11");
	string gnuplotstyle = parameter::getS(mcmc->parmap, "GNUPLOT_TRACESTYLE", "linespoints");
	
	fprintf(gnuplotpipe, "%s\n", ("set term " + gnuplotterm).c_str());
	fprintf(gnuplotpipe, "%s\n", ("set title '" + mcmc->runnickname + "'").c_str());
	fprintf(gnuplotpipe, "%s\n", "set xlabel 'iteration'");
	fprintf(gnuplotpipe, "%s\n", "set ylabel 'parameter value'");
	fflush(gnuplotpipe);
	
	header = "plot ";
	for(int i = 0; i < mcmc->ThetaList->ParamNames.size(); i++){
		header = header + "'-' w "+ gnuplotstyle + " t '"+ mcmc->ThetaList->ParamNames[i];
		if(i != mcmc->ThetaList->ParamNames.size()-1){
			header = header + "', ";
		}
	}
	header = header + "\n";
	
	paramvalues = new string[mcmc->ThetaList->ParamNames.size()];
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
	for(int i = (HighestItnReadIn % mcmc->WRITEOUT); i < mcmc->WRITEOUT; i++){
		if(mcmc->ThetaList->Theta[i]->Used){
			int currentiteration = mcmc->WRITEOUT*(mcmc->ThetaList->WriteOutCounter-1)+i+1;
			if(currentiteration > HighestItnReadIn){
				HighestItnReadIn = currentiteration;
			}
			for(int j =0; j< mcmc->ThetaList->Theta[i]->Values.size(); j++){
				if(paramvalues[j].empty()){
					paramvalues[j] = "";
				}
				ss << paramvalues[j] << mcmc->WRITEOUT*(mcmc->ThetaList->WriteOutCounter-1)+i+1 << " "\
				 << mcmc->ThetaList->Theta[i]->Values[j] << "\n";
				paramvalues[j] = ss.str();
				ss.str(string()); //clears the stringstream.
			}
		}
	}
	string plotcommand = header;
	for(int i = 0; i < mcmc->ThetaList->ParamNames.size(); i++){
		plotcommand = plotcommand + paramvalues[i]+ "e\n";
	}
	// cout << plotcommand << endl;
	fprintf(gnuplotpipe, "%s", plotcommand.c_str());
	fflush(gnuplotpipe);
}

#endif