#ifndef __VISUALIZATION_CC__
#define __VISUALIZATION_CC__

#include "visualization.h"

using namespace std;

VizHandler::VizHandler(MCMC *mcmc_in){
	mcmc = mcmc_in;
	ThetaListSize = mcmc->WRITEOUT;
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
	fprintf(gnuplotpipe, "%s\n", "set key out vert right top");
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
	// cout << "Update trace." << endl;
	int count = 0;
	stringstream ss;
	// cout << "Above for loop." << endl;
	if(!mcmc){
		cout << "MCMC pointer not found." << endl;
	}
	for(int i = 0; i < ThetaListSize; i++){
		if(mcmc->ThetaList->Theta[i]->Used && !(mcmc->ThetaList->Theta[i]->InTrace)){
			count++;
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
	string plotcommand = header;
	for(int i = 0; i < mcmc->ThetaList->ParamNames.size(); i++){
		plotcommand = plotcommand + paramvalues[i]+ "e\n";
	}
	// cout << plotcommand << endl;
	
	fprintf(gnuplotpipe, "%s", plotcommand.c_str());
	fflush(gnuplotpipe);
	// cout << "Done. Added " << count << " iterations to trace. " << endl;
}

#endif