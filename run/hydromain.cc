#include "CHydro.h"

int main (int argc, char * const argv[]) {
	string qualifier;
	if(argc!=2){
		printf("USAGE: hydromain qualifier\n");
		exit(1);
	}
	qualifier=argv[1];
	parameterMap* pMap = new parameterMap();
	string parsfilename="parameters/"+qualifier+"/fixed.param";
	parameter::ReadParsFromFile(*pMap,parsfilename);
	parsfilename="parameters/"+qualifier+"/stats.param";
	parameter::ReadParsFromFile(*pMap,parsfilename);
//	parameter::PrintPars(*pMap);
	CHydro* mHydro = new CHydro(pMap);
	int status = mHydro->runHydro();
	return status;
}

