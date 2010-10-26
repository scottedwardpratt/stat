#include "CHydro.h"

int main (int argc, char * const argv[]) {

	parameterMap* pMap = new parameterMap();
	for (int i=1;i<argc;i++) 
		parameter::ReadParsFromFile(*pMap, argv[i]);
		
//	parameter::PrintPars(*pMap);

	CHydro* mHydro = new CHydro(pMap);
	int status = mHydro->runHydro();

}

