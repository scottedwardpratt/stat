#include "CHydro.h"

int main (int argc, char * const argv[]) {
	string run_name;
	if(argc!=2){
		printf("USAGE: hydromain run_name\n");
		exit(1);
	}
	run_name=string(argv[1]);
	string command = "mkdir -p output";
	system(command.c_str());
	command = "mkdir -p output/" + run_name;
	system(command.c_str());
	
	const int bListSize = 5;
	double bList[bListSize] = {0., 2.2, 3.7, 5.2, 7.0};
	
	parameterMap* pMap = new parameterMap();
	
	string parsfilename="parameters/"+run_name+"/fixed.param";
	parameter::ReadParsFromFile(*pMap,parsfilename);
	parsfilename="parameters/"+run_name+"/stats.param";
	parameter::ReadParsFromFile(*pMap,parsfilename);
	
	CHydro* mHydro;
	int status;
	
	for (int i=0;i<1;i++) {//bListSize;i++) {
		
		char bName[7];
		sprintf(bName,"/b%1.3g/",bList[i]);

		string dataRoot = string("output/") + run_name + string(bName);
			//printf("\n\nmkdir -p %s\n\n",dataRoot.c_str());
		
		command = "mkdir -p " + dataRoot;
		system(command.c_str());		
		
		parameter::set(*pMap,"HYDRO_OUTPUT_DATAROOT",dataRoot);
		parameter::set(*pMap,"GLAUBER_B",bList[i]);
			//		parameter::PrintPars(*pMap);
		
		mHydro = new CHydro(pMap);
		printf("check a\n");
		status = mHydro->runHydro();
		printf("check b\n");
		
		if (status != 0) {
			printf("\n\n*******crash generating %s.....\n\n*******aborting with %d unfinished runs!!!!!!\n\n",
				   dataRoot.c_str(),bListSize-i);
			return 1;
		}
		else 
			printf("\nSuccessfully processed %s....\n",dataRoot.c_str());
		delete mHydro;
	}
	
	return 0;
}

