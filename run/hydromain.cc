#include "CHydro.h"
#include "qualifier.h"

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
	
	parameterMap* pMap = new parameterMap();
	
	string parsfilename="parameters/"+run_name+"/fixed.param";
	parameter::ReadParsFromFile(*pMap,parsfilename);
	parsfilename="parameters/"+run_name+"/stats.param";
	parameter::ReadParsFromFile(*pMap,parsfilename);
	
	CHydro* mHydro;
	int status;
	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.dat");
	
	for (int iqual=0;iqual<qualifiers.nqualifiers;iqual++) {
		
		qualifiers.SetPars(pMap,iqual);
		string dataRoot = string("output/") + run_name + "/"+qualifiers.qualifier[iqual]+"/";
			//printf("\n\nmkdir -p %s\n\n",dataRoot.c_str());
		parameter::set(*pMap,"HYDRO_OUTPUT_DATAROOT",dataRoot);		
		command = "mkdir -p " + dataRoot;
		system(command.c_str());		
		
			//		parameter::PrintPars(*pMap);
		
		mHydro = new CHydro(pMap);
		status = mHydro->runHydro();
		break;
		
		if (status != 0) {
			printf("\n\n*******crash generating %s.....\n\n*******aborting with %d unfinished runs!!!!!!\n\n",
				   dataRoot.c_str(),qualifiers.nqualifiers-iqual);
			return 1;
		}
		else 
			printf("\nSuccessfully processed %s....\n",dataRoot.c_str());
		delete mHydro;
	}
	
	return 0;
}

