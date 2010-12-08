#include "b3d.h"
#include "analyze.h"
#include "qualifier.h"
using namespace std;

int main(int argc, char *argv[]){
	if (argc != 2) {
		printf("Usage: b3d run_name\n");
		exit(-1);
	}
	double dnchdy;
	int nparts;
	int ievent=0,iqual,neventsmax;
	string run_name=argv[1];
	CAnalyze *anal=new CAnalyze(run_name);
	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.dat");
	
	for(int iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		anal->SetQualifier(qualifiers.qualifier[iqual]);
		anal->CALCGARRAYS=false;
		anal->CalcSpectra();
	}
	//anal->CalcGamma();
	//anal->CalcBalance();
	//anal->ReadVizData(20.0);
	return 0;
}
