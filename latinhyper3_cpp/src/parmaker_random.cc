#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <cstring>
#include "coralutils.h"

using namespace std;

int main(int argc, char *argv[]){
	if(argc!=2){
		printf("USAGE: parmaker nruns\n");
		exit(1);
	}
	CRandom *randy=new CRandom(-1234);
	const int NPARSMAX=40;
	int npars=0,ipar;
	char dummy1[100],dummy2[100],dirname[100];
	char name[NPARSMAX][40],type[NPARSMAX][10];
	double MINPAR[NPARSMAX],MAXPAR[NPARSMAX];
	string statsfilename,command;
	double value,min,max;
	int nruns=atoi(argv[1]),irun;
	FILE *input,*output,*latin3;
	input=fopen("ranges.dat","r");
	while(!feof(input)){
		fscanf(input,"%s %s %lf %lf",type[npars],name[npars],&MINPAR[npars],&MAXPAR[npars]);
		npars+=1;
	}
	fclose(input);
	printf("npars=%d\n",npars);
	for(ipar=0;ipar<npars;ipar++)
		printf("%s %s: min=%g, max=%g\n",type[ipar],name[ipar],MINPAR[ipar],MAXPAR[ipar]);
	
		//sprintf(dummy1,"rm -f latin3.dat");
		//system(dummy1);

	for(irun=1;irun<=nruns;irun++){
		sprintf(dirname,"parameters/run%d",irun);
		//printf("making directory %s\n",dirname);
		command="mkdir -p "+string(dirname);
		system(command.c_str());
		command="cp parameters/default/fixed.param "+string(dirname)+"/";
		system(command.c_str());
		statsfilename=string(dirname)+"/stats.param";
		output=fopen(statsfilename.c_str(),"w");
		for(ipar=0;ipar<npars;ipar++){
			value=MINPAR[ipar]+randy->ran()*(MAXPAR[ipar]-MINPAR[ipar]);
			//printf("%s %s=value=%g; min=%g, max=%g\n",type[ipar],name[ipar],value,MINPAR[ipar],MAXPAR[ipar]);
			fprintf(output,"%s %s %g\n",type[ipar],name[ipar],value);
		}
		fclose(output);
	}
  return 0;
}


