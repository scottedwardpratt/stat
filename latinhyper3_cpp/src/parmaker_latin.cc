#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <cstring>
#include "latin_center.h"

const double PI=3.14159265358979323844;
const double HBARC=197.3269602;

using namespace std;

int main(int argc, char *argv[]){
	if(argc!=2){
		printf("USAGE: parmaker nruns\n");
		exit(1);
	}
	const int NPARSMAX=100;
	int npars=0,ipar;
	char dummy1[100],dummy2[100],dirname[100];
	char name[NPARSMAX][40],type[NPARSMAX][10];
	double MINPAR[NPARSMAX],MAXPAR[NPARSMAX];
	string statsfilename,command;
	double value,min,max,**x;
	int nruns=atoi(argv[1]),irun,seed=1234;
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

	x=new double *[nruns];
	for(irun=0;irun<nruns;irun++)
		x[irun]=new double[npars];
	double *xx=new double[npars*nruns];
	latin_center(npars,nruns,&seed,xx);
	
	for(irun=0;irun<nruns;irun++){
		for(ipar=0;ipar<npars;ipar++){
			x[irun][ipar]=xx[ipar*nruns+irun];
			if(ipar*nruns+irun>=npars*nruns){
				printf("out of bounds\n");
			}
		}
		sprintf(dirname,"parameters/run%d",irun+1);
		//printf("making directory %s\n",dirname);
		command="mkdir -p "+string(dirname);
		system(command.c_str());
		command="cp parameters/default/fixed.param "+string(dirname)+"/";
		system(command.c_str());
		statsfilename=string(dirname)+"/stats.param";
		output=fopen(statsfilename.c_str(),"w");
		for(ipar=0;ipar<npars;ipar++){
			value=MINPAR[ipar]+x[irun][ipar]*(MAXPAR[ipar]-MINPAR[ipar]);
			fprintf(output,"%s %s %g\n",type[ipar],name[ipar],value);
			printf("%7.5f ",x[irun][ipar]);
		}
		printf("\n");
		fclose(output);
	}
	for(irun=0;irun<nruns;irun++)
		delete [] x[irun];
	delete [] x;
	delete [] xx;
	return 0;
}


