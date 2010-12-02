#include "b3d.h"
using namespace std;

int main(int argc, char *argv[]){
	if (argc != 2) {
		printf("Usage: b3d run_name\n");
		exit(-1);
	}
	char qchar[120],rchar[120];
	const int bListSize = 5;
	double bList[bListSize] = {0.35, 2.2, 3.7, 5.2, 7.0};
	double dnchdy,N,nparts;
	bool filealive;
	int ievent=0,ib,neventsmax;
	string run_name=argv[1],qualifier;
	CB3D *b3d=new CB3D(run_name);
	neventsmax=parameter::getI(b3d->parmap,"B3D_NEVENTSMAX",1);
	b3d->randy->reset(-time(NULL));

	for(ib=0;ib<bListSize;ib++){
		sprintf(qchar,"%g",bList[ib]);
		qualifier="b"+string(qchar);
		b3d->SetQualifier(qualifier);
		printf("check a\n");
		b3d->hydrotob3d->ReadInput();
		printf("check b\n");
		for(ievent=0;ievent<neventsmax;ievent++){
		//printf("eventlist size=%d, deadparts=%d, liveparts=%d\n",int(b3d->DeadEventMap.size()),int(b3d->DeadPartMap.size()),int(b3d->PartMap.size()));
			printf("check c, ievent=%d\n",ievent);	
			nparts=b3d->hydrotob3d->MakeEvent();
			printf("nparts=%d\n",nparts);
			b3d->PrintPartList();
			printf("check cc\n");
			if(b3d->VIZWRITE){
				for(double tau=1.0;tau<20.01;tau+=0.5) b3d->AddAction_VizWrite(tau);
			}
			printf("check ccc\n");
			b3d->PerformAllActions();
			printf("check d\n");
			dnchdy=b3d->WriteDataH5();
			printf("check e\n");
			N=double(b3d->nactivate);
			printf("##### finished event %d ##### dNch/dy=%g #####\n",b3d->ievent_write,double(dnchdy));
			printf("nactivate=%lld, nexit/N=%g, nscatter/N=%g, nmerge/N=%g, ndecay/N=%g,\n", b3d->nactivate,double(b3d->nexit)/N,double(b3d->nscatter)/N,double(b3d->nmerge)/N,double(b3d->ndecay)/N);
			printf("ncheck=%lld, nscatter=%lld, npass=%lld\n",b3d->ncheck,b3d->nscatter,b3d->npass);

		}
	}

	return 0;
}
