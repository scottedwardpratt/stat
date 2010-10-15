#include "b3d.h"
using namespace std;

int main(int argc, char *argv[]){
	if (argc != 3) {
		printf("Usage: b3d qualifier nevents\n");
		exit(-1);
	}
	char qchar[120];
	double dnchdy,N,nparts;
	bool filealive;
	int ievent=0;
	string qualifier=argv[1];
	int neventsmax=atoi(argv[2]);
	CB3D *b3d=new CB3D("parameters",qualifier);
	//b3d->randy->reset(-time(NULL));
	b3d->randy->reset(-123);
	b3d->bjmaker.Init();
	for(ievent=0;ievent<neventsmax;ievent++){
		//b3d->KillAllParts();
		//printf("eventlist size=%d, deadparts=%d, liveparts=%d\n",int(b3d->DeadEventMap.size()),int(b3d->DeadPartMap.size()),int(b3d->PartMap.size()));
		nparts=b3d->bjmaker.MakeEvent();
		b3d->PerformAllActions();
		dnchdy=b3d->WriteDataH5();
		N=double(b3d->nactivate);
		printf("##### finished event %d ##### dNch/dy=%g #####\n",b3d->ievent_write,double(dnchdy));
		printf("nactivate=%lld, nexit/N=%g, nscatter/N=%g, nmerge/N=%g, ndecay/N=%g,\n", b3d->nactivate,double(b3d->nexit)/N,double(b3d->nscatter)/N,double(b3d->nmerge)/N,double(b3d->ndecay)/N);
		printf("ncheck=%lld, nscatter=%lld, npass=%lld\n",b3d->ncheck,b3d->nscatter,b3d->npass);

	}

	return 0;
}
