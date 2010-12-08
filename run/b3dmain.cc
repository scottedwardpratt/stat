#include "b3d.h"
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
	CB3D *b3d=new CB3D(run_name);
	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.dat");
	neventsmax=parameter::getI(b3d->parmap,"B3D_NEVENTSMAX",1000);
	//b3d->randy->reset(-time(NULL));

	for(int iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		b3d->SetQualifier(qualifiers.qualifier[iqual]);
		qualifiers.SetPars(&(b3d->parmap),iqual);
		b3d->hydrotob3d->ReadInput();
		for(ievent=0;ievent<neventsmax;ievent++){
			nparts=b3d->hydrotob3d->MakeEvent();
			if(b3d->VIZWRITE){
				for(double tau=1.0;tau<20.01;tau+=0.5) b3d->AddAction_VizWrite(tau);
			}
			//printf("DeadActionMap size=%d, deadparts=%d, liveparts=%d\n",int(b3d->DeadActionMap.size()),int(b3d->DeadPartMap.size()),int(b3d->PartMap.size()));
			b3d->PerformAllActions();
			dnchdy=b3d->WriteDataH5();
			printf("##### finished event %d ##### dNch/dy=%g #####\n",b3d->ievent_write,double(dnchdy));
			//printf("nactivate=%lld, nexit/N=%g, nscatter/N=%g, nmerge/N=%g, ndecay/N=%g,\n", b3d->nactivate,double(b3d->nexit)/N,double(b3d->nscatter)/N,double(b3d->nmerge)/N,double(b3d->ndecay)/N);
			//printf("ncheck=%lld, nscatter=%lld, npass=%lld\n",b3d->ncheck,b3d->nscatter,b3d->npass);

		}
	}

	return 0;
}
