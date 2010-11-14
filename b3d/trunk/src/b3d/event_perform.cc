#ifndef __ACTION_PERFORM_CC__
#define __ACTION_PERFORM_CC__
#include "b3d.h"

void CAction::Perform(){
	if(currentmap!=&(b3d->ActionMap)){
		printf("FATAL: trying to perform dead action\n");
		exit(1);
	}

	DeleteFromCurrentMap();
	b3d->tau=tau;
	//if(type==1) printf("tau=%g, nactions=%lld, type=%d\n",tau,b3d->nactions,type);
	//printf("tau=%g, nactions=%lld, type=%d\n",tau,b3d->nactions,type);

	if(tau+1.0E-4<b3d->tau){
		printf("FATAL:: action earlier than tau!!!!, b3d->tau=%15.10e, action tau=%15.10e\n",b3d->tau,tau);
		exit(1);
	}
	//printf("performing action: type=%d\n",type);
	if(type==0) PerformActivate();
	else if(type==1) PerformDecay();
	else if(type==2) PerformCollide();
	else if(type==3) PerformVizWrite();
	//else if(type==3) PerformResetCollisions();
	//else if(type==4) PerformSwallowParticles();
	else if(type==6) PerformExitCell();
	else{
		printf("FATAL: action type = %d is unknown, exiting\n",type);
		exit(1);
	}
	//b3d->CheckParts();

	CleanPartMap();
	tau=0.0;
	AddToMap(&(b3d->DeadActionMap));
	currentmap=&(b3d->DeadActionMap);
	//printf("action finished\n");
}


#endif
