#ifndef __ACTION_PERFORM_COLLIDE_CC__
#define __ACTION_PERFORM_COLLIDE_CC__
#include "b3d.h"

void CAction::PerformExitCell(){
	CPart *part;
	CPartMap::iterator ppos;
	double mt;

	if(partmap.size()!=1){
		printf("FATAL: wrong number of particles in partmap for ExitCell, =%d\n",int(partmap.size()));
		exit(1);
	}
	ppos=partmap.begin();
	part=ppos->second;
	part->Propagate(tau);
	//printf("eta=%g, x=%g, y=%g\n",part->eta,part->r[1],part->r[2]);
	double *r=part->r;
	double *p=part->p;
	
	/*
	if(part->nextcell==NULL)
		printf("EXITING (part %d): tau=%g, old cell: (%d,%d,%d), new cell = NULL\n",part->listid,tau, part->cell->ix,part->cell->iy,part->cell->ieta);
	else
		printf("EXITING (part %d): tau=%g, old cell: (%d,%d,%d), new cell: (%d,%d,%d)\n",part->listid,tau, part->cell->ix,part->cell->iy,part->cell->ieta, part->nextcell->ix,part->nextcell->iy,part->nextcell->ieta);*/

	if(b3d->BJORKEN && part->cell->creflection!=NULL && part->nextcell==part->cell->creflection && fabs(fabs(part->eta)-b3d->ETAMAX)<1.0E-6){
		if(part->y<-b3d->ETAMAX){
			part->y+=2.0*b3d->ETAMAX;
			part->eta=b3d->ETAMAX;
		}
		else if(part->y>b3d->ETAMAX){
			part->y-=2.0*b3d->ETAMAX;
			part->eta=-b3d->ETAMAX;
		}
		else{
			printf("trying to reflect particle not moving fast enough!");
			Print();
			exit(1);
		}
		r[0]=tau*cosh(part->eta);
		r[3]=tau*sinh(part->eta);
		mt=sqrt(p[0]*p[0]-p[3]*p[3]);
		p[3]=mt*sinh(part->y);
		p[0]=mt*cosh(part->y);
		//printf("CHECK: m2=%g == %g\n",p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3],part->resinfo->mass*part->resinfo->mass);
		//part->Print();
		//part->cell->Print();
	}
	part->cell=part->nextcell;
	part->DeleteFromCurrentMap();
	part->KillActions();
	if(part->cell!=NULL){
		part->AddToMap(&(part->cell->partmap));
	//part->FindCollisions();
	}
	else{
		part->AddToMap(&(b3d->FinalPartMap));
		//printf("eta=%g, x=%g, y=%g\n",part->eta,part->r[1],part->r[2]);
	}
	part->Reset();
	b3d->nexit+=1;
	b3d->nactions+=1;
}

#endif
