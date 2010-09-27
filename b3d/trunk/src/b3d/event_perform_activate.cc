#ifndef __ACTION_PERFORM_ACTIVATE_CC__
#define __ACTION_PERFORM_ACTIVATE_CC__
#include "b3d.h"

void CAction::PerformActivate(){
	CPart *part;
	double *r,*p;
	CPartMap::iterator ppos,pend;
	ppos=partmap.begin();
	part=ppos->second;
	part->Propagate(tau);
	double mt;
	if(part->currentmap!=&(b3d->PartMap)){
		printf("FATAL: particles to be activated should be in PartMap\n");
		part->Print();
		exit(1);
	}
	part->cell=part->FindCell();
	p=part->p; r=part->r;
	part->DeleteFromCurrentMap();
	part->active=true;
	part->tau0=tau;
	if(b3d->BJORKEN && fabs(part->eta)>b3d->ETAMAX) part->CyclicReset();
	if(part->cell!=NULL) part->AddToMap(&(part->cell->partmap));
	else{
		part->AddToMap(&b3d->FinalPartMap);
	}
	part->tau_lastint=tau;
	part->actionmother=b3d->nactions;
	b3d->nactions++;
	
	part->cell=part->FindCell();
	//if(part->cell==NULL){
		//printf("part born in NULL space\n");
		//part->Print();
		//Misc::Pause();
	//}
	part->Reset();
	b3d->nactivate+=1;
	//ppos=part->GetPos(&(bb->PartMap));
	//pend=ppos; ++pend;
	//bb->RedoCollisions(ppos,pend);
}
#endif
