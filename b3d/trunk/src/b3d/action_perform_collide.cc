#ifndef __ACTION_PERFORM_COLLIDE_CC__
#define __ACTION_PERFORM_COLLIDE_CC__
#include "b3d.h"

void CAction::PerformCollide(){
	double scompare;
	double tautest;
	int ncolls,colltype;
	CPart *part1,*part2;
	CPartMap::iterator ppos1,ppos2,pbegin,pend,ppos;
	double sigma,einitial,efinal;

	if(partmap.size()!=2){
		printf("FATAL: wrong number of particles in partmap for collision, =%d\n",int(partmap.size()));
		exit(1);
	}
	ppos1=partmap.begin();
	ppos2=ppos1; ++ppos2;
	part1=ppos1->second;
	part2=ppos2->second;
	
	//einitial=part1->p[0]+part2->p[0];

	part1->Propagate(b3d->tau);
	/*
	if(part1->cell!=part1->FindCell() && b3d->tau<b3d->TAUCOLLMAX){
		printf("CAction::PerformCollide, part1: cells different\n");
		part1->cell->Print();
		(part1->FindCell())->Print();
		exit(1);
	}
	*/
	part2->Propagate(b3d->tau);
	/*
	if(part2->cell!=part2->FindCell() && b3d->tau<b3d->TAUCOLLMAX){
		printf("CAction::PerformCollide, part2: cells different\n");
		part2->cell->Print();
		(part2->FindCell())->Print();
		exit(1);
	}
	*/


	if(part1->active && part2->active){
		scompare=b3d->GetPiBsquared(part1,part2);
	/*
	if(scompare<1.5 && (fabs(part1->r[1]-part2->r[1])>b3d->DXY || fabs(part1->r[2]-part2->r[2])>b3d->DXY || fabs(part1->eta-part2->eta)>b3d->DETA)){
	printf("x1-x2=%g or y1-y2=%g > DXY or eta1-eta2=%g > DETA, scompare=%g\n",part1->r[1]-part2->r[1],part1->r[2]-part2->r[2],part1->eta-part2->eta,scompare);
	}
	*/
		//printf("scompare=%g, einitial=%g\n",scompare,part1->p[0]+part2->p[0]);
	colltype=b3d->Collide(part1,part2,scompare);

		//if(colltype==0||colltype==2) efinal=part1->p[0]+part2->p[0];
		//else efinal=part1->p[0];
		//if(fabs(einitial-efinal)>1.0){
			//printf("perform collide, colltype=%d: energy changed, einitial=%g, efinal=%g, colltype=%d\n",colltype,einitial,efinal,colltype);
			//exit(1);
		//}
  }
	else{
		colltype=-1;
		printf("Trying to collide with dead dude\n");
		exit(1);
	}


	if(colltype==0){
		b3d->nactions+=1;
		b3d->npass++;
	}

	if(colltype==1){
		b3d->nactions+=1;
		b3d->nmerge+=1;
		part1->tau_lastint=part1->tau0;
		part1->actionmother=b3d->nactions;
		part1->Reset();
		if(abs(part1->resinfo->code)==2224 || abs(part1->resinfo->code)==2214 || abs(part1->resinfo->code)==2114 || abs(part1->resinfo->code)==1114 || abs(part1->resinfo->code)==2111 || abs(part1->resinfo->code)==2212) b3d->ncheck+=1;
	}
	if(colltype==2){
		b3d->nscatter+=1;
		part1->tau_lastint=part1->tau0;
		part2->tau_lastint=part2->tau0;
		b3d->nactions+=1;
		part1->actionmother=b3d->nactions;
		part2->actionmother=b3d->nactions;
		part1->Reset();
		part2->Reset();
		if(abs(part1->resinfo->code)==2224 || abs(part1->resinfo->code)==2214 || abs(part1->resinfo->code)==2114 || abs(part1->resinfo->code)==1114 || abs(part1->resinfo->code)==2111 || abs(part1->resinfo->code)==2212) b3d->ncheck+=1;
		if(abs(part2->resinfo->code)==2224 || abs(part1->resinfo->code)==2214 || abs(part1->resinfo->code)==2114 || abs(part1->resinfo->code)==1114 || abs(part1->resinfo->code)==2111 || abs(part1->resinfo->code)==2212) b3d->ncheck+=1;
	}
	

}

#endif
