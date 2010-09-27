#ifndef __DECAY_CC__
#define __DECAY_CC__
#include "b3d.h"

using namespace std;

void CB3D::Decay(CPart *&mother,int &nbodies, CPart **&daughter){
	const double HBARC=197.326;
	int ibody,alpha;
	double u[4],mass[4],mtot;
	CPart *dptr;

	double *p[4],kprime[4];
	double e1,e2;
	double q,weight,wmax,sthet,cthet,phi;
	double p3mag,kprimemax,p3max,kprimemax2,kprimemag2,e1prime,e2prime;
	double e12,u12[4];
	//double einitial,efinal;

	mass[0]=mother->GetMass();
 //mass[0]=mother->resinfo->mass;
	p[0]=mother->p;
	//einitial=p[0][0];

	//motherresinfoptr->DecayGetResInfoptr(nbodies,daughterresinfoptr);

	/* Create daughter objects */
	mtot=0.0;
	for(ibody=0;ibody<nbodies;ibody++){
		mass[ibody+1]=daughter[ibody]->resinfo->mass;
		mtot+=mass[ibody+1];
		p[ibody+1]=daughter[ibody]->p;
	}
	if(mtot>mass[0]){
		printf("CB3D::Decay, This mass can't decay, mothermass=%g, mtot of products=%g\n",mass[0],mtot);
		exit(1);
	}

	/* TWO-BODY DECAYS */
	if(nbodies==2){
		cthet=1.0-2.0*randy->ran();
		sthet=sqrt(1.0-cthet*cthet);
		phi=2.0*PI*randy->ran();
		q=sqrt(Misc::triangle(mass[0],mass[1],mass[2]));
		p[1][3]=q*cthet;
		p[1][1]=q*sthet*cos(phi);
		p[1][2]=q*sthet*sin(phi);
		p[2][3]=-p[1][3];
		p[2][2]=-p[1][2];
		p[2][1]=-p[1][1];
		p[1][0]=sqrt(mass[1]*mass[1]+p[1][1]*p[1][1]+p[1][2]*p[1][2]+p[1][3]*p[1][3]);
		p[2][0]=sqrt(mass[2]*mass[2]+p[2][1]*p[2][1]+p[2][2]*p[2][2]+p[2][3]*p[2][3]);
	}
	/* THREE-BODY DECAYS */
	else if(nbodies==3){
		kprimemax2=Misc::triangle(mass[0]-mass[3],mass[1],mass[2]);
		kprimemax=sqrt(kprimemax2);
		p3max=sqrt(Misc::triangle(mass[0],mass[1]+mass[2],mass[3]));
		e1=sqrt(pow(mass[1],2)+p3max*p3max);
		e2=sqrt(pow(mass[2],2)+p3max*p3max);
		wmax=p3max*(e1*e2/(mass[1]*mass[2]))*(mass[1]+mass[2])/(e1+e2);
		do{
			TRY_AGAIN:
			do{
				kprime[1]=kprimemax*(2.0*randy->ran()-1.0);
				kprime[2]=kprimemax*(2.0*randy->ran()-1.0);
				kprime[3]=kprimemax*(2.0*randy->ran()-1.0);
				kprimemag2=kprime[1]*kprime[1]+
					kprime[2]*kprime[2]+kprime[3]*kprime[3];
			} while(kprimemag2>kprimemax2);
			e1prime=sqrt(kprimemag2+mass[1]*mass[1]);
			e2prime=sqrt(kprimemag2+mass[2]*mass[2]);
			if(e1prime+e2prime+mass[3]>mass[0]) goto TRY_AGAIN;
			p3mag=sqrt(Misc::triangle(mass[0],e1prime+e2prime,mass[3]));
			cthet=1.0-2.0*randy->ran();
			sthet=sqrt(1.0-cthet*cthet);
			phi=2.0*PI*randy->ran();
			p[3][3]=p3mag*cthet;
			p[3][1]=p3mag*sthet*cos(phi);
			p[3][2]=p3mag*sthet*sin(phi);
			p[3][0]=sqrt(p3mag*p3mag+mass[3]*mass[3]);
			e12=sqrt(pow(e1prime+e2prime,2)+p3mag*p3mag);
			for(alpha=1;alpha<4;alpha++) u12[alpha]=-p[3][alpha]/(e1prime+e2prime);
			u12[0]=sqrt(1.0+u12[1]*u12[1]+u12[2]*u12[2]+u12[3]*u12[3]);
			kprime[0]=e1prime;
			Misc::lorentz(u12,kprime,p[1]);
			kprime[0]=e2prime;
			for(alpha=1;alpha<=3;alpha++) kprime[alpha]=-kprime[alpha];
			Misc::lorentz(u12,kprime,p[2]);
			weight=p3mag*(p[1][0]*p[2][0]/(e1prime*e2prime))
				*((e1prime+e2prime)/(p[1][0]+p[2][0]));
		} while(randy->ran()<weight/wmax);
	}

	/* Boost the new particles */
	for(alpha=0;alpha<4;alpha++) u[alpha]=mother->p[alpha]/mother->GetMass();
	double pprime[4];
	for(ibody=0;ibody<nbodies;ibody++){
		dptr=daughter[ibody];
		Misc::lorentz(u,p[ibody+1],pprime);
		for(alpha=0;alpha<4;alpha++) dptr->p[alpha]=pprime[alpha];
		dptr->SetY();
		dptr->tau0=mother->tau0;
		for(alpha=0;alpha<4;alpha++) dptr->r[alpha]=mother->r[alpha];
		dptr->eta=mother->eta;
		dptr->Setp0();
		dptr->tau_lastint=daughter[ibody]->tau0;
		dptr->tau0=mother->tau0;
	}
}

#endif
