#ifndef __COLLIDE_CC__
#define __COLLIDE_CC__
#include "b3d.h"
using namespace std;
#define __SCATTERING_ON__

int CB3D::Collide(CPart *part1,CPart *part2,double scompare){
	double sigma=0.0,Gamma,G2,MR,M,m1,m2,b,q2,qR2,tan2delta;
	int ir1,ir2,irflip,alpha;
	CMerge *merge;
	double jR,j1,j2;
	ir1=part1->resinfo->ires; ir2=part2->resinfo->ires;
	if(ir1>ir2){
		irflip=ir1; ir1=ir2; ir2=irflip;
	}
	sigma=SIGMADEFAULT/NSAMPLE; // Use 15 mb for s-wave scattering
	//sigma=0.0;
	if(sigma>scompare){
		//printf("COLLIDE, einitial=%g\n",part1->p[0]+part2->p[0]);
		Scatter(part1,part2);
		return 2;
	}
	else{
		merge=reslist->MergeArray[ir1][ir2];
		if(merge!=NULL){
			j1=part1->resinfo->spin;
			j2=part2->resinfo->spin;
			m1=part1->GetMass();
			m2=part2->GetMass();
			M=pow(part1->p[0]+part2->p[0],2);
			for(alpha=1;alpha<4;alpha++) M-=pow(part1->p[alpha]+part2->p[alpha],2);
			M=sqrt(M);
			if(M<m1+m2){
				printf("SCREWY MASSES: m1=%g, m2=%g, M=%g\n",m1,m2,M);
				part1->Print();
				part2->Print();
				//Misc::Pause();
				q2=1.0E-6;
			}
			else q2=Misc::triangle(M,m1,m2);
		}
		while(merge!=NULL){
			Gamma=merge->resinfo->width;
			G2=0.25*Gamma*Gamma;
			b=merge->branching;
			jR=merge->resinfo->spin;
			MR=merge->resinfo->mass;
			if(m1+m2<MR){
				qR2=Misc::triangle(MR,m1,m2);
				tan2delta=(b*G2/((M-MR)*(M-MR)))*2.0*q2/(qR2+q2);
				sigma+=((4.0*PI*HBARC*HBARC/q2)*(tan2delta/(1.0+tan2delta))
					*((2.0*jR+1.0)/((2.0*j1+1.0)*(2.0*j2+1.0))))/NSAMPLE;
				if(sigma>scompare){
					Merge(part1,part2,merge->resinfo);
					//printf("Merging: sigma=%g, m1=%g, m2=%g, mR=%g\n",sigma,m1,m2,MR);
					return 1;
				}
			}
			merge=merge->next;
		}
	}
	if(part1->r[0]!=part1->r[0] || part2->r[0]!=part2->r[0]){
		printf("NaN problem in CB3D::Collide");
		exit(1);
	}
	return 0;
}

void CB3D::Scatter(CPart *part1,CPart *part2){
	double ptot[4],u[4],q[4],qprime[4],roots=0.0,newroots,taucoll,g[4]={1.0,-1.0,-1.0,-1.0};
	double ctheta,stheta,phi,qmag,vr1,vr2,m1,m2,*p1=part1->p,*p2=part2->p;
	double y1,mt;
	int nparts,ipart,alpha;
	bool flip=false;
	
#ifdef __SCATTERING_ON__
	m1=part1->GetMass();
	m2=part2->GetMass();
	CCell *cell1=part1->cell,*cell2=part2->cell;
	if(BJORKEN && ((cell1->ieta==0 && cell2->ieta==2*NETA-1) || (cell1->ieta==2*NETA-1 && cell2->ieta==0))){
		flip=true;
		y1=part1->y;
		if(cell1->ieta==0) y1+=2.0*ETAMAX;
		else y1-=2.0*ETAMAX;
		mt=sqrt(p1[0]*p1[0]-p1[3]*p1[3]);
		p1[3]=mt*sinh(y1);
		p1[0]=mt*cosh(y1);
	}
	//else printf("b=%g, delrperp=%g, deleta=%g\n",sqrt(GetPiBsquared(part1,part2))/PI,
	//sqrt(pow(part1->r[1]-part2->r[1],2)+pow(part1->r[2]-part2->r[2],2)),fabs(part1->eta-part2->eta));

	for(alpha=0;alpha<4;alpha++){
		ptot[alpha]=p1[alpha]+p2[alpha];
		q[alpha]=0.5*(p1[alpha]-p2[alpha]);
		roots+=g[alpha]*ptot[alpha]*ptot[alpha];
	}
	roots=sqrt(roots);
	//printf("BEFORE: roots=%g, ptot=(%g,%g,%g,%g)\n",roots,ptot[0],ptot[1],ptot[2],ptot[3]);
	for(alpha=0;alpha<4;alpha++) u[alpha]=ptot[alpha]/roots;
	Misc::BoostToCM(u,q,qprime);
	qmag=sqrt(qprime[1]*qprime[1]+qprime[2]*qprime[2]+qprime[3]*qprime[3]);
	
	ctheta=1.0-2.0*randy->ran();
	phi=2.0*PI*randy->ran();
	qprime[3]=qmag*ctheta;
	stheta=sqrt(1.0-ctheta*ctheta);
	qprime[1]=qmag*stheta*cos(phi);
	qprime[2]=qmag*stheta*sin(phi);
	//printf("qprime[0]=%g =? %g\n",qprime[0],0.5*sqrt(m1*m1+qmag*qmag)-0.5*sqrt(m2*m2+qmag*qmag));

	qprime[0]=sqrt(m1*m1+qmag*qmag);
	Misc::Boost(u,qprime,p1);
	part1->SetY();
	for(alpha=1;alpha<4;alpha++) qprime[alpha]=-qprime[alpha];
	qprime[0]=sqrt(m2*m2+qmag*qmag);
	Misc::Boost(u,qprime,p2);
	part2->SetY();
	
  /*
	newroots=0.0;
	for(alpha=0;alpha<4;alpha++){
		ptot[alpha]=p1[alpha]+p2[alpha];
		newroots+=g[alpha]*ptot[alpha]*ptot[alpha];
	}
	newroots=sqrt(newroots);
	//printf("AFTER:  roots=%g, ptot=(%g,%g,%g,%g)\n",newroots,ptot[0],ptot[1],ptot[2],ptot[3]);
	if(fabs(roots-newroots)>1.0){
		printf("in CB3d::Collide oldroots=%g, newroots=%g\n",roots,newroots);
		part1->Print();
		part2->Print();
		exit(1);
	}
	*/
	
	if(flip){
		y1=part1->y;
		if(cell1->ieta==0){
			y1-=2.0*ETAMAX;
		}
		else{
			y1+=2.0*ETAMAX;
		}
		mt=sqrt(p1[0]*p1[0]-p1[3]*p1[3]);
		p1[3]=mt*sinh(y1);
		p1[0]=mt*cosh(y1);
		part1->y=y1;
	}
#endif
	//part1->Setp0();
	//part2->Setp0();
	//if(part1->cell!=part2->cell) ncheck+=1;
}

void CB3D::Merge(CPart *part1,CPart *part2,CResInfo *resinfo){
	int alpha;
	double ptot[4],minv,mt,y1,y2;
	bool flip=false;
	CPart *pswitch;
	CCell *cell1=part1->cell,*cell2=part2->cell;
	double *p1=part1->p,*p2=part2->p;

	if(BJORKEN && ((cell1->ieta==0 && cell2->ieta==2*NETA-1) || (cell1->ieta==2*NETA-1 && cell2->ieta==0))){
		flip=true;
		p2=new double[4];
		y2=part2->y;
		mt=sqrt(p2[0]*p2[0]-p2[3]*p2[3]);
		for(alpha=0;alpha<4;alpha++) p2[alpha]=part2->p[alpha];
		if(cell1->ieta==0){
			y2-=2.0*ETAMAX;
		}
		else{
			y2+=2.0*ETAMAX;
		}
		p2[3]=mt*sinh(part2->y);
		p2[0]=mt*cosh(part2->y);
	}
	for(alpha=0;alpha<4;alpha++) ptot[alpha]=part1->p[alpha]+part2->p[alpha];
	minv=sqrt(ptot[0]*ptot[0]-ptot[1]*ptot[1]-ptot[2]*ptot[2]-ptot[3]*ptot[3]);
	if(minv<resinfo->minmass){
		//printf("CB3D::Merge cancel merge for resonance with ID=%d and nominal mass of %g. No sense with masses, minv=%g, minmass=%g\n",resinfo->code,resinfo->mass,minv,resinfo->minmass);
		part1->actionmother=part2->actionmother=nactions;
	}
	else{
		for(alpha=0;alpha<4;alpha++) part1->p[alpha]=part1->p[alpha]+part2->p[alpha];
		part1->resinfo=resinfo;
		part1->SetY();
		part2->Kill();
		part1->Reset();
	}
	if(flip){
		delete [] p2; 
	}
}

#endif
