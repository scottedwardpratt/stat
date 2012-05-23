
#ifndef __RHICSTATMCMC_CC__
#define __RHICSTATMCMC_CC__
#include "coral.h"
#include "rhicstat.h"
using namespace std;

void CRHICStat::Metropolis(unsigned int nburn,unsigned int nsample){
	int ix,ixx,iz,nwrite,ibin,NBINS=50;
	unsigned int isample,iburn,nsuccess=0,norm=0;
	CRandom *randy=new CRandom(-1234);
	double ll,oldll,*x,*oldx,*bestx,*xbarbar,bestll=-1.0E99,root12=sqrt(12);
	double xpmax=sqrt(3); // sqrt(3) corresponds to xmin and xmax as boundaries
	bool success;
	double **spread;
	double **xdist=new double*[NX];
	for(ix=0;ix<NX;ix++) xdist[ix]=new double[NBINS]();
	oldx=new double[NX];
	x=new double[NX];
	bestx=new double[NX];
	xbarbar=new double[NX]();
	spread=new double *[NX];
	for(ix=0;ix<NX;ix++) spread[ix]=new double[NX]();
	FILE *mcmc=fopen("mcmctrace.dat","w");
	ll=1.0E-99;
	for(ix=0;ix<NX;ix++) oldx[ix]=xmin[ix]+randy->ran()*(xmax[ix]-xmin[ix]);
	oldll=GetLL(oldx);
	iburn=0;
	while(iburn<nburn){
		do{
			iburn+=1;
			success=false;
			for(ix=0;ix<NX;ix++) x[ix]=oldx[ix]+0.1*randy->gauss();
			for(ix=0;ix<NX;ix++) if(fabs(x[ix])>xpmax) x[ix]-=2.0*xpmax*x[ix]/fabs(x[ix]);
			ll=GetLL(x);
				//printf("iburn=%d, ll=%g\n",iburn,oldll);
			if((ll>oldll || randy->ran()<exp(ll-oldll)) && (ll-oldll>-40)){
				success=true;
				oldll=ll;
				for(ix=0;ix<NX;ix++) oldx[ix]=x[ix];
			}
		}while(success==false);
	}
	
	norm=0;
	nwrite=0;
	isample=0;
	while(isample<nsample){
		do{
			isample+=1;
			nwrite+=1;
			success=false;
			for(ix=0;ix<NX;ix++) x[ix]=oldx[ix]+0.1*randy->gauss();
			for(ix=0;ix<NX;ix++){
				while(fabs(x[ix])>xpmax){
					x[ix]-=2.0*xpmax*x[ix]/fabs(x[ix]);
				}
			}
			ll=GetLL(x);
			if((ll>oldll || randy->ran()<exp(ll-oldll)) && (ll-oldll>-40)){
				success=true;
				nsuccess+=1;
				oldll=ll;
				for(ix=0;ix<NX;ix++) oldx[ix]=x[ix];
				if(ll>bestll){
					bestll=ll;
					for(ix=0;ix<NX;ix++) bestx[ix]=x[ix];
					printf("bestll=%g\n",bestll);
				}
			}
			if(nwrite==100){
				for(ix=0;ix<NX;ix++) fprintf(mcmc,"%7.4f ",oldx[ix]/xpmax);
				fprintf(mcmc,"\n");
				nwrite=0;
				norm+=1;
				for(ix=0;ix<NX;ix++){
					xbarbar[ix]+=oldx[ix];
					for(ixx=0;ixx<NX;ixx++){
						spread[ix][ixx]+=oldx[ix]*oldx[ixx];
					}
					norm+=1;
				}
				for(ix=0;ix<NX;ix++){
					ibin=lrint(floor((1.0+(oldx[ix]/xpmax))*0.5*NBINS));
					if(ibin<0 || ibin>=NBINS){
						printf("ibin=%d is out of range\n",ibin);
						exit(1);
					}
					xdist[ix][ibin]+=1.0;
				}
			}
		}while(success==false);
	}
	
	
	printf("----- nsuccess=%d, bestll=%g, bestx'=\n",nsuccess,bestll);
	
	for(ix=0;ix<NX;ix++){
		bestx[ix]=xbar[ix]+(xmax[ix]-xmin[ix])*bestx[ix]/root12;
		printf("%25s= %7.4f\n",xname[ix].c_str(),bestx[ix]);
	}

	printf("  xbar=(",nsuccess);
	for(ix=0;ix<NX;ix++){
		xbarbar[ix]=xbarbar[ix]/(root12*double(norm));
		if(ix!=0) printf(",");
		printf("%7.4f",xbarbar[ix]);
	}
	printf(")\n");
	printf("SPREAD = \n");
	for(ix=0;ix<NX;ix++){
		for(ixx=0;ixx<NX;ixx++){
			spread[ix][ixx]=spread[ix][ixx]/(12.0*norm);
			spread[ix][ixx]-=xbarbar[ix]*xbarbar[ixx];
			if(ixx==0) printf(" ");
			printf("%7.4f",12.0*spread[ix][ixx]);
		}
		printf("\n");
	}
	fclose(mcmc);
	
	FILE *xdfile;
	char filename[80];
	for(ix=0;ix<NX;ix++){
		sprintf(filename,"xdist%d.dat",ix);
		xdfile=fopen(filename,"w");
		fprintf(xdfile,"%s\n",xname[ix].c_str());
		for(ibin=0;ibin<NBINS;ibin++){
			fprintf(xdfile,"%3d %g\n",ibin,NBINS*xdist[ix][ibin]/double(norm));
		}
		fclose(xdfile);
	}
	delete [] oldx;
	delete [] x;
	for(ix=0;ix<NX;ix++){
		delete [] spread[ix];
		delete [] xdist[ix];
	}
	delete [] spread;
	delete [] xdist;
	
}

double CRHICStat::GetLL(double *x){
	int ix,iz,NZ=NX;
	double *z=new double[NZ];
	GetZquad(x,z);
	double ll=0.0;
	for(iz=0;iz<NZ;iz++){
		ll-=0.5*pow(z[iz]-expinfo->z[iz],2);
	}
	delete [] z;
	return ll;
}

#endif