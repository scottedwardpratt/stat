
#ifndef __RHICSTATMCMC_CC__
#define __RHICSTATMCMC_CC__
#include "rhicstat.h"
using namespace std;

void CRHICStat::Metropolis(){
	CRunInfo *bestrun=new CRunInfo(NX,NY);
	int ix,ixx,iz,nwrite,ibin,NBINS=50;
	unsigned long long int isample,iburn,nsuccess=0,norm=0;
	CRandom *randy=new CRandom(-1234);
	double ll,oldll,*x,*oldx,*realx,*oldrealx,*bestx,*besty,*xmcmcbar,bestll=-1.0E99,root12=sqrt(12),*bestz;
	double xpmax=1.0; // 1.0 corresponds to xmin and xmax at boundaries
	bool success;
	double **spread;
	double **xdist=new double*[NX];
	for(ix=0;ix<NX;ix++){
		xdist[ix]=new double[NBINS]();
	}
	if(xmcmc==NULL){
		xmcmc=new double[NX];
		for(ix=0;ix<NX;ix++){
			xmcmc[ix]=xbar[ix];
		}
	}
	x=new double[NX];
	bestx=new double[NX];
	realx=new double[NX];
	oldx=new double[NX];
	oldrealx=new double[NX];
	besty=new double[NY];
	xmcmcbar=new double[NX]();
	spread=new double *[NX];
	for(ix=0;ix<NX;ix++) spread[ix]=new double[NX]();
		FILE *mcmc=fopen("mcmctrace.dat","w");
	ll=1.0E-99;
	for(ix=0;ix<NX;ix++){
		oldx[ix]=(xmcmc[ix]-xbar[ix])*root12/(xmax[ix]-xmin[ix]);
		oldrealx[ix]=xbar[ix]+oldx[ix]*(xmax[ix]-xmin[ix])/root12;
	}
	oldll=GetLL(oldx);
	iburn=0;
	while(iburn<NBURN){
		do{
			iburn+=1;
			success=false;
			for(ix=0;ix<NX;ix++){
				x[ix]=oldx[ix]+0.1*randy->gauss();
			}
			for(ix=0;ix<NX;ix++){
				realx[ix]=xbar[ix]+x[ix]*(xmax[ix]-xmin[ix])/root12;
				while(realx[ix]>xmax[ix]||realx[ix]<xmin[ix]){
					if(realx[ix]>xmax[ix]) realx[ix]=xmax[ix]-(realx[ix]-xmax[ix]);
					if(realx[ix]<xmin[ix]) realx[ix]=xmin[ix]+(xmin[ix]-realx[ix]);
					x[ix]=(realx[ix]-xbar[ix])*root12/(xmax[ix]-xmin[ix]);
				}
			}
			ll=GetLL(x);
				//printf("iburn=%d, ll=%g\n",iburn,oldll);
			if((ll>oldll || randy->ran()<exp(ll-oldll)) && (ll-oldll>-40)){
				success=true;
				oldll=ll;
				for(ix=0;ix<NX;ix++){
					oldx[ix]=x[ix];
					oldrealx[ix]=realx[ix];
				}
			}
		}while(success==false);
	}

	norm=0;
	nwrite=0;
	isample=0;
	while(isample<NMCMC){
		do{
			isample+=1;
			if((10*isample%NMCMC)==0){
				printf("----- finished %g percent of MCMC search -----\n",double(100*isample/NMCMC));
			}
			nwrite+=1;
			success=false;
			for(ix=0;ix<NX;ix++){
				x[ix]=oldx[ix]+0.1*randy->gauss();
			}
			for(ix=0;ix<NX;ix++){
				realx[ix]=xbar[ix]+x[ix]*(xmax[ix]-xmin[ix])/root12;
				while(realx[ix]>xmax[ix]||realx[ix]<xmin[ix]){
					if(realx[ix]>xmax[ix]) realx[ix]=xmax[ix]-(realx[ix]-xmax[ix]);
					if(realx[ix]<xmin[ix]) realx[ix]=xmin[ix]+(xmin[ix]-realx[ix]);
					x[ix]=(realx[ix]-xbar[ix])*root12/(xmax[ix]-xmin[ix]);
				}
			}
			ll=GetLL(x);
			if((ll>oldll || randy->ran()<exp(ll-oldll)) && (ll-oldll>-40)){
				success=true;
				nsuccess+=1;
				oldll=ll;
				for(ix=0;ix<NX;ix++){
					oldx[ix]=x[ix];
					oldrealx[ix]=realx[ix];
				}
				if(ll>bestll){
					bestll=ll;
					for(ix=0;ix<NX;ix++) bestx[ix]=x[ix];
						printf("bestll=%g\n",bestll);
				}
			}
			if(nwrite==10){
				for(ix=0;ix<NX;ix++){
					fprintf(mcmc,"%7.4f ",(oldrealx[ix]-xmin[ix])/(xmax[ix]-xmin[ix]));
				}
				fprintf(mcmc,"\n");
				nwrite=0;
				norm+=1;
				for(ix=0;ix<NX;ix++){
					xmcmcbar[ix]+=oldrealx[ix];
					for(ixx=0;ixx<NX;ixx++){
						spread[ix][ixx]+=oldx[ix]*oldx[ixx];
					}
				}
				for(ix=0;ix<NX;ix++){
					ibin=lrint(floor(NBINS*(oldrealx[ix]-xmin[ix])/(xmax[ix]-xmin[ix])));
					if(ibin<0 || ibin>=NBINS){
						printf("ibin=%d is out of range\n",ibin);
						exit(1);
					}
					xdist[ix][ibin]+=1.0;
				}
			}
		}while(success==false);
	}


	printf("----- nsuccess=%lld, bestll=%g\n",nsuccess,bestll);
	bestz=new double[NY];
	for(ix=0;ix<NX;ix++){
		bestrun->x[ix]=bestx[ix];
			//printf("%25s= %7.4f\n",xname[ix].c_str(),xbar[ix]+(xmax[ix]-xmin[ix])*bestx[ix]/root12);
	}
	zgetter->GetZ(bestx,bestz);
	for(iz=0;iz<NY;iz++){
		bestrun->z[iz]=bestz[iz];
	}
	GetYFromZ(bestrun);
	PrintX(bestrun);
	PrintY(bestrun);
	delete [] bestz;

	printf("  xbar=(");
	for(ix=0;ix<NX;ix++){
		xmcmcbar[ix]=xmcmcbar[ix]/double(norm);
		if(ix!=0) printf(",");
		printf("%7.4f",xmcmcbar[ix]);
	}
	printf(")\n");
	printf("SPREAD = \n");
	for(ix=0;ix<NX;ix++){
		for(ixx=0;ixx<NX;ixx++){
			spread[ix][ixx]=spread[ix][ixx]/(12.0*norm);
			spread[ix][ixx]-=xmcmcbar[ix]*xmcmcbar[ixx];
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
		xmcmc[ix]=oldrealx[ix];
	}
	delete [] oldx;
	delete [] oldrealx;
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
		zgetter->GetZ(x,z);
		double ll=0.0;
		for(iz=0;iz<NZ;iz++){
			ll-=0.5*pow(z[iz]-expinfo->z[iz],2);
		}
		delete [] z;
		ll=ll/(1.0+SIGMA2_EMULATOR);
		return ll;
	}

#endif