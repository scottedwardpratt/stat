#ifndef __QUADFIT_CC__
#define __QUADFIT_CC__
#include "coral.h"
#include "rhicstat.h"
using namespace std;

void CRHICStat::QuadFit(){
	int i,j,i1,i2,j1,j2,iz,irun,n=(NX+1)*(NX+2)/2;
	CRunInfo *ri;
	double **M=new double *[n];
	double *error=new double[NRUNS];
	for(i=0;i<n;i++){
		M[i]=new double[n];
	}
	double *yy=new double[n];
	double *xx=new double[n];
	if(Aquad==NULL){
		Aquad=new double **[NX];
		for(iz=0;iz<NX;iz++){
			Aquad[iz]=new double*[NX];
			for(i1=0;i1<NX;i1++){
				Aquad[iz][i1]=new double[NX];
				for(i2=0;i2<NX;i2++) Aquad[iz][i1][i2]=0.0;
			}
		}
	}
	if(Bquad==NULL){
		Bquad=new double*[NX];
		for(iz=0;iz<NX;iz++){
			Bquad[iz]=new double[NX];
			for(i1=0;i1<NX;i1++) Bquad[iz][i1]=0.0;
		}
	}
	if(Cquad==NULL){
		Cquad=new double[NX];
		for(iz=0;iz<NX;iz++) Cquad[iz]=0.0;
	}
	CGSLMatrix_Real *gslmatrix=new CGSLMatrix_Real(n);
	
	/** for checking
	for(irun=0;irun<NRUNS;irun++){
		ri=runinfo[irun];
		for(iz=0;iz<NX;iz++){
			ri->z[iz]=-1.0+5*iz;
			for(i1=0;i1<NX;i1++){
				for(i2=0;i2<=i1;i2++) ri->z[iz]+=0.5*(iz+1+i1-0.12345*i2)*ri->x[i1]*ri->x[i2];
				ri->z[iz]+=0.6*(2.3352-iz)*(i1-0.25252*i2)*ri->x[i1];
			}
		}
	}
	*/
	
	for(iz=0;iz<NX;iz++){
		for(i=0;i<n;i++){
			yy[i]=xx[i]=0.0;
		}
		for(i=0;i<n;i++){
			for(j=0;j<n;j++) M[i][j]=0.0;
		}
			// from dM/dA_ij=0
		i=0;
		for(i1=0;i1<NX;i1++){
			for(i2=0;i2<=i1;i2++){
				j=0;
				for(j1=0;j1<NX;j1++){
					for(j2=0;j2<=j1;j2++){
						for(irun=0;irun<NRUNS;irun++){
							if(runinfo[irun]->good) M[i][j]+=runinfo[irun]->x[i1]*runinfo[irun]->x[i2]*runinfo[irun]->x[j1]*runinfo[irun]->x[j2];
						}
						j+=1;
					}
				}
				j=NX*(NX+1)/2;
				for(j1=0;j1<NX;j1++){
					for(irun=0;irun<NRUNS;irun++){
						if(runinfo[irun]->good) M[i][j]+=runinfo[irun]->x[i1]*runinfo[irun]->x[i2]*runinfo[irun]->x[j1];
					}
					j+=1;
				}
				j=-1+(NX+1)*(NX+2)/2;
				for(irun=0;irun<NRUNS;irun++){
					if(runinfo[irun]->good) M[i][j]+=runinfo[irun]->x[i1]*runinfo[irun]->x[i2];
				}
				yy[i]=0.0;
				for(irun=0;irun<NRUNS;irun++){
					if(runinfo[irun]->good) yy[i]+=runinfo[irun]->x[i1]*runinfo[irun]->x[i2]*runinfo[irun]->z[iz];
				}
				i+=1;
			}
		}
			// from dM/dB_i=0
		i=NX*(NX+1)/2;
		for(i1=0;i1<NX;i1++){
			j=0;
			for(j1=0;j1<NX;j1++){
				for(j2=0;j2<=j1;j2++){
					for(irun=0;irun<NRUNS;irun++){
						if(runinfo[irun]->good) M[i][j]+=runinfo[irun]->x[i1]*runinfo[irun]->x[j1]*runinfo[irun]->x[j2];
					}
					j+=1;
				}
			}
			j=NX*(NX+1)/2;
			for(j1=0;j1<NX;j1++){
				for(irun=0;irun<NRUNS;irun++){
					if(runinfo[irun]->good) M[i][j]+=runinfo[irun]->x[i1]*runinfo[irun]->x[j1];
				}
				j+=1;
			}
			j=-1+(NX+1)*(NX+2)/2;
			for(irun=0;irun<NRUNS;irun++){
				if(runinfo[irun]->good) M[i][j]+=runinfo[irun]->x[i1];
			}
			yy[i]=0.0;
			for(irun=0;irun<NRUNS;irun++){
				if(runinfo[irun]->good) yy[i]+=runinfo[irun]->x[i1]*runinfo[irun]->z[iz];
			}
			i+=1;
		}
			// from dM/dC=0
		i=-1+(NX+1)*(NX+2)/2;
		j=0;
		for(j1=0;j1<NX;j1++){
			for(j2=0;j2<=j1;j2++){
				for(irun=0;irun<NRUNS;irun++){
					if(runinfo[irun]->good) M[i][j]+=runinfo[irun]->x[j1]*runinfo[irun]->x[j2];
				}
				j+=1;
			}
		}
		j=NX*(NX+1)/2;
		for(j1=0;j1<NX;j1++){
			for(irun=0;irun<NRUNS;irun++){
				if(runinfo[irun]->good) M[i][j]+=runinfo[irun]->x[j1];
			}
			j+=1;
		}
		j=-1+(NX+1)*(NX+2)/2;
		for(irun=0;irun<NRUNS;irun++){
			if(runinfo[irun]->good) M[i][j]+=1.0;
		}
		yy[i]=0.0;
		for(irun=0;irun<NRUNS;irun++){
			if(runinfo[irun]->good) yy[i]+=runinfo[irun]->z[iz];
		}
		
		gslmatrix->SolveLinearEqs(yy,M,xx);
		
		i=0;
		for(i1=0;i1<NX;i1++){
			for(i2=0;i2<=i1;i2++){
				Aquad[iz][i1][i2]=xx[i];
					//Aquad[iz][i2][i1]=xx[i];
				i+=1;
			}
		}
		i=(NX+1)*NX/2;
		for(i1=0;i1<NX;i1++){
			Bquad[iz][i1]=xx[i];
			i+=1;
		}
		i=-1+(NX+1)*(NX+2)/2;
		Cquad[iz]=xx[i];
	}
	for(iz=0;iz<NX;iz++){
		printf("---------------- iz=%d ---------------------\n",iz);
		printf("Aquad\n");
		for(i1=0;i1<NX;i1++){
			for(i2=0;i2<NX;i2++){
				printf("%8.4f ",Aquad[iz][i1][i2]);
			}
			printf("\n");
		}
		printf("Bquad\n");
		for(i2=0;i2<NX;i2++){
			printf("%8.4f ",Bquad[iz][i2]);
		}
		for(i2=0;i2<NX;i2++){
			printf("%8.4f ",Bquad[iz][i2]);
		}
		printf("\n");
		printf("Cquad=%8.4f\n",Cquad[iz]);
	}
	
	
	for(irun=0;irun<NRUNS;irun++){
		ri=runinfo[irun];
		if(ri->good){
			error[irun]=0.0;
			for(iz=0;iz<NX;iz++){
				ri->zquad[iz]=0.0;
				for(i1=0;i1<NX;i1++){
					for(i2=0;i2<=i1;i2++){
						ri->zquad[iz]+=Aquad[iz][i1][i2]*ri->x[i1]*ri->x[i2];
					}
					ri->zquad[iz]+=Bquad[iz][i1]*ri->x[i1];
				}
				ri->zquad[iz]+=Cquad[iz];
				error[irun]+=pow(ri->zquad[iz]-ri->z[iz],2);
			}
			if(error[irun]>3.0){
				printf("----- irun=%3d: error^2=%g -----\n",irun,error[irun]);
				PrintX(ri);
				PrintY(ri);
				printf("yquad = %8.4f=?%8.4f %8.4f=?%8.4f %8.4f=?%8.4f, %8.4f=?%8.4f, %8.4f=?%8.4f, %8.4f=?%8.4f\n",ri->zquad[0],ri->z[0],ri->zquad[1],ri->z[1],ri->zquad[2],ri->z[2],ri->zquad[3],ri->z[3],ri->zquad[4],ri->z[4],ri->zquad[5],ri->z[5]);
			}
		}
	}
	
	delete gslmatrix;
	delete [] error;
}
#endif