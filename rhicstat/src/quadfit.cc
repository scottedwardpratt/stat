#ifndef __QUADFIT_CC__
#define __QUADFIT_CC__
#include "coral.h"
#include "rhicstat.h"
using namespace std;

void CRHICStat::QuadFit(){
	int i,j,i1,i2,j1,j2,iz,irun,n=(NX+1)*(NX+2)/2;
	double **M=new double *[n];
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
							M[i][j]+=runinfo[irun]->x[i1]*runinfo[irun]->x[i2]*runinfo[irun]->x[j1]*runinfo[irun]->x[j2];
						}
						j+=1;
					}
				}
				j=NX*(NX+1)/2;
				for(j1=0;j1<NX;j1++){
					for(irun=0;irun<NRUNS;irun++){
						M[i][j]+=runinfo[irun]->x[i1]*runinfo[irun]->x[i2]*runinfo[irun]->x[j1];
					}
					j+=1;
				}
				j=-1+(NX+1)*(NX+2)/2;
				for(irun=0;irun<NRUNS;irun++){
					M[i][j]+=runinfo[irun]->x[i1]*runinfo[irun]->x[i2];
				}
				yy[i]=0.0;
				for(irun=0;irun<NRUNS;irun++){
					yy[i]+=runinfo[irun]->x[i1]*runinfo[irun]->x[i2]*runinfo[irun]->z[iz];
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
						M[i][j]+=runinfo[irun]->x[i1]*runinfo[irun]->x[j1]*runinfo[irun]->x[j2];
					}
					j+=1;
				}
			}
			j=NX*(NX+1)/2;
			for(j1=0;j1<NX;j1++){
				for(irun=0;irun<NRUNS;irun++){
					M[i][j]+=runinfo[irun]->x[i1]*runinfo[irun]->x[j1];
				}
				j+=1;
			}
			j=-1+(NX+1)*(NX+2)/2;
			for(irun=0;irun<NRUNS;irun++){
				M[i][j]+=runinfo[irun]->x[i1];
			}
			yy[i]=0.0;
			for(irun=0;irun<NRUNS;irun++){
				yy[i]+=runinfo[irun]->x[i1]*runinfo[irun]->z[iz];
			}
			i+=1;
		}
			// from dM/dC=0
		i=-1+(NX+1)*(NX+2)/2;
		j=0;
		for(j1=0;j1<NX;j1++){
			for(j2=0;j2<=j1;j2++){
				for(irun=0;irun<NRUNS;irun++){
					M[i][j]+=runinfo[irun]->x[j1]*runinfo[irun]->x[j2];
				}
				j+=1;
			}
		}
		j=NX*(NX+1)/2;
		for(j1=0;j1<NX;j1++){
			for(irun=0;irun<NRUNS;irun++){
				M[i][j]+=runinfo[irun]->x[j1];
			}
			j+=1;
		}
		j=-1+(NX+1)*(NX+2)/2;
		for(irun=0;irun<NRUNS;irun++){
			M[i][j]+=1.0;
		}
		yy[i]=0.0;
		for(irun=0;irun<NRUNS;irun++) yy[i]+=runinfo[irun]->z[iz];
		
		gslmatrix->SolveLinearEqs(yy,M,xx);
		
		i=0;
		for(i1=0;i1<NX;i1++){
			for(i2=0;i2<i1;i2++){
				Aquad[iz][i1][i2]=xx[i];
				Aquad[iz][i2][i1]=xx[i];
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
	
	for(irun=0;irun<NRUNS;irun++){
		for(iz=0;iz<NX;iz++){
			runinfo[irun]->zquad[iz]=0.0;
			for(i1=0;i1<NX;i1++){
				for(i2=0;i2<NX;i2++){
					runinfo[irun]->zquad[iz]+=Aquad[iz][i1][i2]*runinfo[irun]->x[i1]*runinfo[irun]->x[i2];
				}
				runinfo[irun]->zquad[iz]+=Bquad[iz][i1]*runinfo[irun]->x[i1];
			}
			runinfo[irun]->zquad[iz]+=Cquad[iz];
		}
	}
	delete gslmatrix;

}
#endif