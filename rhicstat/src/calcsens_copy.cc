#ifndef __PCA_CC__
#define __PCA_CC__
#include "coral.h"
#include "rhicstat.h"
using namespace std;

void CRHICStat::CalcSensitivity(){
	printf("check a\n");
	int ix,iy,iyy,ixx,irun;
	/** Calculate PCA values */
	for(iy=0;iy<NY;iy++){
		for(iyy=0;iyy<NY;iyy++) sigmayy[iy][iyy]=0.0;
	}
	for(iy=0;iy<NY;iy++){
		for(iyy=0;iyy<NY;iyy++){
			for(irun=0;irun<NRUNS;irun++){
				sigmayy[iy][iyy]+=runinfo[irun]->y[iy]*runinfo[irun]->y[iyy];
			}
			sigmayy[iy][iyy]=sigmayy[iy][iyy]/double(NRUNS);
		}
	}
	printf("check b\n");
	gslmatrix_NY->EigenFind(sigmayy,Uytoz,eigenvalyy);
	printf("check c\n");
	for(irun=0;irun<NRUNS;irun++){
		for(ix=0;ix<NX;ix++){
			runinfo[irun]->z[ix]=0.0;
			for(iy=0;iy<NY;iy++){
				runinfo[irun]->z[ix]+=runinfo[irun]->y[iy]*Uytoz[iy][ix];
			}
		}
		for(ix=NX;ix<NY;ix++) runinfo[irun]->z[ix]=0.0;
	}
	printf("check z\n");
	printf("-------- PCA values -----------\n");
	for(iy=0;iy<NY;iy++) printf("%2d: %g\n",iy,eigenvalyy[iy]);
	for(iy=0;iy<NY;iy++){
		printf("%2d %40s: %8.5f %8.5f %8.5f %8.5f\n", iy,yname[iy].c_str(),Uytoz[iy][0],Uytoz[iy][1],Uytoz[iy][2],Uytoz[iy][3]);
	}

	/** Check dxdx */
	double **dxdx,**dxdx_inv;
	dxdx=new double *[NX];
	dxdx_inv=new double *[NX];
	for(ix=0;ix<NX;ix++){
		dxdx[ix]=new double[NX];
		for(ixx=0;ixx<NX;ixx++){
			for(irun=0;irun<NRUNS;irun++){
				dxdx[ix][ixx]+=runinfo[irun]->x[ix]*runinfo[irun]->x[ixx];
			}
			dxdx[ix][ixx]=dxdx[ix][ixx]/double(NRUNS);
		}
	}
	//gslmatrix_NX->Invert(dxdx,dxdx_inv);
	printf("unit matrx ????\n");
	gslmatrix_NX->Print(dxdx);


	/** Calculate dxdz */
	for(ix=0;ix<NX;ix++){
		for(ixx=0;ixx<NX;ixx++) dxdz[ix][ixx]=0.0;
	}
	for(ix=0;ix<NX;ix++){
		for(ixx=0;ixx<NX;ixx++){
			for(irun=0;irun<NRUNS;irun++){
				dxdz[ix][ixx]+=runinfo[irun]->x[ix]*runinfo[irun]->z[ixx];
			}
			dxdz[ix][ixx]=dxdz[ix][ixx]/double(NRUNS);
		}
	}
	gslmatrix_NX->Invert_NonSymm(dxdz,dxdz_inv);

	for(ix=0;ix<NX;ix++){
		for(ixx=0;ixx<NX;ixx++){
			sigmaxx[ix][ixx]=0.0;
			for(iy=0;iy<NX;iy++){
				sigmaxx[ix][ixx]+=dxdz_inv[iy][ix]*dxdz_inv[iy][ixx];
			}
		}
	}
	//gslmatrix_NX->Print(sigmaxx);
	gslmatrix_NX->EigenFind(sigmaxx,Uxtow,eigenvalxx);
	
	
	printf("\n-------------------------------------------------\n");
	printf("error matrix eigenval.s: ");
	for(ix=0;ix<NX;ix++) printf("%9.5f ",eigenvalxx[ix]);
	printf("\n-------------------------------------------------\n");
	gslmatrix_NX->Invert_NonSymm(Uxtow,Uxtow_inv);
	
	double totalerror;
	printf("$$$$$$$$$$ UNCERTAINTIES $$$$$$$$$$\n");
	double *uncertainty=new double[NX];
	for(ix=0;ix<NX;ix++){
		uncertainty[ix]=0.0;
		for(ixx=0;ixx<NX;ixx++){
			totalerror=eigenvalxx[ixx]/(eigenvalxx[ixx]+1.0);
			uncertainty[ix]+=Uxtow_inv[ixx][ix]*Uxtow_inv[ixx][ix]*totalerror;
		}
		printf("%20s: %8.4f\n",xname[ix].c_str(),uncertainty[ix]);
	}
	
	printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
	int iz;
	/** double **mtest=new double *[NX];
	for(ix=0;ix<NX;ix++) mtest[ix]=new double[NX];
	for(ix=0;ix<NX;ix++){
		for(ixx=0;ixx<NX;ixx++){
			mtest[ix][ixx]=0.0;
			for(iz=0;iz<NX;iz++) mtest[ix][ixx]+=dxdz[ix][iz]*dxdz_inv[iz][ixx];
			printf("%8.4f ",mtest[ix][ixx]);
		}
		printf("\n");
	}
	for(ix=0;ix<NX;ix++) delete [] mtest[ix];
	delete [] mtest; */
	
	double chiperfreedom=0.0;
	for(irun=0;irun<NRUNS;irun++){
		totalerror=0.0;
		
		
		//PrintX(runinfo[irun]);
		GetZlinearFromX(runinfo[irun]);
		//GetXlinearFromZlinear(runinfo[irun]);
		//PrintXlinear(runinfo[irun]);
		
		for(iz=0;iz<NX;iz++){
			totalerror+=(runinfo[irun]->z[iz]-runinfo[irun]->zlinear[iz])
				*(runinfo[irun]->z[iz]-runinfo[irun]->zlinear[iz]);
		}
		chiperfreedom+=totalerror;
		if(totalerror>5.0){
			printf("------------------------  irun=%d ------------------------------\n",irun);
			printf("%7.3f: ",totalerror);
			for(iz=0;iz<NX;iz++) printf("%6.3f ",runinfo[irun]->zlinear[iz]-runinfo[irun]->z[iz]);
			printf("\nXvalues: ");
			for(ix=0;ix<NX;ix++) printf("%6.3f ",runinfo[irun]->x[ix]/sqrt(12.0));
			printf("\n");
			//PrintXlinear(runinfo[irun]);
			//PrintX(runinfo[irun]);
		}
	}
	//printf("chiperfreedom=%g\n",chiperfreedom/double(NRUNS*3)); // roughly 3 PCA values 

}

void CRHICStat::GetZFromY(CRunInfo *runinfo){
	int iy,iz;
	for(iz=0;iz<NX;iz++){
		runinfo->z[iz]=0.0;
		for(iy=0;iy<NY;iy++){
			runinfo->z[iz]+=runinfo->y[iy]*Uytoz[iy][iz];
		}
	}
}

void CRHICStat::GetYFromZ(CRunInfo *runinfo){
	int iy,iz;
	double *zz=new double[NY];
	for(iz=0;iz<NX;iz++) zz[iz]=runinfo->z[iz];
	for(iz=NX;iz<NY;iz++) zz[iz]=0.0;
	for(iy=0;iy<NY;iy++){
		runinfo->y[iy]=0.0;
		for(iz=0;iz<NY;iz++){
			runinfo->y[iy]+=zz[iz]*Uytoz_inv[iz][iy];
		}
	}
	delete [] zz;
}

void CRHICStat::GetYlinearFromZlinear(CRunInfo *runinfo){
	int iy,iz;
	double *zz=new double[NY];
	for(iz=0;iz<NX;iz++) zz[iz]=runinfo->zlinear[iz];
	for(iz=NX;iz<NY;iz++) zz[iz]=0.0;
	for(iy=0;iy<NY;iy++){
		runinfo->y[iy]=0.0;
		for(iz=0;iz<NY;iz++){
			runinfo->ylinear[iy]+=zz[iz]*Uytoz_inv[iz][iy];
		}
	}
	delete [] zz;
}



void CRHICStat::GetXlinearFromZ(CRunInfo *runinfo){
	int ix,iz;
	for(ix=0;ix<NX;ix++){
		runinfo->xlinear[ix]=0.0;
		for(iz=0;iz<NX;iz++){
			runinfo->xlinear[ix]+=runinfo->z[iz]*dxdz_inv[iz][ix];
		}
	}
}

void CRHICStat::GetXlinearFromY(CRunInfo *runinfo){
	GetZFromY(runinfo);
	GetXlinearFromZ(runinfo);
}

void CRHICStat::GetYlinearFromX(CRunInfo *runinfo){
	GetZlinearFromX(runinfo);
	GetYlinearFromZlinear(runinfo);
}

void CRHICStat::PrintXlinear(CRunInfo *runinfo){
	int ix;
	printf("      parameter name  : xlinearprime xlinear\n");
	for(ix=0;ix<NX;ix++){
		printf("%25s: %9.4f  %g\n",xname[ix].c_str(),
			runinfo->xlinear[ix],0.5*(xmin[ix]+xmax[ix])+runinfo->xlinear[ix]*(xmax[ix]-xmin[ix])/sqrt(12.0));
	}
}
void CRHICStat::PrintX(CRunInfo *runinfo){
	int ix;
	printf("       parameter name  :       xprime      x\n");
	for(ix=0;ix<NX;ix++){
		printf("%25s: %9.4f  %g\n",xname[ix].c_str(),
			runinfo->x[ix],0.5*(xmin[ix]+xmax[ix])+runinfo->x[ix]*(xmax[ix]-xmin[ix])/sqrt(12.0));
	}
}

void CRHICStat::GetZlinearFromX(CRunInfo *runinfo){
	int ix,iz;
	for(iz=0;iz<NX;iz++){
		runinfo->zlinear[iz]=0.0;
		for(ix=0;ix<NX;ix++){
			runinfo->zlinear[iz]+=runinfo->x[ix]*dxdz[ix][iz];
		}
		//for(iz=NX;iz<NY;iz++) runinfo->zlinear[iz]=0.0;
	}
}

void CRHICStat::GetXlinearFromZlinear(CRunInfo *runinfo){
	int ix,iz;
	for(ix=0;ix<NX;ix++){
		runinfo->xlinear[ix]=0.0;
		for(iz=0;iz<NX;iz++){
			runinfo->xlinear[ix]+=runinfo->zlinear[iz]*dxdz_inv[iz][ix];
		}
	}
}



#endif