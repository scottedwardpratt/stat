
void CRHICStat::AltCalcSensitivity(){
	int ix,iy,iyy,ixx,iz,izz,irun;
	for(ix=0;ix<NX;ix++){
		xbar[ix]=0.0;
		for(irun=0;irun<NRUNS;irun++){
			xbar[ix]+=x[ix][irun];
		}
		xbar[ix]=xbar[ix]/double(NRUNS);
		for(irun=0;irun<NRUNS;irun++){
			xprime[ix][irun]=sqrt(12.0)*(x[ix][irun]-xbar[ix])/(xmax[ix]-xmin[ix]);
		}
	}
	for(iy=0;iy<NY;iy++){
		ybar[iy]=0.0;
		for(irun=0;irun<NRUNS;irun++){
			ybar[iy]+=y[iy][irun];
		}
		ybar[iy]=ybar[iy]/double(NRUNS);
		for(irun=0;irun<NRUNS;irun++){
			yprime[iy][irun]=(y[iy][irun]-ybar[iy])/sigmay[iy];
			//printf("yprime[%d][%d]=%g, y=%g,ybar=%g, sigma=%g\n",iy,irun,yprime[iy][irun],y[iy][irun],ybar[iy],sigmay[iy]);
		}
	}
	// ***************************************************** */
	double **z=new double *[NX];
	double **dxdz=new double *[NX];
	double **dxdx=new double *[NX];
	double **dzdz=new double *[NX];
	double **dydy=new double *[NY];
	for(ix=0;ix<NX;ix++){
		dxdx[ix]=new double[NX];
		dxdz[ix]=new double[NX];
		dzdz[ix]=new double[NX];
		for(ixx=0;ixx<NX;ixx++) dxdx[ix][ixx]=dxdz[ix][ixx]=dzdz[ix][ixx]=0.0;
		z[ix]=new double[NRUNS];
	}
	for(iy=0;iy<NY;iy++){
		dydy[iy]=new double[NY];
		for(iyy=0;iyy<NY;iyy++) dydy[iy][iyy]=0.0;
	}
	
	for(iy=0;iy<NY;iy++){
		for(iyy=0;iyy<NY;iyy++){
			dydy[iy][iyy]=0.0;
			for(irun=0;irun<NRUNS;irun++){
				dydy[iy][iyy]+=yprime[iy][irun]*yprime[iyy][irun];
			}
			dydy[iy][iyy]=dydy[iy][iyy]/double(NRUNS);
		}
	}
	double *eigenvalyy=new double[NY];
	double **evecyy=new double *[NY];
	for(iy=0;iy<NY;iy++) evecyy[iy]=new double[NY];
	CGSLMatrix_Real *gslmatrix=new CGSLMatrix_Real(NY);
	gslmatrix->EigenFind(dydy,evecyy,eigenvalyy);
	printf("-------- PCA values -----------\n");
	for(iy=0;iy<NY;iy++) printf("%2d: %g\n",iy,eigenvalyy[iy]);
	
	for(irun=0;irun<NRUNS;irun++){
		for(iy=0;iy<NX;iy++){
			z[iy][irun]=0.0;
			for(iyy=0;iyy<NY;iyy++){
				z[iy][irun]+=yprime[iyy][irun]*evecyy[iyy][iy];
			}
		}
		for(ix=0;ix<NX;ix++){
			for(ixx=0;ixx<NX;ixx++){
				dzdz[ix][ixx]+=z[ix][irun]*z[ixx][irun];
			}
		}
	}
	for(ix=0;ix<NX;ix++){
		for(ixx=0;ixx<NX;ixx++) dzdz[ix][ixx]=dzdz[ix][ixx]/double(NRUNS);
	}
	printf("-------- PCA values -----------\n");
	for(iy=0;iy<NY;iy++) printf("%2d: %g\n",iy,eigenvalyy[iy]);
	
	delete gslmatrix;
	gslmatrix=new CGSLMatrix_Real(NX);
	printf("---------- dzdz ------------\n");
	gslmatrix->Print(dzdz);
	
	
	for(ix=0;ix<NX;ix++){
		for(ixx=0;ixx<NX;ixx++){
			for(irun=0;irun<NRUNS;irun++){
				dxdz[ix][ixx]+=xprime[ix][irun]*z[ixx][irun];
			}
			dxdz[ix][ixx]=dxdz[ix][ixx]/double(NRUNS);
		}
	}
	printf("00000000 dxdz 00000000\n");
	gslmatrix->Print(dxdz);
	
	double **dxdzinv=new double *[NX];
	double **error=new double *[NX];
	double **error_evec=new double *[NX];
	double *error_eigenval=new double[NX];
	double **ee_inv=new double *[NX];
	for(ix=0;ix<NX;ix++){
		dxdzinv[ix]=new double[NX];
		error[ix]=new double[NX];
		error_evec[ix]=new double [NX];
		ee_inv[ix]=new double[NX];
		for(ixx=0;ixx<NX;ixx++){
			dxdzinv[ix][ixx]=error[ix][ixx]=0.0;
			error_evec[ix][ixx]=ee_inv[ix][ixx]=0.0;
		}
	}
	gslmatrix->Invert_NonSymm(dxdz,dxdzinv);
	for(ix=0;ix<NX;ix++){
		for(ixx=0;ixx<NX;ixx++){
			for(iy=0;iy<NX;iy++){
				error[ix][ixx]+=dxdzinv[ix][iy]*dxdzinv[ixx][iy];
			}
		}
	}
	for(ix=0;ix<NX;ix++) printf("%2d: %20s\n",ix,xname[ix].c_str());
	printf("Error Matrix\n");
	gslmatrix->Print(error);
	gslmatrix->EigenFind(error,error_evec,error_eigenval);
	printf("\n-------------------------------------------------\n");
	printf("         eigenvalues: ");
	for(ix=0;ix<NX;ix++) printf("%9.5f ",error_eigenval[ix]);
	printf("\n-------------------------------------------------\n");
	for(ix=0;ix<NX;ix++){
		printf("%20s: ",xname[ix].c_str());
		for(ixx=0;ixx<NX;ixx++){
			printf("%9.5f ",error_evec[ix][ixx]);
		}
		printf("\n");
	}
	gslmatrix->Invert_NonSymm(error_evec,ee_inv);
	
	double totalerror;
	printf("$$$$$$$$$$ UNCERTAINTIES $$$$$$$$$$\n");
	double *uncertainty=new double[NX];
	for(ix=0;ix<NX;ix++){
		uncertainty[ix]=0.0;
		for(ixx=0;ixx<NX;ixx++){
			totalerror=error_eigenval[ixx]/(error_eigenval[ixx]+1.0);
			uncertainty[ix]+=ee_inv[ixx][ix]*ee_inv[ixx][ix]*totalerror;
		}
		printf("%20s: %8.4f\n",xname[ix].c_str(),uncertainty[ix]);
	}

}
