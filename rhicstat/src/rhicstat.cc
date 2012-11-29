#ifndef __RHICSTAT_CC__
#define __RHICSTAT_CC__
#include "rhicstat.h"
using namespace std;

CRHICStat *CRunInfo::rhicstat=NULL;

CRHICStat::CRHICStat(){
	parameter::ReadParsFromFile(parmap,"statinfo/statpars.dat");
	NRUNS=parameter::getI(parmap,"NRUNS",729);
	NTESTRUNS=parameter::getI(parmap,"NTESTRUNS",32);
	FIT_TYPE=parameter::getS(parmap,"FIT_TYPE","GP");
	NBURN=parameter::getD(parmap,"NBURN",10000);
	NMCMC=parameter::getD(parmap,"NMCMC",1000000);
	SIGMA2_EMULATOR=parameter::getD(parmap,"SIGMA2_EMULATOR",0.1);
//
	randy=new CRandom(-1234);
	int irun;
	InitArrays();
	bestinfo=new CRunInfo(NX,NY);
	runinfo=new CRunInfo *[NRUNS];
	for(irun=0;irun<NRUNS;irun++)
		runinfo[irun]=new CRunInfo(NX,NY);
	CRunInfo::rhicstat=this;
	testinfo=new CRunInfo *[NTESTRUNS];
	for(irun=0;irun<NTESTRUNS;irun++)
		testinfo[irun]=new CRunInfo(NX,NY);
	ReadAllX();
	
	ReadAllY();
	ScaleXY();
	PCA();
	GetZFromY(expinfo);
	if(FIT_TYPE=="GP"){
		zgetter=new CZGetter_GP(this);
		PerformFits();
	}
	if(FIT_TYPE=="QUAD"){
		zgetter=new CZGetter_QuadFit(this);
		PerformFits();
	}
	if(FIT_TYPE=="LOCALLINEAR"){
		zgetter =new CZGetter_LocalLinear(this);
		PerformFits();
	}
}

void CRHICStat::InitArrays(){
	InitX();
	InitY();
	int ix,iy,iz,irun1,irun2;
	eigenvalxx=new double[NX];
	uncertainty=new double[NX];
	eigenvalyy=new double[NY];
	sigmaybar=new double[NY];
	sigmaxx=new double *[NX];
	sigmayy=new double *[NY];
	Uytoz=new double *[NY];
	Uytoz_inv=new double *[NY];
	Uxtow=new double*[NX];
	Uxtow_inv=new double *[NX];
	dxdz=new double *[NX];
	dxdz_inv=new double *[NX];
	for(iy=0;iy<NY;iy++){
		sigmayy[iy]=new double[NY];
		Uytoz[iy]=new double[NY];
		Uytoz_inv[iy]=new double[NY];
	}
	for(ix=0;ix<NX;ix++){
		sigmaxx[ix]=new double[NX];
		Uxtow[ix]=new double[NX];
		Uxtow_inv[ix]=new double[NX];
		dxdz[ix]=new double[NX];
		dxdz_inv[ix]=new double[NX];
	}
	gslmatrix_NY=new CGSLMatrix_Real(NY);
	gslmatrix_NX=new CGSLMatrix_Real(NX);
	expinfo=new CRunInfo(NX,NY);
	fitinfo=new CRunInfo(NX,NY);

}

void CRHICStat::InitX(){
	char dummy[120];
	string filename;
	filename="ranges.dat";
	FILE *fptr=fopen(filename.c_str(),"r");
	/** First Get NX */
	NX=0;
	while(!feof(fptr)){
		fscanf(fptr,"%s",dummy);
		if(dummy[0]!='#' && !feof(fptr)){
			NX+=1;
		}
		fgets(dummy,200,fptr);
	}
	fclose(fptr);
	
	xname=new string[NX];
	xbar=new double[NX];
	xmin=new double[NX];
	xmax=new double[NX];
	xmcmc=NULL;
	fclose(fptr);
	/** Now go back and read in names and ranges */
	fptr=fopen(filename.c_str(),"r");
	int ix=0;
	while(!feof(fptr)){
		fscanf(fptr,"%s",dummy); // should just be reading in "double"
		if(dummy[0]!='#' && !feof(fptr)){
			if(string(dummy)!="double"){
				printf("in reading %s, expected to read \"double\"",filename.c_str());
				exit(1);
			}
			fscanf(fptr,"%s",dummy);
			xname[ix]=dummy;
			fscanf(fptr,"%lf %lf",&xmin[ix],&xmax[ix]);
			//if(ix==1){
				//xmin[ix]=0.0;
				//xmax[ix]=8.0;
			//}
			ix+=1;
		}
		fgets(dummy,200,fptr);
	}
	fclose(fptr);
	NZ=NX;
}

void CRHICStat::ReadAllX(){
	/** Read in parameter files */
	char runchars[10];
	string filename;
	int irun;
	for(irun=0;irun<NRUNS;irun++){
		sprintf(runchars,"%d",irun+1);
		filename="parameters/run"+string(runchars)+"/stats.param";
		ReadX(filename,runinfo[irun]);
	}
	for(irun=0;irun<NTESTRUNS;irun++){
		sprintf(runchars,"%d",irun+1);
		filename="testpars/run"+string(runchars)+"/stats.param";
		ReadX(filename,testinfo[irun]);
	}

}

void CRHICStat::ReadX(string filename,CRunInfo *runinfo){
	char dummy[200],runchars[5];
	int ix;
	FILE *fptr=fopen(filename.c_str(),"r");
	ix=0;
	while(!feof(fptr)){
		fscanf(fptr,"%s",dummy);
		if(dummy[0]!='#' && !feof(fptr)){
			fscanf(fptr,"%s",dummy);
			if(string(dummy)!=xname[ix]){
				printf("In reading file %s, xname[%d]=%s != %s\n",filename.c_str(),ix,xname[ix].c_str(),dummy);
				exit(1);
			}
			fscanf(fptr,"%lf",&(runinfo->x[ix]));
			ix+=1;
		}
		fgets(dummy,200,fptr);
	}
	fclose(fptr);
}

void CRHICStat::InitY(){
	string filename;
	char dummy[200],runchars[5];
	int i,irun,iy;
	filename="pcanames.dat";
	FILE *fptr=fopen(filename.c_str(),"r");
	/** First Get NY */
	NY=0;
	while(!feof(fptr)){
		fscanf(fptr,"%s",dummy);
		if(dummy[0]!='#' && !feof(fptr)){
			NY+=1;
		}
		fgets(dummy,200,fptr);
	}
	fclose(fptr);
	ybar=new double[NY];
	yname=new string[NY];
	/** Now get ynames */
	fptr=fopen(filename.c_str(),"r");
	iy=0;
	while(!feof(fptr)){
		fscanf(fptr,"%s",dummy);
		if(dummy[0]!='#' && !feof(fptr)){
			printf("yname[%d]=%s\n",iy,dummy);
			yname[iy]=dummy;
			iy+=1;
		}
		fgets(dummy,200,fptr);
	}
	fclose(fptr);
}

void CRHICStat::ReadAllY(){
	int irun,iy,ngood;
	char runchars[5];
	string filename;
		// Read Training Data
	NGOODRUNS=0;
	for(irun=1;irun<=NRUNS;irun++){
		sprintf(runchars,"%d",irun);
		filename="model_results/run"+string(runchars)+"/results.dat";
		ReadY(filename,runinfo[irun-1]);
		if(runinfo[irun-1]->good==false)
			printf("%d,",irun);
	}
	printf("\n");
	printf("NGOODRUNS=%d\n",NGOODRUNS);
	ngood=NGOODRUNS;
		//Read Testing Data
	for(irun=0;irun<NTESTRUNS;irun++){
		sprintf(runchars,"%d",irun+1);
		filename="test_results/run"+string(runchars)+"/results.dat";
		ReadY(filename,testinfo[irun]);
	}
		// Read Experimental Data
	ReadY("exp_data/results.dat",expinfo);
	NGOODRUNS=ngood; // don't count from test runs
	for(iy=0;iy<NY;iy++){
		sigmaybar[iy]=0.0;
		for(irun=0;irun<NRUNS;irun++){
			if(runinfo[irun]->good) sigmaybar[iy]+=runinfo[irun]->sigmay[iy];
		}
		sigmaybar[iy]=sigmaybar[iy]/double(NGOODRUNS);
	}	
}

void CRHICStat::ReadY(string filename,CRunInfo *runinfo){
	char dummy[200];
	runinfo->good=false;
	int iy;
	FILE *fptr=fopen(filename.c_str(),"r");
	iy=0;
	while(!feof(fptr)){
		fscanf(fptr,"%s",dummy);
		if(dummy[0]!='#' && !feof(fptr)){
			fscanf(fptr,"%s",dummy);
			if(string(dummy)==yname[iy]){
				fscanf(fptr,"%lf",&(runinfo->y[iy]));
				fscanf(fptr,"%lf",&(runinfo->sigmay[iy]));
				if(string(dummy)=="cent0to5_PHENIX_SPECTRA_PION_YIELD"){
					if(runinfo->y[iy]>0.0 && runinfo->y[iy]<100000.0){
						runinfo->good=true;
					}
					else{
						runinfo->good=false;
						printf("run from filename %s is useless, pion yield for central coll.s = %g\n",filename.c_str(),runinfo->y[iy]);
					}
				}
				iy+=1;
			}
		}
		fgets(dummy,200,fptr);
	}
	fclose(fptr);
	if(runinfo->good) NGOODRUNS+=1;
}

void CRHICStat::ScaleXY(){
	int ix,iy,irun;
	double ybartest;
	for(ix=0;ix<NX;ix++){
		xbar[ix]=0.0;
		for(irun=0;irun<NRUNS;irun++){
			xbar[ix]+=runinfo[irun]->x[ix];
		}
		xbar[ix]=xbar[ix]/double(NGOODRUNS);
		for(irun=0;irun<NRUNS;irun++){
			runinfo[irun]->x[ix]=sqrt(12.0)*(runinfo[irun]->x[ix]-xbar[ix])/(xmax[ix]-xmin[ix]);
		}
		for(irun=0;irun<NTESTRUNS;irun++){
			testinfo[irun]->x[ix]=sqrt(12.0)*(testinfo[irun]->x[ix]-xbar[ix])/(xmax[ix]-xmin[ix]);
		}
	}
	for(iy=0;iy<NY;iy++){
		ybar[iy]=0.0;
		for(irun=0;irun<NRUNS;irun++){
			if(runinfo[irun]->good){
				ybar[iy]+=runinfo[irun]->y[iy];
			}
		}
		ybar[iy]=ybar[iy]/double(NGOODRUNS);
		ybartest=0.0;
		for(irun=0;irun<NRUNS;irun++){
			runinfo[irun]->y[iy]=(runinfo[irun]->y[iy]-ybar[iy])/sigmaybar[iy];
			if(runinfo[irun]->good){
				ybartest+=runinfo[irun]->y[iy];
			}
		}
			//printf("ybartest=%12.5f=?0, ybar=%12.5f, sigmaybar=%g\n",ybartest,ybar[iy],sigmaybar[iy]);
		for(irun=0;irun<NTESTRUNS;irun++){
			testinfo[irun]->y[iy]=(testinfo[irun]->y[iy]-ybar[iy])/sigmaybar[iy];
		}
	}
	for(iy=0;iy<NY;iy++){
		expinfo->y[iy]=(expinfo->y[iy]-ybar[iy])/sigmaybar[iy];
		fitinfo->y[iy]=expinfo->y[iy];
	}
}

void CRHICStat::WritePars(string filename){
	int ix;
	FILE *fptr=fopen(filename.c_str(),"w");
	for(ix=0;ix<NX;ix++){
		fprintf(fptr,"double %s %g\n",xname[ix].c_str(),xmcmc[ix]);
	}
	fclose(fptr);
}

void CRHICStat::PerformFits(){
	CRunInfo *ri;
	int irun,iz;
	for(irun=0;irun<NRUNS;irun++){
		ri=runinfo[irun];
		zgetter->GetZ(ri->x,ri->zfit);
		GetYFromZ(ri);
		GetYfitFromZfit(ri);
		CalcNetDiffFit(ri);
		CalcNetDiffExp(ri);
		CalcNetDiffFitExp(ri);
		for(iz=0;iz<NZ;iz++){
			if(fabs(runinfo[irun]->z[iz]-runinfo[irun]->zfit[iz])>0.02)
				printf("z[%d]=%g =? %g\n",iz,runinfo[irun]->z[iz],runinfo[irun]->zfit[iz]);
		}
	}
	for(irun=0;irun<NTESTRUNS;irun++){
		ri=testinfo[irun];
		zgetter->GetZ(ri->x,ri->zfit);
		GetYFromZ(ri);
		GetYfitFromZfit(ri);
		CalcNetDiffFit(ri);
		CalcNetDiffExp(ri);
		CalcNetDiffFitExp(ri);
	}
}

#endif