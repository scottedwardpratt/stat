#ifndef __PCA_CC__
#define __PCA_CC__
#include "coral.h"
#include "rhicstat.h"
using namespace std;

CRHICStat::CRHICStat(int nruns_set,int ntestruns_set){
	NRUNS=nruns_set;
	NTESTRUNS=ntestruns_set;
	InitArrays();
	runinfo=new CRunInfo *[NRUNS];
	for(int irun=0;irun<NRUNS;irun++) runinfo[irun]=new CRunInfo(NX,NY);
	testinfo=new CRunInfo *[NTESTRUNS];
	for(int irun=0;irun<NTESTRUNS;irun++) testinfo[irun]=new CRunInfo(NX,NY);
	ReadAllX();
	ReadAllY();
	ScaleXY();
}

CRunInfo::CRunInfo(int NX,int NY){
	x=new double[NX];
	w=new double[NX];
	y=new double[NY];
	z=new double[NY];
	sigmay=new double[NY];
	xlinear=new double[NX];
	ylinear=new double[NY];
	zlinear=new double[NY];
	zquad=new double[NY]();
	yquad=new double[NY]();
	xquad=new double[NX]();
	good=true;
}

void CRHICStat::InitArrays(){
	InitX();
	InitY();
	int ix,iy;
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
	Aquad=NULL; Bquad=NULL; Cquad=NULL;
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
	printf("Useless runs: ");
	for(irun=0;irun<NRUNS;irun++){
		sprintf(runchars,"%d",irun+1);
		filename="model_results/run"+string(runchars)+"/results.dat";
		ReadY(filename,runinfo[irun]);
		if(runinfo[irun]->good==false) printf("%d,",irun+1);
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
	NGOODRUNS=ngood; // don't count from test runs
		// Read Experimental Data
	ReadY("exp_data/results.dat",expinfo);
	for(iy=0;iy<NY;iy++){
		sigmaybar[iy]=0.0;
		for(irun=0;irun<NRUNS;irun++){
			if(runinfo[irun]->good) sigmaybar[iy]+=runinfo[irun]->sigmay[iy];
		}
		sigmaybar[iy]=sigmaybar[iy]/double(NGOODRUNS);
		for(irun=0;irun<NRUNS;irun++) runinfo[irun]->sigmay[iy]=sigmaybar[iy];
		for(irun=0;irun<NTESTRUNS;irun++) testinfo[irun]->sigmay[iy]=sigmaybar[iy];
		fitinfo->sigmay[iy]=expinfo->sigmay[iy]=sigmaybar[iy];
		/** if(yname[iy]=="cent20to30_STAR_V2_PION_PTWEIGHT" || yname[iy]=="cent20to30_STAR_V2_KAON_PTWEIGHT" ||yname[iy]=="cent20to30_STAR_V2_PROTON_PTWEIGHT"){
			expinfo->y[iy]*=0.9;
		}*/
	}	
	
}
	
void CRHICStat::ReadY(string filename,CRunInfo *runinfo){
	char dummy[200];
	int iy;
	FILE *fptr=fopen(filename.c_str(),"r");
	iy=0;
	while(!feof(fptr)){
		fscanf(fptr,"%s",dummy);
		if(dummy[0]!='#' && !feof(fptr)){
			fscanf(fptr,"%s",dummy);
			if(string(dummy)==yname[iy]){
				fscanf(fptr,"%lf",&(runinfo->y[iy]));
				fscanf(fptr,"%lf",&runinfo->sigmay[iy]);
				if(string(dummy)=="cent0to5_PHENIX_SPECTRA_PION_YIELD"){
					if(runinfo->y[iy]>325.0 && runinfo->y[iy]<500.0){
						runinfo->good=true;
					}
					else{
						runinfo->good=false;
							//printf("run from filename %s is useless, pion yield for central coll.s = %g\n",filename.c_str(),runinfo->y[iy]);
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
	for(ix=0;ix<NX;ix++){
		xbar[ix]=0.0;
		for(irun=0;irun<NRUNS;irun++){
			if(runinfo[irun]->good) xbar[ix]+=runinfo[irun]->x[ix];
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
			if(runinfo[irun]->good) ybar[iy]+=runinfo[irun]->y[iy];
		}
		ybar[iy]=ybar[iy]/double(NGOODRUNS);
		for(irun=0;irun<NRUNS;irun++){
			runinfo[irun]->y[iy]=(runinfo[irun]->y[iy]-ybar[iy])/runinfo[irun]->sigmay[iy];
		}
		for(irun=0;irun<NTESTRUNS;irun++){
			testinfo[irun]->y[iy]=(testinfo[irun]->y[iy]-ybar[iy])/testinfo[irun]->sigmay[iy];
		}
	}
	for(iy=0;iy<NY;iy++){
		expinfo->y[iy]=(expinfo->y[iy]-ybar[iy])/expinfo->sigmay[iy];
		fitinfo->y[iy]=expinfo->y[iy];
	}
}

#endif