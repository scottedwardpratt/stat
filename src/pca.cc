#ifndef __PCA_CC__
#define __PCA_CC__
#include "coral.h"
#include "pca.h"
using namespace std;

CPCA::CPCA(int nruns_set){
	string filename;
	char dummy[100],type[30];
	double dummyvalue;
	FILE *fptr;
	nruns=nruns_set;
	qualifiers.Read("qualifiers.dat");
	int iy,iqual;
	ny=0;
	for(iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		filename="analysis/run1/"+qualifiers.qualifier[iqual]+"/results.dat";
		printf("reading %s\n",filename.c_str());
		fptr=fopen(filename.c_str(),"r");
		do{
			fscanf(fptr,"%s",type);
			if(!feof(fptr)){
				if(type[0]!='#'){
					fscanf(fptr,"%s",dummy);
					yname[ny]=string(dummy)+"_"+qualifiers.qualifier[iqual];
					fscanf(fptr,"%lf",&dummyvalue);
					fscanf(fptr,"%lf",&sigmay[ny]);
					printf("yname[%d]=%s\n",ny,yname[ny].c_str());
					ny+=1;
				}
				else fgets(dummy,80,fptr);
			}
		}while(!feof(fptr));
		fclose(fptr);
	}
	ybar=new double[ny];
	value=new double[ny];
	spread=new double*[ny];
	for(iy=0;iy<ny;iy++) spread[iy]=new double[ny];
	gslmatrix=new CGSLMatrix_Real(ny);
	printf("CPCA initialized\n");
}

void CPCA::ReadResults(){
	FILE *fptr;
	string filename;
	char type[30],dummy[100];
	double dummysigma;
	int iy,jy,irun,iqual;
	for(irun=0;irun<nruns;irun++){
		iy=0;
		for(iqual=0;iqual<qualifiers.nqualifiers;iqual++){
			sprintf(dummy,"%d",irun+1);
			filename="analysis/run"+string(dummy)+"/"+qualifiers.qualifier[iqual]+"/results.dat";
			printf("reading %s\n",filename.c_str());
			fptr=fopen(filename.c_str(),"r");
			do{
				fscanf(fptr,"%s",type);
				if(!feof(fptr)){
					if(type[0]!='#'){
						fscanf(fptr,"%s",dummy);
						fscanf(fptr,"%lf",&value[iy]);
						fscanf(fptr,"%lf",&dummysigma);
						ybar[iy]+=value[iy];
						for(jy=0;jy<=iy;jy++){
							spread[iy][jy]+=value[iy]*value[jy];
						}
						iy+=1;
					}
					else fgets(dummy,80,fptr);
				}
			}while(!feof(fptr));
			fclose(fptr);
		}
	}
	for(iy=0;iy<ny;iy++){
		ybar[iy]=ybar[iy]/double(nruns);
		for(jy=0;jy<=iy;jy++){
			spread[iy][jy]=spread[iy][jy]/double(nruns);
			spread[iy][jy]=spread[iy][jy]-ybar[iy]*ybar[jy];
			spread[iy][jy]=spread[iy][jy]/(sigmay[iy]*sigmay[jy]);
			if(iy!=jy) spread[jy][iy]=spread[iy][jy];
		}
	}
	printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	for(iy=0;iy<ny;iy++) printf("ybar[%d]=%g, name=%s\n",iy,ybar[iy],yname[iy].c_str());
}

void CPCA::Calc(){
	int iy,jy;
	eigenval=new double[ny];
	evec=new double*[ny];
	for(iy=0;iy<ny;iy++) evec[iy]=new double[ny];
	
	CGSLMatrix_Real *gslmatrix=new CGSLMatrix_Real(ny);
	gslmatrix->EigenFind(spread,evec,eigenval);
	for(iy=0;iy<ny;iy++){
		printf("eigenval[%d]=%g\n",iy,eigenval[iy]);
	}
	
	/*
	printf("0000000000000 checking 0000000000000000\n");
	double *yy=new double[ny];
	for(int ky=0;ky<ny;ky++){
		printf("------- checking eigenvector %d -----------------\n",ky);
		for(iy=0;iy<ny;iy++){
			yy[iy]=0.0;
			for(jy=0;jy<ny;jy++){
				yy[iy]+=spread[iy][jy]*evec[jy][ky];
			}
			printf("%g =? %g\n",yy[iy],eigenval[ky]*evec[iy][ky]);
		}
	}
	*/
	
	double *respower=new double[ny];
	for(iy=0;iy<ny;iy++){
		respower[iy]=0.0;
		for(jy=0;jy<ny;jy++) respower[iy]+=eigenval[jy]*pow(evec[iy][jy],2);
		printf("respower[%d]=%10.3e, %10.3e, %s\n",iy,respower[iy],spread[iy][iy],yname[iy].c_str());
	}
	delete [] respower;
}

#endif