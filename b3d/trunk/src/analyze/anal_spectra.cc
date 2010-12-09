#ifndef __ANALYZE_SPECTRA_CC__
#define __ANALYZE_SPECTRA_CC__
#include "analyze.h"
void CAnalyze::CalcSpectra(){
	const int NSPECIES=4;
	int ievent=1,ipart,nparts,nevents,ispecies,ID;
	double pt2bar[NSPECIES]={0.0},pt3bar[NSPECIES]={0.0},pt4bar[NSPECIES]={0.0},ptbar[NSPECIES]={0.0},ptnorm[NSPECIES]={0.0},pt;
	CPartH5 *ph5;
	printf("CAnalyze::CalcSpectra, opening %s\n",h5_infilename.c_str());
	h5file=new H5File(h5_infilename,H5F_ACC_RDONLY);
	string outfilename=output_dataroot+"/spectra.dat";
	FILE *fptr=fopen(outfilename.c_str(),"w");
	
	nevents=int(h5file->getNumObjs());
	if(nevents>neventsmax) nevents=neventsmax;
	printf("nevents=%d\n",nevents);
	for(ievent=1;ievent<=nevents;ievent++){
		nparts=ReadDataH5(ievent);
		for(ipart=0;ipart<nparts;ipart++){
			ph5=&partH5[ipart];
			pt=sqrt(ph5->px*ph5->px+ph5->py*ph5->py);
			//printf("pt=%g\n",pt);
			//printf("listid=%d, ID=%d, pt=%g, rapidity=%g\n",ph5->listid,ph5->ID,pt,ph5->rapidity);
			ID=ph5->ID;
			ispecies=-1;
			if(ID==211 || ID==-211 || ID==111) ispecies=0;
			else if(ID==310 || ID==130 || ID==321 || ID==-321) ispecies=1;
			else if(ID==2112 || ID==2212 || ID==-2112 || ID==-2212) ispecies=2;
			else if(ID==3334 || ID==-3334) ispecies=3;
			if(ispecies>=0){
				ptbar[ispecies]+=pt;
				pt2bar[ispecies]+=pt*pt;
				pt3bar[ispecies]+=pt*pt*pt;
				pt4bar[ispecies]+=pt*pt*pt*pt;
				ptnorm[ispecies]+=1.0;
			}
		}
	}
	for(ispecies=0;ispecies<NSPECIES;ispecies++){
		ptbar[ispecies]=ptbar[ispecies]/ptnorm[ispecies];
		pt2bar[ispecies]=pt2bar[ispecies]/ptnorm[ispecies];
		pt3bar[ispecies]=pt3bar[ispecies]/ptnorm[ispecies];
		pt4bar[ispecies]=pt4bar[ispecies]/ptnorm[ispecies];
		printf("ptbar=%g, ptnorm=%g\n",ptbar[ispecies],ptnorm[ispecies]);
		fprintf(fptr,"ptbar=%g, ptnorm=%g\n",ptbar[ispecies],ptnorm[ispecies]);
	}
	delete h5file;
	fclose(fptr);
	//Misc::Pause();
}
#endif