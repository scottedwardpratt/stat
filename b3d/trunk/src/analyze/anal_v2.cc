#ifndef __ANALYZE_V2_CC__
#define __ANALYZE_V2_CC__
#include "analyze.h"
void CAnalyze::CalcV2(){
	const int NSPECIES=4;
	int ievent=1,ipart,nparts,nevents,ispecies,ID;
	double v21bar[NSPECIES]={0.0},v22bar[NSPECIES]={0.0},v23bar[NSPECIES]={0.0},v24bar[NSPECIES]={0.0},v2bar[NSPECIES]={0.0};
	double v2norm[NSPECIES]={0.0},v2norm1[NSPECIES]={0,0},v2norm2[NSPECIES]={0,0},v2norm3[NSPECIES]={0,0},v2norm4[NSPECIES]={0,0};
	double pt,v2;
	CPartH5 *ph5;
	printf("CAnalyze::ReadData, opening %s\n",h5_infilename.c_str());
	h5file = new H5File(h5_infilename,H5F_ACC_RDONLY);
	nevents=int(h5file->getNumObjs());
	if(nevents>neventsmax) nevents=neventsmax;
	for(ievent=1;ievent<=nevents;ievent++){
		nparts=ReadDataH5(ievent);
		for(ipart=0;ipart<nparts;ipart++){
			ph5=&partH5[ipart];
			pt=sqrt(ph5->px*ph5->px+ph5->py*ph5->py);
			v2=ph5->px*ph5->px-ph5->py*ph5->py;
			v2=v2/(pt*pt);
			ID=ph5->ID;
			ispecies=-1;
			if(ID==211 || ID==-211 || ID==111) ispecies=0;
			else if(ID==310 || ID==130 || ID==321 || ID==-321) ispecies=1;
			else if(ID==2112 || ID==2212 || ID==-2112 || ID==-2212) ispecies=2;
			else if(ID==3334 || ID==-3334) ispecies=3;
			if(ispecies>=0){
				v2bar[ispecies]+=v2;
				v21bar[ispecies]+=v2*pt;
				v22bar[ispecies]+=v2*pt*pt;
				v23bar[ispecies]+=v2*pt*pt*pt;
				v24bar[ispecies]+=v2*pt*pt*pt*pt;
				v2norm[ispecies]+=1.0;
				v2norm1[ispecies]+=pt;
				v2norm2[ispecies]+=pt*pt;
				v2norm3[ispecies]+=pt*pt*pt;
				v2norm4[ispecies]+=pt*pt*pt*pt*pt;
			}
		}
	}
	delete h5file;
	for(ispecies=0;ispecies<NSPECIES;ispecies++){
		v2bar[ispecies]=v2bar[ispecies]/v2norm[ispecies];
		v21bar[ispecies]=v21bar[ispecies]/v2norm1[ispecies];
		v22bar[ispecies]=v22bar[ispecies]/v2norm2[ispecies];
		v23bar[ispecies]=v23bar[ispecies]/v2norm3[ispecies];
		v24bar[ispecies]=v24bar[ispecies]/v2norm4[ispecies];
		printf("v2bar=%g, v2norm=%g\n",v2bar[ispecies],v2norm[ispecies]);
	}
}
#endif