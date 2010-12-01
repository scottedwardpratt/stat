#ifndef __ANALYZE_HBT_CC__
#define __ANALYZE_HBT_CC__
#include "analyze.h"
void CAnalyze::CalcRx(){
	const int NSPECIES=4,NPTBINS=12;
	double taucutoff=40.0;
	int ievent=1,ipart,nparts,nevents,ispecies,ID,ipt;
	double Rx2bar[NSPECIES][NPTBINS]={0.0},Ry2bar[NSPECIES][NPTBINS]={0.0},Rz2bar[NSPECIES][NPTBINS]={0.0};
	double weight[NSPECIES][NPTBINS]={0.0};
	double pt,x,y,z,t,mass,tau;
	CPartH5 *ph5;
	string infilename=input_dataroot+"/"+qualifier+"/"+h5_infilename;
	printf("CAnalyze::ReadData, opening %s\n",infilename.c_str());
	h5file = new H5File(infilename,H5F_ACC_RDONLY);
	nevents=int(h5file->getNumObjs());
	if(nevents>neventsmax) nevents=neventsmax;
	for(ievent=1;ievent<=nevents;ievent++){
		nparts=ReadDataH5(ievent);
		for(ipart=0;ipart<nparts;ipart++){
			ph5=&partH5[ipart];
			ID=ph5->ID;
			//printf("listid=%d, ID=%d, Rx=%g, rapidity=%g\n",ph5->listid,ph5->ID,Rx,ph5->rapidity);
			ispecies=-1;
			if(ID==211 || ID==-211 || ID==111) ispecies=0;
			else if(ID==310 || ID==130 || ID==321 || ID==-321) ispecies=1;
			else if(ID==2112 || ID==2212 || ID==-2112 || ID==-2212) ispecies=2;
			else if(ID==3334 || ID==-3334) ispecies=3;
			if(ispecies>=0){
				pt=sqrt(ph5->px*ph5->px+ph5->py*ph5->py);
				GetHBTpars(ph5,tau,rout,rside,rlong);
				Rxbar[ispecies]+=Rx;
				Rx1bar[ispecies]+=Rx*pt;
				Rx2bar[ispecies]+=Rx*pt*pt;
				Rx3bar[ispecies]+=Rx*pt*pt*pt;
				Rx4bar[ispecies]+=Rx*pt*pt*pt*pt;
				Rxnorm[ispecies]+=1.0;
				Rxnorm1[ispecies]+=pt;
				Rxnorm2[ispecies]+=pt*pt;
				Rxnorm3[ispecies]+=pt*pt*pt;
				Rxnorm4[ispecies]+=pt*pt*pt*pt*pt;
			}
		}
	}
	delete h5file;
	printf("nevents=%d\n",nevents);
	for(ispecies=0;ispecies<NSPECIES;ispecies++){
		Rxbar[ispecies]=Rxbar[ispecies]/Rxnorm[ispecies];
		Rx1bar[ispecies]=Rx1bar[ispecies]/Rxnorm1[ispecies];
		Rx2bar[ispecies]=Rx2bar[ispecies]/Rxnorm2[ispecies];
		Rx3bar[ispecies]=Rx3bar[ispecies]/Rxnorm3[ispecies];
		Rx4bar[ispecies]=Rx4bar[ispecies]/Rxnorm4[ispecies];
		printf("Rxbar=%g, Rxnorm=%g\n",Rxbar[ispecies],Rxnorm[ispecies]);
	}
}

void CAnalyze::GetHBTpars(CPartH5 *ph5,double &tau,double &rout,double &rside,double &rlong){
	double taucompare=10.0;
	double x=ph5->x,y=ph5->y,px=ph5->px,py=ph5->py;
	tau=ph5->tau;
	rlong=tau*sinh(ph5->eta-ph5->rapidity);
	pt=sqrt(px*px+py*py);
	rout=(px*x+py*y)/pt;
	rside=(px*y-py*x)/pt;
}
#endif