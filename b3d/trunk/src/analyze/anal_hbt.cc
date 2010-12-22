#ifndef __ANALYZE_HBT_CC__
#define __ANALYZE_HBT_CC__
#include "analyze.h"
void CAnalyze::CalcHBT(){
	const int NSPECIES=3,NPTBINS=12;
	double taucutoff=20.0,wtau,delpt=50.0;
	int ievent=1,ipart,nparts,nevents,ispecies,ID,ipt;
	double Rx2bar[NSPECIES][NPTBINS]={0.0},Ry2bar[NSPECIES][NPTBINS]={0.0},Rz2bar[NSPECIES][NPTBINS]={0.0};
	double Rxbar[NSPECIES][NPTBINS]={0.0},Rybar[NSPECIES][NPTBINS]={0.0},Rzbar[NSPECIES][NPTBINS]={0.0};
	double Rnorm[NSPECIES][NPTBINS]={0.0},rootlambda[NSPECIES][NPTBINS]={0.0};
	double weight[NSPECIES][NPTBINS]={0.0};
	double pt,x,y,z,t,mass,tau,rout,rside,rlong,Rout,Rside,Rlong;
	CPartH5 *ph5;
	printf("CAnalyze::CalcHBT, opening %s\n",h5_infilename.c_str());
	h5file = new H5File(h5_infilename,H5F_ACC_RDONLY);
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
			//else if(ID==3334 || ID==-3334) ispecies=3;
			if(ispecies>=0){
				pt=sqrt(ph5->px*ph5->px+ph5->py*ph5->py);
				ipt=rint(floor(pt/delpt));
				if(ipt<NPTBINS){
					GetHBTpars(ph5,tau,rout,rside,rlong);
					wtau=exp(-0.5*pow((tau-10.0)/taucutoff,2));
					wtau*=exp(-0.5*pow((ph5->eta-ph5->rapidity)/2.0,2));
					Rnorm[ispecies][ipt]+=wtau;
					Rxbar[ispecies][ipt]+=wtau*rout;
					Rx2bar[ispecies][ipt]+=wtau*rout*rout;
					Rybar[ispecies][ipt]+=wtau*rside;
					Ry2bar[ispecies][ipt]+=wtau*rside*rside;
					Rzbar[ispecies][ipt]+=wtau*rlong;
					Rz2bar[ispecies][ipt]+=wtau*rlong*rlong;
					rootlambda[ispecies][ipt]+=1.0;
					if(ispecies==2 && ipt==0){
						printf("tau=%g, rout=%g, rside=%g, rlong=%g, wtau=%g\n",tau,rout,rside,rlong,wtau);
					}
				}
			}
		}
	}
	delete h5file;
	printf("nevents=%d\n",nevents);
	for(ispecies=0;ispecies<NSPECIES;ispecies++){
		printf("------------- ispecies=%d ---------------------\n",ispecies);
		printf("pt lambda    norm      <x>   <y>    <z>     Rout    Rside   Rlong Rout/Rside\n");
		for(ipt=0;ipt<NPTBINS;ipt++){
			pt=(0.5+ipt)*delpt;
			rootlambda[ispecies][ipt]=Rnorm[ispecies][ipt]/rootlambda[ispecies][ipt];
			Rxbar[ispecies][ipt]=Rxbar[ispecies][ipt]/Rnorm[ispecies][ipt];
			Rout=sqrt(fabs((Rx2bar[ispecies][ipt]/Rnorm[ispecies][ipt])-Rxbar[ispecies][ipt]*Rxbar[ispecies][ipt]));
			Rybar[ispecies][ipt]=Rybar[ispecies][ipt]/Rnorm[ispecies][ipt];
			Rside=sqrt(fabs((Ry2bar[ispecies][ipt]/Rnorm[ispecies][ipt])-Rybar[ispecies][ipt]*Rybar[ispecies][ipt]));
			Rzbar[ispecies][ipt]=Rzbar[ispecies][ipt]/Rnorm[ispecies][ipt];
			Rlong=sqrt(fabs((Rz2bar[ispecies][ipt]/Rnorm[ispecies][ipt])-Rzbar[ispecies][ipt]*Rzbar[ispecies][ipt]));
			printf("%3.0f %4.3f %8.2f %7.2f %6.3f %6.3f %7.2f %7.2f %7.2f %7.3f\n",
				(0.5+ipt)*delpt,rootlambda[ispecies][ipt]*rootlambda[ispecies][ipt],Rnorm[ispecies][ipt],
				Rxbar[ispecies][ipt],Rybar[ispecies][ipt],Rzbar[ispecies][ipt],
				Rout,Rside,Rlong,Rout/Rside);
		}
	}
}

void CAnalyze::GetHBTpars(CPartH5 *ph5,double &tau,double &rout,double &rside,double &rlong){
	double taucompare=10.0;
	double x=ph5->x,y=ph5->y,px=ph5->px,py=ph5->py,m=ph5->mass,pt,et,t;
	tau=ph5->tau;
	rlong=tau*sinh(ph5->eta-ph5->rapidity);
	t=tau*cosh(ph5->eta-ph5->rapidity);
	pt=sqrt(px*px+py*py);
	et=sqrt(pt*pt+m*m);
	rout=(px*x+py*y)/pt;
	rout=rout-(pt/et)*(t-taucompare);
	rside=(px*y-py*x)/pt;
}
#endif