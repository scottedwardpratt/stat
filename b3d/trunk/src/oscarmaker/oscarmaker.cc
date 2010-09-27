#ifndef __oscarmaker_cc__
#define __oscarmaker_cc__
#define __JOSH_FORMAT__

#include "oscarmaker.h"
using namespace std;

COSCARmaker::COSCARmaker(double T_in,double etamax_in,string inputfilename_in,string oscarfilename_in){
	int ispecies;
	inputfilename=inputfilename_in;
	oscarfilename=oscarfilename_in;
	T=T_in;
	ETA_MAX=etamax_in;
	randy=new CRandom(-time(NULL));
	MC_Ntarget=0.0;
	MC_Nbar=0.0;
	MC_NWrite=0;
	MC_NSample=1.0;
	nevents=0;
	species=new CSpecies_StandardHadrons_Equil();
	//species=new CSpecies_PionsOnly();
	CIntrinsic *intrinsic=new CIntrinsic(species);
	intrinsic->T=T; intrinsic->ZeroMu();
	CEqofst *eqofst=new CEqofst();
	eqofst->FreeGasCalc_of_TMu(intrinsic);
	intrinsic->Print();
	//printf("epsilon_H=%g\n",intrinsic->epsilon);
	lambdafact=eqofst->GetLambdaFact(intrinsic);
	density=new double[species->Nspecies];
	//	inputfilename="data/freezeout.dat";
	//	oscarfilename="oscar.dat";
	tmpfilename="oscar.tmp";
	oscarfile=fopen(oscarfilename.c_str(),"w");
	fprintf(oscarfile,"OSC1997A\n");
	fprintf(oscarfile,"final_id_p_x\n");
	fprintf(oscarfile,"     HydroOut      tag    197     79    197     79 eqsp 0.2130E+05        1\n");
	//      	fprintf(oscarfile,"vh2 1.0\n");
	fclose(oscarfile);
	earlier=new CRing();
	later=new CRing();
	ptbar=0;
	ptbar_norm=0;

	for(ispecies=0;ispecies<species->Nspecies;ispecies++){
		density[ispecies]=intrinsic->density[ispecies];
		//printf("m=%g, density=%g\n",species->m[ispecies],intrinsic->density[ispecies]);
	}
	delete(eqofst);
}

int COSCARmaker::MakeEvent(){
	int nparts=0;
	MC_NWrite=0;
	MC_Nbar=0.0;
	Ncheck=0.0;
	MC_Ntarget=-log(randy->ran())/MC_NSample;
	CRing *ringswitch;
	tmpfile=fopen(tmpfilename.c_str(),"w");
	input=fopen(inputfilename.c_str(),"r");
#ifdef __JOSH_FORMAT__
	later->ReadHeader_Josh(input);
	later->Read_Josh(input);
#else
	later->Read_VH2(input);
#endif
	do{
		ringswitch=earlier;
		earlier=later;
		later=ringswitch;
#ifdef __JOSH_FORMAT__
		later->Read_Josh(input);
		//fscanf(input,"%lf",&later->tau);
#else
		later->Read_VH2(input);
#endif
		if(later->nread==0){
			later->tau=earlier->tau;
		}
		nparts+=GenerateParticles();
	}while(later->nread!=0);
	printf("MC_Nbar=%g\n",MC_Nbar);
	fclose(input);

	nevents+=1;
	oscarfile=fopen(oscarfilename.c_str(),"a");
	fprintf(oscarfile,"%d %d 0.0 0.0\n",nevents,MC_NWrite);
	fclose(oscarfile);
	fclose(tmpfile);

	string command="cat "+tmpfilename+" >> "+oscarfilename;
	system(command.c_str());
	command="rm -f "+tmpfilename;
	system(command.c_str());
	printf("Ncheck=%g\n",Ncheck);

	return nparts;
}

int COSCARmaker::GenerateParticles(){
	double dNbarmax;
	//MC_Ntarget-=log(randy->ran())/MC_NSample;
	int ispecies,iid,iquad;
	int iphi1,iphi2,alpha,beta,nphi=later->nphi,nparts=0;
	double V0,Vx,Vy,V,V2,Vmag,x1[3],x2[3],y1[3],y2[3],a[3],b[3];
	double smallestr,biggestr,cphi1,sphi1,cphi2,sphi2,phi1,phi2;
	double zsize,lambda[4][4],u[4]={1.0,0.0,0.0,0.0},pdotV,pdotu,udotV,wmax,gamma;
	double ptilde[4],p[4],pt,y,et,weight,minv,mass,eta,tau,wx,wy,x[4];
	int ID;
	for(iphi1=0;iphi1<nphi;iphi1++){
		iphi2=iphi1+1;
		//printf("-----------------------grep : iphi1=%d, iphi2=%d, nphi=%d\n",iphi1,iphi2,nphi);
#ifdef __JOSH_FORMAT__
		phi1=0.5*PI*iphi1/double(nphi);
		phi2=0.5*PI*iphi2/double(nphi);
#else
		if(iphi2==nphi) iphi2=0;
		phi1=2.0*PI*iphi1/double(nphi);
		phi2=2.0*PI*iphi2/double(nphi);
#endif
		cphi1=cos(phi1);
		sphi1=sin(phi1);
		cphi2=cos(phi2);
		sphi2=sin(phi2);
		x1[0]=later->tau; x1[1]=later->r[iphi1]*cphi1; x1[2]=later->r[iphi1]*sphi1;
		x2[0]=later->tau; x2[1]=later->r[iphi2]*cphi2; x2[2]=later->r[iphi2]*sphi2;
		y1[0]=earlier->tau; y1[1]=earlier->r[iphi1]*cphi1; y1[2]=earlier->r[iphi1]*sphi1;
		y2[0]=earlier->tau; y2[1]=earlier->r[iphi2]*cphi2; y2[2]=earlier->r[iphi2]*sphi2;
		for(alpha=1;alpha<4;alpha++){
			for(beta=1;beta<4;beta++)
				lambda[alpha][beta]=(0.25/lambdafact)*(later->lambda[iphi1][alpha][beta]+later->lambda[iphi2][alpha][beta]
				+earlier->lambda[iphi1][alpha][beta]+earlier->lambda[iphi2][alpha][beta]);
			lambda[alpha][alpha]+=1.0;
		}
		for(alpha=0;alpha<3;alpha++){
			a[alpha]=0.5*(x2[alpha]-x1[alpha]+y2[alpha]-y1[alpha]);
			b[alpha]=0.5*(x2[alpha]-y2[alpha]+x1[alpha]-y1[alpha]);
		}
		//a[0]=b[0]=0.0;
		zsize=ETA_MAX*(earlier->tau+later->tau);
		u[1]=0.25*(earlier->ux[iphi1]+earlier->ux[iphi2]+later->ux[iphi1]+later->ux[iphi2]);
		u[2]=0.25*(earlier->uy[iphi1]+earlier->uy[iphi2]+later->uy[iphi1]+later->uy[iphi2]);
		//u[1]=u[2]=0.0;
		u[0]=sqrt(1.0+u[1]*u[1]+u[2]*u[2]);
		V0=(a[1]*b[2]-b[1]*a[2])*zsize;
		Vx=-(a[2]*b[0]-b[2]*a[0])*zsize;
		Vy=-(a[0]*b[1]-b[0]*a[1])*zsize;
#ifdef __JOSH_FORMAT__
		V0*=4.0; Vx*=4.0; Vy*=4.0;
#endif
		V2=V0*V0-Vx*Vx-Vy*Vy;
		udotV=u[0]*V0-u[1]*Vx-u[2]*Vy;
		V=sqrt(fabs(V2));
		//printf("phi1=%g, phi2=%g, V=%g\n",phi1*180/PI,phi2*180/PI,V);
		//printf("iphi1=%d, V=%g=(%g,%g,%g), u=(%g,%g,%g)\n",iphi1,V,V0,Vx,Vy,u[0],u[1],u[2]);
		//printf("a=(%g,%g,%g), b=(%g,%g,%g)\n",a[0],a[1],a[2],b[0],b[1],b[2]);
		//printf("x1=(%g,%g,%g)\n",x1[0],x1[1],x1[2]);
		//printf("y1=(%g,%g,%g)\n",y1[0],y1[1],y1[2]);
		//printf("_________________________________\n");
		if(V2>0){
			gamma=fabs(udotV/V);
		}
		else{
			Vmag=sqrt(Vx*Vx+Vy*Vy);
			gamma=(u[0]*Vmag-(V0/Vmag)*(Vx*u[1]+Vy*u[2]))/V;
		}
		wmax=V*(gamma+sqrt(gamma*gamma-1.0));
		// For MC procedure, weights cannot exceed unity. These factors are added into densities and weights
		// to ensure MC weights do not exceed unity for any p. wmax is calculated so that the maximum weight
		// will be unity for a given u and V as p is varied. If wmax were larger, the answer would not change
		// but the efficiency would suffer. We check to make sure weights never exceed unity
		if(gamma<1.0){
			printf("Disaster: gamma =%g, but should be >1\n",gamma);
			printf("x1=(%g,%g,%g), x2=(%g,%g,%g), y1=(%g,%g,%g), y2=(%g,%g,%g)\n",
				x1[0],x1[1],x1[2],x2[0],x2[1],x2[2],y1[0],y1[2],y1[2],y2[0],y2[1],y2[2]);
			printf("V0=%g, Vx=%g, Vy=%g\n",V0,Vx,Vy);
			printf("V2=%g, V=%g, gamma=%g, udotV=%g wmax=%g\n",V2,V,gamma,udotV,wmax);
			printf("fabs(MC_Nbar)=%g\n",fabs(MC_Nbar));
			exit(1);
		}
		for(ispecies=0;ispecies<species->Nspecies;ispecies++){
			dNbarmax=density[ispecies]*wmax;
			//printf("dNbarmax=%g,wmax=%g\n",dNbarmax,wmax);
			//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
			mass=species->m[ispecies];
			randy->generate_boltzmann(mass,T,ptilde);
			for(alpha=1;alpha<4;alpha++){
				p[alpha]=0.0;
				for(beta=1;beta<4;beta++) p[alpha]+=ptilde[beta]*lambda[alpha][beta];
				p[0]=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
			}
			for(alpha=0;alpha<4;alpha++) ptilde[alpha]=p[alpha];
			Misc::Boost(u,ptilde,p);
			pdotu=p[0]*u[0]-p[1]*u[1]-p[2]*u[2];
			pdotV=p[0]*V0-p[1]*Vx-p[2]*Vy;
			if(pdotV>0.0)	Ncheck+=density[ispecies]*pdotV/pdotu;
			//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
			MC_Nbar+=dNbarmax;
			//printf("ispecies=%d, MC_Nbar=%g, dens=%g, wmax=%g\n",ispecies,MC_Nbar,density[ispecies],wmax);
			while(MC_Nbar>MC_Ntarget){
				mass=species->m[ispecies];
				randy->generate_boltzmann(mass,T,ptilde);
				for(alpha=1;alpha<4;alpha++){
					p[alpha]=0.0;
					for(beta=1;beta<4;beta++) p[alpha]+=ptilde[beta]*lambda[alpha][beta];
					p[0]=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
				}
				for(alpha=0;alpha<4;alpha++) ptilde[alpha]=p[alpha];
				minv=sqrt(ptilde[0]*ptilde[0]-ptilde[1]*ptilde[1]-ptilde[2]*ptilde[2]-ptilde[3]*ptilde[3]);
				Misc::Boost(u,ptilde,p);
				pt=sqrt(p[1]*p[1]+p[2]*p[2]);
				minv=sqrt(p[0]*p[0]-pt*pt-p[3]*p[3]);

				pdotu=p[0]*u[0]-p[1]*u[1]-p[2]*u[2];
				pdotV=p[0]*V0-p[1]*Vx-p[2]*Vy;
				if(pdotV>0.0){
					pdotu=p[0]*u[0]-p[1]*u[1]-p[2]*u[2];
					weight=(pdotV/pdotu)/wmax;
					//printf("weight=%g\n",weight);
					if(weight>1.00000001){
						printf("DISASTER, weight=%g\n",weight);
						printf("V0=%g, Vx=%g, Vy=%g\n",V0,Vx,Vy);
						printf("p=(%g,%g,%g,%g), V2=%g, pdotV=%g, V/udotV=%g\n",p[0],p[1],p[2],p[3],V2,pdotV,(V/udotV));
						exit(1);
					}
					if(randy->ran()<weight){
						MC_NWrite+=1;
						nparts+=1;
						//printf("success: p=(%g,%g,%g,%g)\n",p[0],p[1],p[2],p[3]);
						iid=rint(floor(species->NID[ispecies]*randy->ran()));
						wx=randy->ran();
						wy=1.0-wx;
						for(alpha=1;alpha<3;alpha++){
							x[alpha]=0.5*(wx*x1[alpha]+wx*x2[alpha]+wy*y1[alpha]+wy*y2[alpha]);
							x[alpha]+=(0.5-randy->ran())*a[alpha];
						}
						eta=ETA_MAX-2.0*randy->ran()*ETA_MAX;
						tau=wx*later->tau+wy*earlier->tau;
						x[0]=tau*cosh(eta);
						x[3]=tau*sinh(eta);
						et=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]);
						y=0.5*log((p[0]+p[3])/(p[0]-p[3]));
						y+=eta;
						p[0]=et*cosh(y);
						p[3]=et*sinh(y);
						if(fabs(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3]-mass*mass)>1.0E-2){
							printf("invariant mass screwed up, =%g !=%g\n",sqrt(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3]),mass);
							printf("p=(%g,%g,%g,%g)\n",p[0],p[1],p[2],p[3]);
							exit(1);
						}
						//printf("pid=%d, pt=%g, ux=%g, uy=%g, uhat=%g\n",species->ID[ispecies][iid],pt,u[1],u[2],(u[1]*p[1]+u[2]*p[2])/pt);
#ifdef __JOSH_FORMAT__
						iquad=lrint(floor(4*randy->ran()));
						if(iquad==1 || iquad==2){
							x[1]=-x[1];
							p[1]=-p[1];
						}
						if(iquad==2 || iquad==3){
							x[2]=-x[2];
							p[2]=-p[2];
						}
#endif
						fprintf(tmpfile,"%d %d %g %g %g %g %g %g %g %g %g\n", MC_NWrite,species->ID[ispecies][iid],p[1]*0.001,p[2]*0.001,p[3]*0.001,p[0]*0.001,mass*0.001,x[1],x[2],x[3],x[0]);
						ID=species->ID[ispecies][iid];
						if(ID==211 || ID==-211 || ID==111){
							ptbar+=sqrt(p[1]*p[1]+p[2]*p[2]);
							ptbar_norm+=1;
						}
							//printf("y=%g, eta=%g, y-eta=%g\n",y,eta,y-eta);
					}
				}
				//else{
					//printf("pdotV < 0, =%g\n",pdotV);
				//}
				MC_Ntarget-=log(randy->ran())/MC_NSample;
			}
		}
	}
	return nparts;
}

CRing::CRing(){
	nphi=36;
	int iphi,alpha,beta;
	r=new double[nphi+1];
	ux=new double[nphi+1];
	uy=new double[nphi+1];
	lambda=new double**[nphi+1];
	for(iphi=0;iphi<=nphi;iphi++){
		lambda[iphi]=new double*[4];
		for(alpha=0;alpha<4;alpha++){
			lambda[iphi][alpha]=new double[4];
			for(beta=0;beta<4;beta++){
				lambda[iphi][alpha][beta]=0.0;
			}
		}
	}
}

void CRing::Read_VH2(FILE *fptr){
	int i;
	string dummystring;
	char filename[100];
	char dummy[120];
	double *x,*y,*uxi,*uyi;
	x=new double[300]; y=new double[300]; uxi=new double[300]; uyi=new double[300];
	double u[4];
	double **PiOverh,**PiOverhtilde,Pixxoverh,Pixyoverh,Piyyoverh;
	double lambda_xx[300],lambda_xy[300],lambda_yy[300];
	double phi,cphi,sphi,dphi,dphi1,dphi2,dphi1_min,dphi2_min,xx,yy,r1,r2,w1,w2;
	PiOverhtilde=new double*[4];
	PiOverh=new double*[4];
	for(int alpha=0;alpha<4;alpha++){
		PiOverhtilde[alpha]=new double[4];
		PiOverh[alpha]=new double[4];
	}

	i=0;
	do{
		fscanf(fptr,"%s",dummy);
		dummystring.assign(dummy);
		if(dummystring=="TIME"){
			fscanf(fptr,"%lf",&tau);
		}
		else if(dummystring.length()>0 &&feof(fptr)==false){
			x[i]=atof(dummy);
			fscanf(fptr,"%lf %lf %lf %lf %lf %lf",&y[i],&uxi[i],&uyi[i],&Pixxoverh,&Pixyoverh,&Piyyoverh);
			FillOutPi(PiOverh,Pixxoverh,Pixyoverh,Piyyoverh,uxi[i],uyi[i]);
			u[1]=uxi[i]; u[2]=uyi[i]; u[0]=sqrt(1.0+u[1]*u[1]+u[2]*u[2]); u[3]=0.0;
			Misc::BoostToCM(u,PiOverh,PiOverhtilde);
			lambda_xx[i]=PiOverhtilde[1][1];
			lambda_xy[i]=PiOverhtilde[1][2];
			lambda_yy[i]=PiOverhtilde[2][2];
			i+=1;
		}
	}while(dummystring!="TIME" && feof(fptr)==false);
	nread=i;
	if(nread>0) FillOutArrays_VH2(x,y,uxi,uyi,lambda_xx,lambda_xy,lambda_yy);
	//if(nread>0) printf("TIME=%g, %d points for Cooper-Frye ring\n",tau,nread);
}

void CRing::ReadHeader_Josh(FILE *fptr){
	char dummy[180];
	string dummystring;
	do{
		fgets(dummy,120,fptr);
		dummystring.assign(dummy);
		//printf("%s",dummy);
	}while(dummystring!="END_OF_HEADER\n");
	fscanf(fptr,"%s",dummy);
	dummystring.assign(dummy);
	if(dummystring!="time:"){
		printf("YIKES, dummystring should have been TIME\n");
		exit(1);
	}
}

void CRing::Read_Josh(FILE *fptr){
	FILE *ringinfo;
	int i;
	string dummystring;
	char dummy[120];
	double *x,*y,*uxi,*uyi;
	x=new double[300]; y=new double[300]; uxi=new double[300]; uyi=new double[300];
	double **PiOverh,**PiOverhtilde,Pixxoverh,Pixyoverh,Piyyoverh;
	double lambda_xx[300],lambda_xy[300],lambda_yy[300];
	int i1,i2,iphi,alpha,beta;
	double nt,nx,ny,a1,a2,a3,a4,a5,b=0.0,epsilon,P,Rqgp,T;
	const double root3=sqrt(3.0);
	PiOverhtilde=new double*[4];
	PiOverh=new double*[4];
	for(alpha=0;alpha<4;alpha++){
		PiOverhtilde[alpha]=new double[4];
		PiOverh[alpha]=new double[4];
	}
	i=-1;		
	fscanf(fptr,"%lf",&tau);
	//printf("XXXXXXXXXX tau=%g XXXXXXXXXXXX\n",tau); //Misc::Pause();
	fscanf(fptr,"%lf",&x[0]);
	do{
		i+=1;
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &y[i],&epsilon,&P,&T,&Rqgp,&uxi[i],&uyi[i],&nt,&nx,&ny,&a1,&a2,&a3,&a4,&a5);
		fgets(dummy,100,fptr);
		/*
		printf("x=%g, y=%g\n",x[i],y[i]);
		printf("epsilon=%g, P=%g, T=%g, Rqgp=%g\n",epsilon,P,T,Rqgp);
		printf("u=(%g,%g)\n",uxi[i],uyi[i]);
		printf("n=(%g,%g,%g)\n",nt,nx,ny);
		*/
		//printf("uxi[%d]=%g, uyi[%d]=%g\n",i,uxi[i],i,uyi[i]);
		//FillOutPi(PiOverh,lambda_xx[i],lambda_xy[i],lambda_yy[i],uxi[i],uyi[i]);
		lambda_xx[i]=(b+a1+a2/root3)/(epsilon+P);
		lambda_xy[i]=a3/(epsilon+P);
		lambda_yy[i]=(b-a1+a2/root3)/(epsilon+P);
		//printf("lambda=(%g,%g,%g)\n",lambda_xx[i],lambda_xy[i],lambda_yy[i]);
		fscanf(fptr,"%s",dummy);
		dummystring.assign(dummy);
		if(dummystring!="time:") x[i+1]=atof(dummy);
	}while(dummystring!="time:" && feof(fptr)==false);
	nread=i;
	FillOutArrays_Josh(x,y,uxi,uyi,lambda_xx,lambda_xy,lambda_yy);
	//if(nread>0) printf("TIME=%g, %d points for Cooper-Frye ring\n",tau,nread);
}

void CRing::FillOutArrays_Josh(double *x,double *y,double *uxi,double *uyi,double *lambda_xx,double *lambda_xy,double *lambda_yy){
	char filename[120];
	FILE *ringinfo;
	double phi,cphi,sphi,dphi,dphi1,dphi2,dphi1_min,dphi2_min,xx,yy,r1,r2,w1,w2;
	int i,i1,i2,iphi,alpha,beta;
	if(nread>0){
		//sprintf(filename,"plotdata/ring_tau%g.dat",tau);
		//ringinfo=fopen(filename,"w");
		for(iphi=0;iphi<=nphi;iphi++){
			phi=iphi*0.5*PI/nphi;
			cphi=cos(phi);
			sphi=sin(phi);
			i1=i2=0;
			dphi1=10000*PI;
			dphi2=-10000*PI;
			for(i=0;i<nread;i++){
				xx=x[i]*cphi+y[i]*sphi;
				yy=x[i]*sphi-y[i]*cphi;
				dphi=atan2(yy,xx);
				if(dphi<0.0) dphi+=2.0*PI;
				if(dphi<dphi1){
					dphi1=dphi;
					i1=i;
				}
				dphi-=2.0*PI;
				if(dphi>dphi2){
					dphi2=dphi;
					i2=i;
				}
			}
			w1=fabs(dphi2)/(dphi1-dphi2);
			w2=1.0-w1;
			r1=sqrt(x[i1]*x[i1]+y[i1]*y[i1]);
			r2=sqrt(x[i2]*x[i2]+y[i2]*y[i2]);
			r[iphi]=w1*r1+w2*r2;
			ux[iphi]=w1*uxi[i1]+w2*uxi[i2];
			uy[iphi]=w1*uyi[i1]+w2*uyi[i2];
			lambda[iphi][1][1]=w1*lambda_xx[i1]+w2*lambda_xx[i2];
			lambda[iphi][1][2]=w1*lambda_xy[i1]+w2*lambda_xy[i2];
			lambda[iphi][2][1]=lambda[iphi][1][2];
			lambda[iphi][2][2]=w1*lambda_yy[i1]+w2*lambda_yy[i2];
			lambda[iphi][3][3]=-lambda[iphi][2][2]-lambda[iphi][1][1];
			//printf("iphi=%d, phi=%g, r=%g, ux=%g, uy=%g\n",iphi,phi*180.0/PI,r[iphi],ux[iphi],uy[iphi]);
			//fprintf(ringinfo,"%7.4f %8.5f %8.4f %8.4f %10.3e %10.3e %10.3e\n",phi*180.0/PI,r[iphi],ux[iphi],uy[iphi],lambda[iphi][1][1],lambda[iphi][2][2],lambda[iphi][1][2]);
		}
		//fclose(ringinfo);
		//exit(1);
	}
	else{
		for(iphi=0;iphi<=nphi;iphi++){
			ux[iphi]=uy[iphi]=0.0;
			r[iphi]=0.0;
			lambda[iphi][1][1]=0.0;
			lambda[iphi][1][2]=lambda[iphi][2][1]=0.0;
			lambda[iphi][2][2]=0.0;
			lambda[iphi][3][3]=0.0;
		}
	}
}

void CRing::FillOutArrays_VH2(double *x,double *y,double *uxi,double *uyi,double *lambda_xx,double *lambda_xy,double *lambda_yy){
	char filename[120];
	FILE *ringinfo;
	double phi,cphi,sphi,dphi,dphi1,dphi2,dphi1_min,dphi2_min,xx,yy,r1,r2,w1,w2;
	int i,i1,i2,iphi,alpha,beta;
	if(nread>0){
		//sprintf(filename,"plotdata/ring_tau%g.dat",tau);
		//ringinfo=fopen(filename,"w");
		for(iphi=0;iphi<=nphi;iphi++){
			phi=iphi*2.0*PI/nphi;
			cphi=cos(phi);
			sphi=sin(phi);
			i1=i2=0;
			dphi1=10000*PI;
			dphi2=-10000*PI;
			for(i=0;i<nread;i++){
				xx=x[i]*cphi+y[i]*sphi;
				yy=x[i]*sphi-y[i]*cphi;
				dphi=atan2(yy,xx);
				if(dphi<0.0) dphi+=2.0*PI;
				if(dphi<dphi1){
					dphi1=dphi;
					i1=i;
				}
				dphi-=2.0*PI;
				if(dphi>dphi2){
					dphi2=dphi;
					i2=i;
				}
			}
			w1=fabs(dphi2)/(dphi1-dphi2);
			w2=1.0-w1;
			r1=sqrt(x[i1]*x[i1]+y[i1]*y[i1]);
			r2=sqrt(x[i2]*x[i2]+y[i2]*y[i2]);
			r[iphi]=w1*r1+w2*r2;
			ux[iphi]=w1*uxi[i1]+w2*uxi[i2];
			uy[iphi]=w1*uyi[i1]+w2*uyi[i2];
			lambda[iphi][1][1]=w1*lambda_xx[i1]+w2*lambda_xx[i2];
			lambda[iphi][1][2]=w1*lambda_xy[i1]+w2*lambda_xy[i2];
			lambda[iphi][2][1]=lambda[iphi][1][2];
			lambda[iphi][2][2]=w1*lambda_yy[i1]+w2*lambda_yy[i2];
			lambda[iphi][3][3]=-lambda[iphi][2][2]-lambda[iphi][1][1];
			//printf("iphi=%d, phi=%g, r=%g, ux=%g, uy=%g\n",iphi,phi*180.0/PI,r[iphi],ux[iphi],uy[iphi]);
			//fprintf(ringinfo,"%7.4f %8.5f %8.4f %8.4f %10.3e %10.3e %10.3e\n",phi*180.0/PI,r[iphi],ux[iphi],uy[iphi],lambda[iphi][1][1],lambda[iphi][2][2],lambda[iphi][1][2]);
		}
		//fclose(ringinfo);
		//exit(1);
	}
	else{
		for(iphi=0;iphi<=nphi;iphi++){
			ux[iphi]=uy[iphi]=0.0;
			r[iphi]=0.0;
			lambda[iphi][1][1]=0.0;
			lambda[iphi][1][2]=lambda[iphi][2][1]=0.0;
			lambda[iphi][2][2]=0.0;
			lambda[iphi][3][3]=0.0;
		}
	}
}

void FillOutPi(double **Pi,double pixx,double pixy,double piyy,double ux,double uy){
	double u0=sqrt(1.0+ux*ux+uy*uy);
	Pi[0][0]=(ux*ux*pixx+2*ux*uy*pixy+uy*uy*piyy)/(u0*u0);
	Pi[0][1]=(ux*pixx+uy*pixy)/u0;
	Pi[0][2]=(ux*pixy+uy*piyy)/u0;
	Pi[0][3]=0.0;
	Pi[1][0]=Pi[0][1];
	Pi[1][1]=pixx;
	Pi[1][2]=pixy;
	Pi[1][3]=0.0;
	Pi[2][0]=Pi[0][2];
	Pi[2][1]=Pi[1][2];
	Pi[2][2]=piyy;
	Pi[2][3]=0.0;
	Pi[3][0]=0.0;
	Pi[3][1]=0.0;
	Pi[3][2]=0.0;
	Pi[3][3]=Pi[0][0]-Pi[1][1]-Pi[2][2];
	//printf("Trace Pi=%g =? 0\n",Pi[0][0]-Pi[1][1]-Pi[2][2]-Pi[3][3]);
}

#endif
