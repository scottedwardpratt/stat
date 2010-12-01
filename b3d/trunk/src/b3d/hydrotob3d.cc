#ifndef __hydrotob3d_cc__
#define __hydrotob3d_cc__
#define __JOSH_FORMAT__

#include "b3d.h"
using namespace std;

void CHYDROtoB3D::Init(){
	int ires,iring;
	string inputfilename;
	double pi,ei,sigma2i,dedti,degen;
	T=parameter::getD(b3d->parmap,"HYDRO_TF",160.0);
	ETAMAX=b3d->ETAMAX;
	randy=b3d->randy;
	MC_Ntarget=0.0;
	MC_Nbar=0.0;
	MC_NWrite=0;
	nsample=b3d->NSAMPLE;
	reslist=b3d->reslist;
	nres=reslist->NResonances;
	//printf("epsilon_H=%g\n",intrinsic->epsilon);
	density=new double[nres];
	ring=new CRing[b3d->NRINGSMAX+1];
	epsilon=P=0.0;
	CResInfo *resinfo=reslist->GfirstResInfoptr;
	for(ires=0;ires<nres;ires++){
		if(resinfo->code!=22){
			degen=2.0*resinfo->spin+1.0;
			//printf("ires=%d, ID=%d, degen=%g, mass=%g, T=%g\n",ires,resinfo->code,degen,resinfo->mass,T);
			b3d->freegascalc_onespecies(resinfo->mass,T,pi,ei,density[ires],sigma2i,dedti);
			//printf("mass=%g, T=%g, pi=%g, ei=%g, di=%g\n",resinfo->mass,T,pi,ei,density[ires]);
			density[ires]*=degen;
			P+=degen*pi;
			epsilon+=degen*ei;
		}
		resinfo=resinfo->nextResInfoptr;
	}
	GetLambdaFact();
	initialization=true;
}

void CHYDROtoB3D::ReadInput(){
	if(!initialization) Init();
	inputfilename="output/"+b3d->run_name+"/"+b3d->qualifier+"/freezeout.dat";
	printf("freezeout info from %s\n",inputfilename.c_str());
	input=fopen(inputfilename.c_str(),"r");
	ReadHeader(input);
	nrings=0;
	do{
		ring[nrings].Read(input);
		
			//fscanf(input,"%lf",&later->tau);
		if(ring[nrings].nread==0){
			ring[nrings].tau=ring[nrings-1].tau;
		}
		nrings++;
	}while(ring[nrings-1].nread!=0);
	fclose(input);
}

int CHYDROtoB3D::MakeEvent(){
	b3d->Reset();
	int nparts=0,iring;
	MC_NWrite=0;
	MC_Nbar=0.0;
	Ncheck=0.0;
	nparts=0;
	MC_Ntarget=-log(randy->ran())/nsample;
	for(iring=1;iring<nrings;iring++){
		nparts+=GenerateParticles(iring-1,iring);
	}
	//printf("nparts=%d, MC_Nbar=%g\n",nparts,MC_Nbar);
	return nparts;
}

int CHYDROtoB3D::GenerateParticles(int iring1,int iring2){
	CRing *earlier=&ring[iring1];
	CRing *later=&ring[iring2];
	//printf("iring1=%d, tau1=%g, iring2=%d, tau2=%g\n",iring1,earlier->tau,iring2,later->tau);
	CResInfo *resinfo;
	double dNbarmax;
	//MC_Ntarget-=log(randy->ran())/nsample;
	int ires,iquad;
	int iphi1,iphi2,alpha,beta,nphi=later->nphi,nparts=0;
	double V0,Vx,Vy,V,V2,Vmag,x1[3],x2[3],y1[3],y2[3],a[3],b[3];
	double smallestr,biggestr,cphi1,sphi1,cphi2,sphi2,phi1,phi2;
	double zsize,lambda[4][4],u[4]={1.0,0.0,0.0,0.0},pdotV,pdotu,udotV,wmax,gamma;
	double ptilde[4],p[4],pt,y,et,weight,minv,mass,eta,tau,wx,wy,x[4];
	int ID,ipart;
	for(iphi1=0;iphi1<nphi;iphi1++){
		//printf("iphi1=%d, nphi=%d, ",iphi1,nphi);
		//printf("r1=%g r2=%g\n",earlier->r[iphi1],later->r[iphi1]);
		iphi2=iphi1+1;
		//printf("-----------------------grep : iphi1=%d, iphi2=%d, nphi=%d\n",iphi1,iphi2,nphi);
		phi1=0.5*PI*iphi1/double(nphi);
		phi2=0.5*PI*iphi2/double(nphi);
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
		zsize=ETAMAX*(earlier->tau+later->tau);
		u[1]=0.25*(earlier->ux[iphi1]+earlier->ux[iphi2]+later->ux[iphi1]+later->ux[iphi2]);
		u[2]=0.25*(earlier->uy[iphi1]+earlier->uy[iphi2]+later->uy[iphi1]+later->uy[iphi2]);
		//u[1]=u[2]=0.0;
		u[0]=sqrt(1.0+u[1]*u[1]+u[2]*u[2]);
		V0=(a[1]*b[2]-b[1]*a[2])*zsize;
		Vx=-(a[2]*b[0]-b[2]*a[0])*zsize;
		Vy=-(a[0]*b[1]-b[0]*a[1])*zsize;
		V0*=4.0; Vx*=4.0; Vy*=4.0;
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
		resinfo=reslist->GfirstResInfoptr;
		for(ires=0;ires<nres;ires++){
			if(resinfo->code!=22){
				dNbarmax=density[ires]*wmax;
				MC_Nbar+=dNbarmax;
				while(MC_Nbar>MC_Ntarget){
					mass=resinfo->mass;
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
							wx=randy->ran();
							wy=1.0-wx;
							for(alpha=1;alpha<3;alpha++){
								x[alpha]=0.5*(wx*x1[alpha]+wx*x2[alpha]+wy*y1[alpha]+wy*y2[alpha]);
								x[alpha]+=(0.5-randy->ran())*a[alpha];
							}
							eta=ETAMAX-2.0*randy->ran()*ETAMAX;
							tau=wx*later->tau+wy*earlier->tau;
							x[0]=tau*cosh(eta);
							x[3]=tau*sinh(eta);
							//printf("success: p=(%g,%g,%g,%g), x=(%g,%g,%g,%g), tau=%g, eta=%g\n",p[0],p[1],p[2],p[3],x[0],x[1],x[2],x[3],tau,eta);
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
							iquad=lrint(floor(4*randy->ran()));
							if(iquad==1 || iquad==2){
								x[1]=-x[1];
								p[1]=-p[1];
							}
							if(iquad==2 || iquad==3){
								x[2]=-x[2];
								p[2]=-p[2];
							}
							ID=resinfo->code;
							ipart=int(b3d->PartMap.size());
							//printf("%d: ID=%d, x=%g, y=%g, tau=%g, eta=%g, px=%g, py=%g, mass=%g, y=%g\n",ipart,ID,x[1],x[2],tau,eta,p[1],p[2],mass,y);
							b3d->part[ipart]->Init(ID,x[1],x[2],tau,eta,p[1],p[2],mass,y);
					//printf("y=%g, eta=%g, y-eta=%g\n",y,eta,y-eta);
						}
					}
				//else{
					//printf("pdotV < 0, =%g\n",pdotV);
				//}
					MC_Ntarget-=log(randy->ran())/nsample;
				}
			}
			resinfo=resinfo->nextResInfoptr;
		}
	}
	return nparts;
}

void CHYDROtoB3D::freegascalc_onespecies(double m,double t,double &p,double &e,double &dens,double &sigma2,double &dedt){
	const double prefactor=1.0/(2.0*PI*PI*pow(HBARC,3));
	double k0,k1,z,k0prime,k1prime,m2,m3,m4,t2,t3,I1,I2,Iomega;
	m2=m*m;
	m3=m2*m;
	m4=m2*m2;
	t2=t*t;
	t3=t2*t;
	z=m/t;
	if(z>1000.0){
		p=e=dens=dedt=0.0;
		printf("z is huge=%g, m=%g, t=%g\n",z,m,t);
	}
	else{
		if(z<0.0){
			printf("___z=%g,m=%g,T=%g ___\n",z,m,t);
			exit(1);
		}
		k0=Bessel::K0(z);
		k1=Bessel::K1(z);
		p=prefactor*(m2*t2*k0+2.0*m*t3*k1);
		e=prefactor*(3.0*m2*t2*k0+(m3*t+6.0*m*t3)*k1);
		dens=p/t;
		k0prime=-k1;
		k1prime=-k0-k1/z;
		dedt=prefactor*(6.0*m2*t*k0+(m3+18.0*m*t2)*k1-3.0*m3*k0prime-((m4/t)+6.0*m2*t)*k1prime);
		Iomega=exp(-m/t)/(30.0*PI*PI*HBARC*HBARC*HBARC);
		I1=pow(m,1.5)*pow(t,3.5)*7.5*sqrt(2.0*PI);
		I2=24.0*pow(t,5);
		sigma2=Iomega*(I1+I2+0.5*sqrt(I1*I2));  // this is an approximation (+/-2%) to messy integral
	}
}

// \lambda_{ij}=\delta T_{ij}/(\epsilon+P) \frac{1}{\lambda_{\rm fact}}

double CHYDROtoB3D::GetLambdaFact(){
	int ires,iQ,n,i;
	const int nmax=70;
	double G[nmax+5];
	double m,degen,z,Ipp=0.0,dIpp,J,nfact,sign,alpha;
	CResInfo *resinfo=reslist->GfirstResInfoptr;
	for(ires=0;ires<nres;ires++){
		if(resinfo->code!=22){
			m=resinfo->mass;
			degen=(2.0*resinfo->spin+1);
			z=m/T;
			alpha=0.0;

			G[0]=gsl_sf_gamma_inc(5,z)*pow(z,-5);
			for(int i=1;i<nmax+5;i++){
				n=5.0-2.0*i;
				if(n!=-1)	G[i]=(-exp(-z)/n)+(G[i-1]*z*z-z*exp(-z))/((n+1.0)*n);
				else G[i]=gsl_sf_gamma_inc(-1,z)*z;
			}
			J=0.0;
			nfact=1.0;
			sign=1.0;
			for(n=0;n<nmax;n+=1.0){
				if(n>0) sign=-1.0;
				J+=sign*nfact*(G[n]-2.0*G[n+1]+G[n+2]);
				nfact=nfact*0.5/(n+1.0);
				if(n>0) nfact*=(2.0*n-1.0);
			}
			Ipp+=degen*exp(alpha)*pow(m,4)*(-z*J+15.0*gsl_sf_bessel_Kn(2,z)/(z*z));
		}
		resinfo=resinfo->nextResInfoptr;
	}
	lambdafact=Ipp/(60.0*PI*PI*HBARC*HBARC*(P+epsilon));
}

void CHYDROtoB3D::ReadHeader(FILE *fptr){
	char dummy[180];
	string dummystring;
	do{
		fgets(dummy,120,fptr);
		dummystring.assign(dummy);
		//printf("%s",dummy);
	}while(dummystring!="END_OF_HEADER\n");
	fscanf(input,"%s",dummy);
	dummystring.assign(dummy);
	if(dummystring!="time:"){
		printf("YIKES, dummystring should have been TIME\n");
		exit(1);
	}
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

void CRing::Read(FILE *fptr){
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
	fscanf(fptr,"%lf",&x[0]);
	do{
		i+=1;
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &y[i],&epsilon,&P,&T,&Rqgp,&uxi[i],&uyi[i],&nt,&nx,&ny,&a1,&a2,&a3,&a4,&a5);
		fgets(dummy,100,fptr);
		//FillOutPi(PiOverh,lambda_xx[i],lambda_xy[i],lambda_yy[i],uxi[i],uyi[i]);
		lambda_xx[i]=(b+a1+a2/root3)/(epsilon+P);
		lambda_xy[i]=a3/(epsilon+P);
		lambda_yy[i]=(b-a1+a2/root3)/(epsilon+P);
		fscanf(fptr,"%s",dummy);
		dummystring.assign(dummy);
		if(dummystring!="time:") x[i+1]=atof(dummy);
	}while(dummystring!="time:" && feof(fptr)==false);
	nread=i;
	FillOutArrays(x,y,uxi,uyi,lambda_xx,lambda_xy,lambda_yy);
	//if(nread>0) printf("TIME=%g, %d points for Cooper-Frye ring\n",tau,nread);
}

void CRing::FillOutArrays(double *x,double *y,double *uxi,double *uyi,double *lambda_xx,double *lambda_xy,double *lambda_yy){
	//char filename[120];
	//FILE *ringinfo;
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

#endif
