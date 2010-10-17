#ifndef __ANALYZE_GAMMA_CC__
#define __ANALYZE_GAMMA_CC__
#include "analyze.h"
void CAnalyze::CalcGamma(){
	long long int naccept=0,dnaccept=0,npairs=0,MULT=0;
	const int NETA=20,NPHI=24;
	double gss[NETA]={0.0},gssprime[NETA]={0.0},gb[NETA]={0.0},gbprime[NETA]={0.0};
	double css[NETA],cssprime[NETA],cb[NETA],cbprime[NETA];
	long long int denom[NETA]={0};
	int ievent=1,ipart,jpart,jpartmax,nparts,nevents,ieta;
	double ptbar=0.0,pt,pt2,mt,m,pz;
	double sphi,cphi,sphitot,cphitot,sphi2tot,cphi2tot;
	double v2=0.0,v2prime=0.0,gammass=0.0,gammassprime=0.0,gammab=0.0,gammabprime=0.0;
	double sphij,cphij,deleta=0.1;
	double pxtot,px2tot,pytot,py2tot,px2totprime,pytot2prime,y,pmag,eta,pt2bar=0.0,etaj;
	CPartH5 *parti,*partj;
	string infilename=input_dataroot+"/"+qualifier+"/"+h5_infilename;
	printf("CAnalyze::ReadData, opening %s\n",infilename.c_str());
	h5file = new H5File(infilename,H5F_ACC_RDONLY);
	nevents=int(h5file->getNumObjs());
	if(nevents>neventsmax) nevents=neventsmax;
	for(ievent=1;ievent<=nevents;ievent++){
		cphitot=cphi2tot=sphitot=sphi2tot=0.0;
		pxtot=px2tot=pytot=py2tot=0.0;
		nparts=ReadDataH5(ievent);
		dnaccept=0;
		for(ipart=0;ipart<nparts;ipart++){
			parti=&partH5[ipart];
			pt2=parti->px*parti->px+parti->py*parti->py;
			pt=sqrt(pt2);
			m=parti->mass;
			mt=sqrt(pt2+m*m);
			y=parti->rapidity;
			pz=mt*sinh(y);
			pmag=sqrt(pt2+pz*pz);
			eta=atanh(pz/pmag);
			if(fabs(eta)<1.0) MULT+=1;
			if((!STAR_ACCEPTANCE) || (pt<2000 && fabs(eta)<1.0 && pt>200.0)){
	//if(fabs(y)<1.0 && pt>0.15 && pt<2.0){
				dnaccept+=1;
				ptbar+=pt;
				pt2bar+=pt2;
	//ptbar+=pt;
				cphi=parti->px/pt; sphi=parti->py/pt;
				cphitot+=cphi;
				sphitot+=sphi;
				cphi2tot+=cphi*cphi;
				sphi2tot+=sphi*sphi;
				pxtot+=parti->px;
				pytot+=parti->py;
				px2tot+=parti->px*parti->px;
				py2tot+=parti->py*parti->py;
				v2+=cphi*cphi-sphi*sphi;
				v2prime+=parti->px*parti->px-parti->py*parti->py;

				if(CALCGARRAYS){
					jpartmax=ipart-1;
					if(ipart%2==1) jpartmax=ipart-2;
					for(jpart=0;jpart<=jpartmax;jpart++){
						partj=&partH5[jpart];
						pt2=partj->px*partj->px+partj->py*partj->py;
						m=partj->mass;
						mt=sqrt(m*m+pt2);
						pz=mt*sinh(partj->rapidity);
						pt=sqrt(pt2);
						pmag=sqrt(pt2+pz*pz);
						etaj=atanh(pz/pmag);
						if((!STAR_ACCEPTANCE) || (pt<2000 && fabs(etaj)<1.0 && pt>200.0)){
							cphij=partj->px/pt;
							sphij=partj->py/pt;
							ieta=lrint(floor(fabs((eta-etaj)/deleta)));
							if(ieta<NETA){
								gss[ieta]+=(cphi*cphij-sphi*sphij);
								denom[ieta]+=1;
								gssprime[ieta]+=(parti->px*partj->px-parti->py*partj->py);
							}
						}
					}
				}
				if(ipart%2==1){
					jpart=ipart-1;
					partj=&partH5[jpart];
					pt2=partj->px*partj->px+partj->py*partj->py;
					m=partj->mass;
					mt=sqrt(m*m+pt2);
					pz=mt*sinh(partj->rapidity);
					pt=sqrt(pt2);
					pmag=sqrt(pt2+pz*pz);
					etaj=atanh(pz/pmag);
					if((!STAR_ACCEPTANCE) || (pt<2000 && fabs(etaj)<1.0 && pt>200.0)){
						cphij=partj->px/pt;
						sphij=partj->py/pt;
				//printf("sphij=%g, cphij=%g, pt=%g, m=%g\n",sphij,cphij,pt,m);
						ieta=lrint(floor(fabs((eta-etaj)/deleta)));
						if(ieta<NETA){
							gb[ieta]+=cphi*cphij-sphi*sphij;
							gbprime[ieta]+=parti->px*partj->px-parti->py*partj->py;
						}
						gammab+=cphi*cphij-sphi*sphij;
						gammabprime+=parti->px*partj->px-parti->py*partj->py;
					}
				}
			}
		}
		naccept+=dnaccept;
		npairs+=dnaccept*(dnaccept-1);
		gammass+=(cphitot*cphitot-cphi2tot-sphitot*sphitot+sphi2tot);
		gammassprime+=(pxtot*pxtot-px2tot-pytot*pytot+py2tot);

		if((10*(ievent+1))%nevents==0) printf("---------- Finished %g percent ------------\n",100.0*double(ievent+1)/double(nevents));
	}
	delete h5file;

	ptbar=ptbar/double(naccept);
	pt2bar=pt2bar/double(naccept);
	gammass=(gammass-2.0*gammab)*double(nsample)/(double(npairs));
	gammab=4.0*gammab*double(nsample)/(double(npairs));
	v2=v2/double(naccept);
	v2prime=v2prime/(pt2bar*double(naccept));
	gammassprime=(gammassprime-2.0*gammabprime)*double(nsample)/(double(npairs)*pt2bar);
	gammabprime=4.0*gammabprime*double(nsample)/(double(npairs)*pt2bar);
	for(ieta=0;ieta<NETA;ieta++){
		css[ieta]=gss[ieta]*double(nsample)/double(denom[ieta]);
		cssprime[ieta]=gssprime[ieta]*double(nsample)/(pt2bar*double(denom[ieta]));
		cb[ieta]=2.0*gb[ieta]*double(nsample)/double(denom[ieta]);
		cbprime[ieta]=2.0*gbprime[ieta]*double(nsample)/(pt2bar*double(denom[ieta]));

		gss[ieta]=2.0*gss[ieta]*double(nsample)/(deleta*double(npairs));
		gssprime[ieta]=2.0*gssprime[ieta]*double(nsample)/(pt2bar*deleta*double(npairs));
		gb[ieta]=4.0*gb[ieta]*double(nsample)/(deleta*double(npairs));
		gbprime[ieta]=4.0*gbprime[ieta]*double(nsample)/(pt2bar*deleta*double(npairs));
	}


	string filename;
	if(STAR_ACCEPTANCE)	filename="results/"+qualifier+"_STAR.dat";
	else filename="results/"+qualifier+".dat";
	FILE *fptr=fopen(filename.c_str(),"w");

	// Now fixing for 1/3 neutrals
	double M=2.0*double(naccept)/(3.0*double(nsample*nevents));
	printf("Total Multiplicity within acceptance = M=%g\n",M);
	gammab*=1.5;
	gammabprime*=1.5;
	for(ieta=0;ieta<NETA;ieta++){
		cb[ieta]*=1.5;
		cbprime[ieta]*=1.5;
		gb[ieta]*=1.5;
		gbprime[ieta]*=1.5;
	}

	printf("___________ MULTIPLICITY for |eta|<1.0 =%g _____________\n",double(MULT)/double(nevents*nsample));
	printf("ptbar=%g, nevents=%d, naccept=%g, npairs=%g\n",ptbar,nevents,double(naccept),double(npairs));
	printf("v2=%g, v2prime=%g\n",v2,v2prime);
	printf("gamma_ss=%g, gamma_ss'=%g\n",gammass,gammassprime);
	printf("gamma_b=%g, gamma_b'=%g\n",gammab,gammabprime);
	printf("gss*M/v2=%g, gss'*M/v2'=%g, gb*M/v2=%g, gb'*M/v2'=%g, p-sumrule=%g, p-sumrule'=%g=?1.0\n",
		gammass*M/v2,gammassprime*M/v2prime,gammab*M/v2,gammabprime*M/v2prime,
		(1.5*gammass+0.5*gammab)*M/v2,(1.5*gammassprime+0.5*gammabprime)*M/v2prime);
	printf("# eta   gamma_ss    gamma_b   gamma_ss'   gamma_b' |     c_ss       c_b       c_ss'      c_b'\n");

	fprintf(fptr,"#___________ MULTIPLICITY for |eta|<1.0 =%g _____________\n",double(MULT)/double(nevents*nsample));
	fprintf(fptr,"#ptbar=%g, nevents=%d, naccept=%g, npairs=%g\n",ptbar,nevents,double(naccept),double(npairs));
	fprintf(fptr,"#v2=%g, v2prime=%g\n",v2,v2prime);
	fprintf(fptr,"#gamma_ss=%g, gamma_ss'=%g\n",gammass,gammassprime);
	fprintf(fptr,"#gamma_b=%g, gamma_b'=%g\n",gammab,gammabprime);
	fprintf(fptr,"#gss*M/v2=%g, gss'*M/v2'=%g, gb*M/v2=%g, gb'*M/v2'=%g, p-sumrule=%g, p-sumrule'=%g=?1.0\n",
		gammass*M/v2,gammassprime*M/v2prime,gammab*M/v2,gammabprime*M/v2prime,
		(1.5*gammass+0.5*gammab)*M/v2,(1.5*gammassprime+0.5*gammabprime)*M/v2prime);
	fprintf(fptr,"# eta   gamma_ss    gamma_b   gamma_ss'   gamma_b' |     c_ss       c_b       c_ss'      c_b'\n");


	for(ieta=0;ieta<NETA;ieta++){
		printf("%4.2f | %10.3e %10.3e %10.3e %10.3e | %10.3e %10.3e %10.3e %10.3e\n",(0.5+ieta)*deleta,gss[ieta],gb[ieta],gssprime[ieta],gbprime[ieta],
			css[ieta],cb[ieta],cssprime[ieta],cbprime[ieta]);
		fprintf(fptr,"%4.2f   %10.3e %10.3e %10.3e %10.3e   %10.3e %10.3e %10.3e %10.3e\n",(0.5+ieta)*deleta,gss[ieta],gb[ieta],gssprime[ieta],gbprime[ieta],
			css[ieta],cb[ieta],cssprime[ieta],cbprime[ieta]);
	}
	fclose(fptr);
}

#endif