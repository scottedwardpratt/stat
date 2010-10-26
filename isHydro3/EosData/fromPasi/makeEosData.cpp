//  runHydro.cpp
//  Created by Joshua Vredevoogd on 8/30/10.

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>
#define PRECISION "%0.6g"

int main(int argc, char * const argv[]){

	double a1 = 4.654;
	double a2 = -879.;
	double a3 = 8081.;
	double a4 = -7039000.;

	double d2, d4, c1, c2, T0;
	int n1, n2;

	std::string mEos(argv[1]);

	if (mEos == "s95p") {
		d2 = 0.266;
		d4 = 0.002403;
		c1 = -2.809E-7;
		c2 = 6.073E-23;
		n1 = 10;
		n2 = 30;
		T0 = 0.1838;
	} 
	else if (mEos == "s95n") {
		d2 = 0.2654;
		d4 = 0.006563;
		c1 = -4.37E-5;
		c2 = 5.774E-6;
		n1 = 8;
		n2 = 9;
		T0 = 0.1718;
	} 
	else if (mEos == "s90f") {
		d2 = 0.2495;
		d4 = 0.01355;
		c1 = -3.237E-3;
		c2 = 1.439E-14;
		n1 = 5;
		n2 = 18;
		T0 = 0.170;
	}
	else{
		printf("\nunknown eos requested.... terminating....\n");
		exit(1);
	}

	double dT = 0.001;
	const int tSize = 501;
	double temp[tSize], inter[tSize], pres[tSize], ed[tSize], s[tSize], cs2[tSize];
	
	temp[0]  = 0.07;	
	inter[0] = a1*temp[0] + a2*pow(temp[0],3) + a3*pow(temp[0],4) + a4*pow(temp[0],10);
	pres[0]   = 0.1661;
	ed[0]     = inter[0] + 3*pres[0];
	s[0]      = (ed[0] + pres[0])/temp[0];

	std::string fOutName(argv[1]);
	fOutName.append(".dat");
	FILE* fOut = fopen(fOutName.c_str(),"w");
	fprintf(fOut,"T I E P S E/P CS2 (GeV^n)\n");

	for(int i=1;i<tSize;i++) {
		temp[i]  = temp[i-1] + dT;
		
		if (temp[i] > T0)
			inter[i] = d2/pow(temp[i],2.) + d4/pow(temp[i],4.) + c1/pow(temp[i],n1) + c2/pow(temp[i],n2);
		else 
			inter[i] = a1*temp[i] + a2*pow(temp[i],3) + a3*pow(temp[i],4) + a4*pow(temp[i],10);
			
		pres[i]  = pres[i-1] + 0.5* dT * ( inter[i-1]/temp[i-1] + inter[i]/temp[i]);
		ed[i]    = inter[i] + 3*pres[i];
		s[i]     = (ed[i] + pres[i])/temp[i];
		
		if (i>1)  cs2[i-1] =  (pres[i]*pow(temp[i],4) - pres[i-2]*pow(temp[i-2],4))
							/ (ed[i]*pow(temp[i],4) - ed[i-2]*pow(temp[i-2],4));
	}
	
// linear extrapolations
	cs2[0]   = 2*cs2[1] - cs2[2];
	cs2[tSize-1] = 2*cs2[tSize-2] - cs2[tSize-3];

/*
// low temp stuff
	double lowTemp[71], lowInter[71], lowPres[71], lowEd[71], lowS[71], lowCs2[71];
	lowTemp[70]  = temp[0];
	lowInter[70] = inter[0];
	lowPres[70]  = pres[0];
	lowEd[70]    = ed[0];
	lowS[70]     = s[0];
	lowCs2[70]   = cs2[0];
	
	for (int i=69;i>=0;i--) {
		lowTemp[i] = dT*i;// lowTemp[i+1] - dT;
		lowInter[i] = a1*lowTemp[i] + a2*pow(lowTemp[i],3) + a3*pow(lowTemp[i],4) + a4*pow(lowTemp[i],10);
		lowPres[i] = lowPres[i+1] - 0.5*dT* (lowInter[i]/lowTemp[i] + lowInter[i+1]*lowTemp[i+1]);
		lowEd[i] = lowInter[i] + 3*lowPres[i];
		lowS[i] = (lowEd[i]+lowPres[i])/lowTemp[i];
	}
	
	for (int i=1;i<70;i++) 
		lowCs2[i] =  (lowPres[i+1]*pow(lowTemp[i+1],4) - lowPres[i-1]*pow(lowTemp[i-1],4))
					/(  lowEd[i+1]*pow(lowTemp[i+1],4) -   lowEd[i-1]*pow(lowTemp[i-1],4));
	
	lowCs2[0] = 2*lowCs2[1] - lowCs2[2];
	
// print it	
	for (int i=40;i<=70;i++)
		fprintf(fOut,"%0.9g %0.9g %0.9g %0.9g %0.9g %0.9g %0.9g\n",
					lowTemp[i],lowInter[i],lowEd[i],lowPres[i],lowS[i],lowPres[i]/lowEd[i],lowCs2[i]);
*/					
	for (int i=0;i<tSize;i++) 
		fprintf(fOut,"%0.12g %0.12g %0.12g %0.12g %0.12g %0.12g %0.12g\n",
					temp[i],inter[i],ed[i],pres[i],s[i],pres[i]/ed[i],cs2[i]);
}