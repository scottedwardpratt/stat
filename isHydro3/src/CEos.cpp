/*
 *  CEosIdeal.cpp
 *  Created by Joshua Vredevoogd on 3/4/09.
 */

#include "CEos.h"
#include "hydroDef.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string.h>

//constructor - read from EOS interpolation
CEos::CEos(parameterMap* pM) {
 lastAccess=10;
 pMap = pM;
 mLatEos  = parameter::getB(*pMap,"EQOFST_LATEOS",true);
 mSVRatio = parameter::getD(*pMap,"HYDRO_SVRATIO",0.0);
 mBVRatio = parameter::getD(*pMap,"HYDRO_BVRATIO",0.0);

 if (mLatEos) {
 
	 const int size=2400;
  
	 temp = new double[size];
	 tA   = new double[size];
	 ed   = new double[size]; 
	 pr   = new double[size]; 
	 sd   = new double[size]; 
	 cs2  = new double[size];
	 
	 if (parameter::getS(*pMap,"EQOFST_LATDATA","null") == string("null")){
		 printf("\n** CEos requires parameterMap variable EQOFST_LATDATA if EQOFST_LATEOS == true **\nAborting...\n\n");
		 exit(1);
	 }
		 
	 std::ifstream latFile(parameter::getS(*pMap,"EQOFST_LATDATA","none").c_str());
	 
	 string hrgFileName = parameter::getS(*pMap,"EQOFST_HRGDATA","none");
	 std::ifstream hrgFile;
	 if (hrgFileName != string("none"))
		 hrgFile.open(hrgFileName.c_str());
  
	 double lowT  = parameter::getD(*pMap,"EQOFST_LOW_CROSS_T",0.15);
	 double highT = parameter::getD(*pMap,"EQOFST_HIGH_CROSS_T",0.2);
  
  double lT[size], lE[size], lP[size], lS[size];
  double hT[size], hE[size], hP[size], hS[size];
  int lSize=0;
  int hSize=0;
  
  char test;
  latFile >> test;
  latFile.ignore(1000,'\n');
  
	 if (test == 'e') { // EOS from Pasi
		 double junk;
		 int i;
		 for(i=0; !latFile.eof() && i<size; i++) {
			 latFile >> lE[i] >> lP[i] >> lS[i] >> junk;
			 lT[i] = (lE[i] + lP[i])/lS[i];
		 }
		 lSize = i-1;
	
		 if (hrgFileName == string("none")) {
			 printf(" *** No HRG EOS data provided for merging - continuity not ensured... ***\n *** Proceeding with lattice data only... ***\n\n");
			 aSize = lSize;
			 for (int j=0;j<aSize;j++) {
				 temp[j] = lT[j]; 
				 sd[j] = lS[j];
				 ed[j] = lE[j];
				 pr[j] = lP[j];
			 } 
			 return;
		 }
	  
	  
	lowT  += (lT[1] - lT[0])/10.;
	highT += (lT[1] - lT[0])/10.;
	
  	for (int j=0;j<8;j++) 
		hrgFile.ignore(1000,'\n');
		
	for(i=0; !hrgFile.eof() && i<size; i++) {
		hrgFile >> hT[i] >> hS[i] >> hP[i] >> hE[i];
		hT[i] /= 1000.; hP[i] /= 1000.; hE[i] /= 1000.;
		if (hT[i] < 0.02) i--;
	}
	hSize = i-1;	

	int hlowTIndex=0, llowTIndex=0, highTIndex=0;
	while (hT[hlowTIndex] < lowT ) hlowTIndex++;
	while (lT[llowTIndex] < lowT ) llowTIndex++; 
	while (lT[highTIndex] < highT) highTIndex++;
//	highTIndex--;
	
	aSize=0;
	for (int j=0; j<=hlowTIndex; j++) {
		sd[j] = hS[j];
		temp[j] = hT[j];
		aSize++;
	}

	for (int j=1;;j++) {
		temp[aSize] = lT[llowTIndex+j];
		if (temp[aSize] >= highT) break;
		double w = (tanh( tan((PI/2.)*( 2.*(highT-temp[aSize])/(highT-lowT) - 1.)))+1.)/2. ;
		sd[aSize] = w*hS[hlowTIndex+j] + (1.-w)*lS[llowTIndex+j];
		aSize++;
	}
	aSize--;
	
	for (int j=highTIndex;j<lSize;j++) {
		sd[aSize] = lS[j];
		temp[aSize] = lT[j];
		aSize++;
	}
	aSize--;

	for (int j=0;j<=aSize;j++) 
		if (j<= hlowTIndex)
			pr[j] = hP[j];
		else 
			pr[j] = pr[j-1] + 0.5*(temp[j]-temp[j-1])*(sd[j]+sd[j-1]);
	
	for (int j=0;j<=aSize;j++) 
		ed[j] = temp[j]*sd[j] - pr[j];
		
	for (int j=1;j<aSize;j++)
		cs2[j] = (sd[j]/temp[j])*(temp[j+1]-temp[j-1])/(sd[j+1]-sd[j-1]);
	cs2[0] = 2*cs2[1] - cs2[2];
	cs2[aSize] = 2*cs2[aSize-1] - cs2[aSize-2];
	
//	for (int j=0;j<=aSize;j++)
//		printf("%0.6g %0.6g %0.6g %0.6g %0.6g\n",temp[j],ed[j],sd[j],pr[j],cs2[j]);
  } 
  else { // EOS from Ron
    double junk;
    int i;
    for(i=0; !latFile.eof() && i<1200; i++) {
      latFile >> temp[i] >>  tA[i] >> ed[i] >> pr[i] >> sd[i] >> junk >> cs2[i];
	  tA[i]   *= pow(temp[i],4.)/pow(JHBARC,3.);
	  ed[i]   *= pow(temp[i],4.)/pow(JHBARC,3.);
	  pr[i]   *= pow(temp[i],4.)/pow(JHBARC,3.);
	  sd[i]   *= pow(temp[i],3.)/pow(JHBARC,3.);
  }	
    aSize = i-1;
  }
 }
}
CEos::~CEos() {
  delete [] temp;
  delete [] tA;
  delete [] ed;
  delete [] pr;
  delete [] sd;
  delete [] cs2;
}
double CEos::getCs2(double e){
  if (!mLatEos) return (1./3.);
  else {
    if (e < ed[0]) {
      lastAccess = 0;
//      return cs2[0];
	  return pr[0]/ed[0];
    }

    if (e > ed[lastAccess]) {
      for (;lastAccess<aSize;lastAccess++) if (ed[lastAccess] > e) break;
      return ((ed[lastAccess]-e)*cs2[lastAccess-1] 
			+ (e-ed[lastAccess-1])*cs2[lastAccess])/(ed[lastAccess] - ed[lastAccess-1]);
	}
	else {
	  for (;lastAccess>0;lastAccess--)     if (ed[lastAccess] < e) break;
      return ((ed[lastAccess+1]-e)*cs2[lastAccess] 
			+ (e-ed[lastAccess])*cs2[lastAccess+1])/(ed[lastAccess+1] - ed[lastAccess]);
	}
  }
}
double CEos::getP(double e) {
  if (!mLatEos) return (e/3.);
  else {
    if (e < ed[0]) {
      lastAccess = 0;
	  return cs2[0]*e;
//      return pr[0] * e/ed[0];
    }

    if (e > ed[lastAccess]) {
      for (;lastAccess<aSize;lastAccess++) if (ed[lastAccess] > e) break;
      return ((ed[lastAccess]-e)*pr[lastAccess-1] 
			+ (e-ed[lastAccess-1])*pr[lastAccess])/(ed[lastAccess] - ed[lastAccess-1]);
	}
	else {
	  for (;lastAccess>0;lastAccess--) if (ed[lastAccess] < e) break;
      return ((ed[lastAccess+1]-e)*pr[lastAccess] 
			+ (e-ed[lastAccess])*pr[lastAccess+1])/(ed[lastAccess+1] - ed[lastAccess]);
	}
  }
}
double CEos::getS(double e)  {
  if (!mLatEos) return (e + getP(e))/getT(e);
  else {
    if (e < ed[0]) {
      lastAccess = 0;
      return sd[0] * pow(e/ed[0],0.75);
    }

    if (e > ed[lastAccess]) {
      for (;lastAccess<aSize;lastAccess++) if (ed[lastAccess] > e) break;
      return ((ed[lastAccess]-e)*sd[lastAccess-1] 
			+ (e-ed[lastAccess-1])*sd[lastAccess])/(ed[lastAccess] - ed[lastAccess-1]);
	}
	else {
	  for (;lastAccess>0;lastAccess--) if (ed[lastAccess] < e) break;
      return ((ed[lastAccess+1]-e)*sd[lastAccess] 
			+ (e-ed[lastAccess])*sd[lastAccess+1])/(ed[lastAccess+1] - ed[lastAccess]);
	}
  }
}
double CEos::getTIS(double e) {
//	return 2.;
	return 3.*mSVRatio*JHBARC/getT(e);
}
double CEos::getTISB(double e) {
	return (3.*mBVRatio*JHBARC)/getT(e);
}
double CEos::getSV(double e) {
//	return mSVRatio*e;
	return (JHBARC*mSVRatio*getS(e));
}
double CEos::getBV(double e) {
//	return mBVRatio*e;
	return JHBARC*mBVRatio*getS(e);
}
double CEos::getT(double e) {
  if (!mLatEos) return 0.5179033247*pow(e,0.25)*pow(JHBARC,0.75);  // e = (16 + 21/2 * NDOF)* PI^2/30 * T^4/HBARC^3 for NDOF=2.5
  else {
    if (e < ed[0]) {
      lastAccess = 0;
      return temp[0] * pow(e/ed[0],0.25);
    }

    if (e > ed[lastAccess]) {
	  for (;lastAccess<aSize;lastAccess++) 	if (ed[lastAccess] > e) break;
      return ((ed[lastAccess]-e)*temp[lastAccess-1] 
			+ (e-ed[lastAccess-1])*temp[lastAccess])/(ed[lastAccess] - ed[lastAccess-1]);
	}
	else {
	  for (;lastAccess>0;lastAccess--) 	if (ed[lastAccess] < e) break;
      return ((ed[lastAccess+1]-e)*temp[lastAccess] 
			+ (e-ed[lastAccess])*temp[lastAccess+1])/(ed[lastAccess+1] - ed[lastAccess]);
	}
  }
}
double CEos::getSigmaA(double e) {
	return getISAlpha(e)/sqrt(getS(e));
}
double CEos::getSigmaB(double e) {
	return getISGamma(e)/sqrt(getS(e));
}
double CEos::getISAlpha(double e){
  if	  (ISSCALE == 'j') return sqrt(getT(e)*getSV(e)/getTIS(e));
  else if (ISSCALE == 'e') return e;
  else if (ISSCALE == 's') return getS(e);
  else return 1.;
}
double CEos::getISGamma(double e){
  if	  (ISSCALE == 'j') return sqrt( 0.8*getP(e)*getT(e)*getS(e));
  else if (ISSCALE == 'e') return e;
  else if (ISSCALE == 's') return getS(e);
  else return 1.;
}
double CEos::getDISAlphaDE(double e) {
	if      (ISSCALE == 'j') return (5./(8.*e));
	else if (ISSCALE == 'e') return 1.;
	else if (ISSCALE == 's') return (1./(4.*e));
	else if (ISSCALE == 'c') return 0.;
//  return ( getISAlpha(1.001*e) - getISAlpha(0.999*e))/(0.002*e);
}
double CEos::getDISGammaDE(double e) {
  return ( getISGamma(1.001*e) - getISGamma(0.999*e))/(0.002*e);
}
double CEos::getISBeta(double e) {
	return (e+getP(e))/(3.*getISAlpha(e));
}
double CEos::getTA(double e) {
  if (mLatEos) {
    if (e < ed[0]) {
      lastAccess = 0;
      return tA[0];
    }

    if (e > ed[lastAccess]) {
      for (;lastAccess<aSize;lastAccess++) if (ed[lastAccess] > e) break;
      return ((ed[lastAccess]-e)*tA[lastAccess-1] 
			+ (e-ed[lastAccess-1])*tA[lastAccess])/(ed[lastAccess] - ed[lastAccess-1]);
	}
	else {
	  for (;lastAccess>0;lastAccess--) if (ed[lastAccess] < e) break;
      return ((ed[lastAccess+1]-e)*tA[lastAccess] 
			+ (e-ed[lastAccess])*tA[lastAccess+1])/(ed[lastAccess+1] - ed[lastAccess]);
	}
  }
  else return 0.;
}
double CEos::getEGivenT(double T) {
  if (!mLatEos) return 1811.507976 * pow(T,4);
  else {
    int i=1;
	for (;i<aSize;i++) if (temp[i] >= T) break;
	if (T==temp[i]) return ed[i];
	return ((temp[i] - T)*ed[i-1] + (T-temp[i-1])*ed[i])/(temp[i] - temp[i-1]);
  }
}

bool CEos::mLatEos;
double CEos::mSVRatio, CEos::mBVRatio;
