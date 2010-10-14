/*
 *  CEosIdeal.cpp
 *  Created by Joshua Vredevoogd on 3/4/09.
 */

#include "CEos.h"
#include "def.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string.h>

//constructor - read from EOS interpolation
CEos::CEos(parameterMap* pM) {
 lastAccess=0;
 pMap = pM;
 mLatEos  = parameter::getB(*pMap,"EQOFST_LATEOS",true);
 mSVRatio = parameter::getD(*pMap,"HYDRO_SVRATIO",0.0);
 mBVRatio = parameter::getD(*pMap,"HYDRO_BVRATIO",0.0);

 if (mLatEos) {
  temp = new double[1200];
  tA   = new double[1200];
  ed   = new double[1200]; 
  pr   = new double[1200]; 
  sd   = new double[1200]; 
  cs2  = new double[1200];

//  std::ifstream mFile("hrg000p4_0MeV_eos.dat");
  std::ifstream mFile(parameter::getS(*pMap,"EQOFST_LATDATA","").c_str());
  char test;
  mFile >> test;
  mFile.ignore(1000,'\n');
  
  if (test == 'e') { // EOS from Pasi
    double junk;
    int i;
    for(i=0; !mFile.eof() && i<1200; i++) {
	  mFile >> ed[i] >> pr[i] >> sd[i] >> junk;
      temp[i] = (ed[i] + pr[i])/sd[i];
	  tA[i] = (ed[i] - 3*pr[i])/pow(temp[i],4);
	}
	aSize = i-1;
	for (int j=1;j+1<aSize;j++)
	  cs2[j] = (pr[j+1] - pr[j-1])/(ed[j+1] - ed[j-1]);
	cs2[0] = 2*cs2[1] - cs2[2];
	cs2[aSize-1] = 2*cs2[aSize-2] - cs2[aSize-3];
  } 
  else { // EOS from Ron
    double junk;
    int i;
    for(i=0; !mFile.eof() && i<1200; i++) {
      mFile >> temp[i] >>  tA[i] >> ed[i] >> pr[i] >> sd[i] >> junk >> cs2[i];
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
      return cs2[0];
//	  return pr[0]/ed[0];
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
      for (;lastAccess<aSize;lastAccess++) if (ed[lastAccess] > e) break;
      return ((ed[lastAccess]-e)*temp[lastAccess-1] 
			+ (e-ed[lastAccess-1])*temp[lastAccess])/(ed[lastAccess] - ed[lastAccess-1]);
	}
	else {
	  for (;lastAccess>0;lastAccess--)     if (ed[lastAccess] < e) break;
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
