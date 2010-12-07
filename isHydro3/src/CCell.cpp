/*  CCell.cpp
 *  isHydro
 *  Created by Joshua Vredevoogd on 2/11/09.
 */

#include "CCell.h"
#include "CEos.h"
#include <iostream>
#include <fstream>
#include "coral.h"

CCell::CCell(parameterMap* pM){
	paramFill(pM);
}

CCell::CCell(double t0,double n0, double x0, double y0) { 
  if (mLinT) x[0] = t0;
  else if (mLogT) x[0] = log(t0/mT0);
  else if (mLogSinhT) x[0] = log( sinh(t0));
//  else printf("whoa!!!\n");
  x[3] = n0;
  x[1] = x0;
  x[2] = y0;
  
//  printf("%f %f %f %f\n",x[0],x[1],x[2],x[3]);
}

CCell::CCell(CCell* cell) {
	for (int i=0;i<4;i++) x[i] = cell->x[i];
	for (int i=0;i<11;i++) s[i] = cell->s[i];
	active = cell->active;
}  

CCell::~CCell() {
	//noop
}
  
void CCell::paramFill(parameterMap* pM) {
  pMap = pM;
  
  mDebug  = parameter::getB(*pMap,"HYDRO_DEBUG",false);
  mSVTrim = parameter::getB(*pMap,"HYDRO_SVTRIM",false);
  mViscNS = parameter::getB(*pMap,"HYDRO_VISCNS",false);
  mPureBjorken = parameter::getB(*pMap,"HYDRO_PURE_BJORKEN",false);
  mBjorken = parameter::getB(*pMap,"HYDRO_BJORKEN",true);
  mLinT = parameter::getB(*pMap,"HYDRO_LINT",false);
  mLogT = parameter::getB(*pMap,"HYDRO_LOGT",false);
  mLogSinhT = parameter::getB(*pMap,"HYDRO_LOGSINHT",true);
  mISVort = parameter::getB(*pMap,"HYDRO_IS_VORT",false);
  mISMax = parameter::getB(*pMap,"HYDRO_IS_MAX",false);

  mT0 = parameter::getD(*pMap,"HYDRO_T0",1.0);
  mSVRatio = parameter::getD(*pMap,"HYDRO_SVRATIO",0.0);
  mBVRatio = parameter::getD(*pMap,"HYDRO_BVRATIO",0.0);
  mISAMax  = parameter::getD(*pMap,"HYDRO_IS_AMAX",0.0);
  mISBMax  = parameter::getD(*pMap,"HYDRO_IS_BMAX",0.0);
  
//  dx[0] = parameter::getD(*pMap,"DT");
  dx[1] = parameter::getD(*pMap,"HYDRO_DX",0.1);
  dx[2] = parameter::getD(*pMap,"HYDRO_DY",0.1);
  dx[3] = parameter::getD(*pMap,"HYDRO_DN",0.1);
}

void CCell::calcDeriv() {
	
	for (int j=0;j<10;j++)
		for (int i=0;i<3;i++)
			if (neighbors[i][1]->getActive())
				dS[j][i] = (neighbors[i][1]->s[j+1] - neighbors[i][0]->s[j+1])/(2.*dx[i+1]);
			else 
				dS[j][i] = (3.*s[j+1] - 4.*neighbors[i][0]->s[j+1] 
							+ neighbors[i][0]->neighbors[i][0]->s[j+1]) / (2.*dx[i+1]);
	
	for (int i=0;i<3;i++)
		dS[3][i] *= getE();
	
	if (mPureBjorken)
		for (int j=0;j<10;j++)
			dS[j][2] = 0.;
	else if (mBjorken)
		for (int j=0;j<10;j++)
			dS[j][2] /= getTau();
	
	if (mBjorken) dS[2][2] += s[0]/getTau();
	
	for (int j=0;j<4;j++) 
		for (int i=0;i<4;i++)  
			dULocal[j][i] = getSpLocDS(j,i);
	
	if (ISSCALE == 'c') for (int i=0;i<3;i++) dTISS[i] = 0.;
	if (ISSCALE == 'e') for (int i=0;i<3;i++) dTISS[i] = dS[3][i];
	if (ISSCALE == 'j') for (int i=0;i<3;i++) dTISS[i] = (5./8.)*dS[3][i]*exp(-s[4]);  //FIXME	
	
	for (int j=0;j<3;j++)
		for (int k=0;k<3;k++) 
			for (int i=0;i<3;i++)
				DiPixyLocal[i][j][k] = getDiPixyLocal(i,j+1,k+1);
	
	for (int i=0;i<3;i++) {
		DPixy[i] = 0.;
		for (int j=0;j<3;j++) 
			DPixy[i] += DiPixyLocal[j][j][i];
	}
}

// Fills the matrix inverted in forward
double CCell::getM(int i, int j) { 
  double value=0.;
  
  if (j==10) {
    if (i==0) {
	  for (int m=0;m<3;m++) {
	    value += s[m+1] * dS[3][m];
		value += (getE() + getP()) * dS[m][m];
		for (int n=0;n<3;n++) 
		  value += Pixy[m][n]*dULocal[m+1][n+1];
	  }
	  return (-value);
    } 
	else if (i<4) {  
	  value = getCs2() * dS[3][i-1];
 	  value += (getE() + getP()) * dULocal[0][i];
	  value += DPixy[i-1];
	  for (int m=0;m<3;m++) {
		value += Pixy[i-1][m] * dULocal[0][m+1];
		value += s[i]*s[m+1]/(1+s[0]) * getCs2() * dS[3][m];
		for (int n=0;n<3;n++) 
		  value += s[m+1]*s[n+1]/(1+s[0]) * DiPixyLocal[m][n][i-1];
	  }
	  return (-value);
    } 
	else {					// i>3
	  if ( (i<9 && mSVRatio == 0.) || (i==9 && mBVRatio == 0.) ) return 0.;
	  if (mViscNS) return 0.;
	  
	  if (i<9) {
		value -= s[i+1]/getTIS();
		for (int m=0;m<3;m++) value -= s[m+1]*dS[i][m];
		
		if (ISRESCALE) {
		  for (int m=0;m<3;m++) value -= (4./3.)*s[i+1]*dS[m][m];
		}
	  }
	  
	  if(!mISMax) {
		if (i==4) {  //IS equations
		  value -= getISBeta() * ( dULocal[1][1] - dULocal[2][2]);
		  return value;
		} 
		else if (i==5) {
		  value -= getISBeta()/ROOT3 * ( dULocal[1][1] + dULocal[2][2] - 2.*dULocal[3][3]);
		  return value;
		} 
		else if (i==6) {
		  value -= getISBeta() * (dULocal[2][1] + dULocal[1][2]);
		  return value;
		} 
		else if (i==7) {
		  value -= getISBeta() * (dULocal[1][3] + dULocal[3][1]);
		  return value;
		} 
		else if (i==8) {
		  value -= getISBeta() * (dULocal[2][3] + dULocal[3][2]);
			return value;
		} 
		else {  //i==9 FIXME
		  return 0.;
		}
	  }
	  else {  //mISMax FIXME
		return 0.;
	  }
	}
  } //end if(j==10)

  if (i>3 && j>2 && i!=j) return 0.;

  if (i==0) {  
	if (j<3) return ((getE()+getP())*s[j+1]/s[0]) + up[j] - (uup*s[j+1]/(s[0]*(1+s[0])));
	else if (j==3) return s[0];
	else return 0.;
  } 
  else if (i<4) { // d_k T_ik = 0
    if (j<3) {
	  if (j==i-1) value += s[0]*(getE() + getP());

	  value -= s[i]*s[j+1]/(1+s[0]) * (getE() + getP());
	  value += s[0]*Pixy[i-1][j];
	  value -= up[i-1]*s[j+1]/(1+s[0]);
	  return value;
	}
	else if (j==3) {
	  value += s[i]*getCs2();
	  value += getDISAlphaDE()*up[i-1]/getISAlpha();
	  return value;
	}
	else if (j==4){
	  if      (i==1) value += getISAlpha()*s[1];
	  else if (i==2) value -= getISAlpha()*s[2];
	  return value;
	} 
	else if (j==5){
	  if     (i==1) value += getISAlpha()*s[1]/ROOT3;
	  else if(i==2) value += getISAlpha()*s[2]/ROOT3;
	  else          value -= getISAlpha()*2*s[3]/ROOT3;
	  return value;
	}
	else if (j==6){
	  if     (i==1) value += getISAlpha()*s[2];
	  else if(i==2) value += getISAlpha()*s[1];
	  return value;
	}
	else if (j==7){
	  if     (i==1) value += getISAlpha()*s[3];
	  else if(i==3) value += getISAlpha()*s[1];
	  return value;
	}
	else if (j==8){
	  if     (i==3) value += getISAlpha()*s[2];
	  else if(i==2) value += getISAlpha()*s[3];
	  return value;
	}
	else return 0.;
  } 
  else {
	//IS equations
    if (mViscNS) {
	  if (i==j) return 1.;
	  else return 0.;
	}
	if (i==j) return s[0]; 
	if( (i<9 && mSVRatio == 0.) || (i==9 && mBVRatio == 0.) ) return 0.;
	
	if (ISRESCALE) value = (4./3.)*s[i+1]*s[j+1]/s[0];
	if (i==4) {
	  value -= getISBeta()*(s[1]*s[1] - s[2]*s[2]) * s[j+1]/(s[0]*(1+s[0]));
	  
	  if (j==0) {
		value += getISBeta()*s[1];
		if (mISVort) value += s[7]*s[2] + 0.5*s[8]*s[3];
		return value;
	  } 
	  else if (j==1) {
	    value -= getISBeta()*s[2];
		if (mISVort) value += s[7]*s[1] + 0.5*s[9]*s[3];
		return value;
	  } 
	  else if (j==2) {
		if (mISVort) value += 0.5 * ( s[8]*s[1] - s[9]*s[2]);
		return value;
	  } 
	}
	else if(i==5) {
	  value -= (getISBeta()/ROOT3) *(s[1]*s[1] + s[2]*s[2] - 2*s[3]*s[3]) * s[j+1]/(s[0]*(1.+s[0]));

	  if(j==0) {
	    value += (getISBeta()/ROOT3) * s[1];
  		if (mISVort) value -= s[8] * (ROOT3/2.) * s[3];
		return value;
	  } 
	  else if (j==1) { 
		value += (getISBeta()/ROOT3) * s[2];
		if (mISVort) value -= (ROOT3/2.) * s[9]*s[3];
		return value;
	  } 
	  else if (j==2) {
		value -= (getISBeta()/ROOT3) * 2.*s[3];
		if (mISVort) value += (ROOT3/2.) * ( s[8]*s[1] + s[9]*s[2]);
		return value;
	  }
	}  
	else if(i==6) {
	  value -= getISBeta()*(2.*s[1]*s[2]) * s[j+1]/(s[0]*(1+s[0]));
	  if(j==0) {
		value += getISBeta() * s[2];
		if (mISVort) 
		  value += s[5]*s[2] - 0.5*s[9]*s[3];
		return value;
	  } 
	  else if (j==1) {
	    value += getISBeta() * s[1];
		if (mISVort) 
		  value = -(s[5]*s[1] + 0.5*s[8]*s[3]);
		return value;
	  } 
	  else if (j==2) {
		if (mISVort) 
		  value += 0.5*(s[8]*s[2] + s[9]*s[1]);
		return value;	
	  }
	} //i6
	else if(i==7) {
	  value -= getISBeta()*(2*s[1]*s[3]) * s[j+1]/(s[0]*(1+s[0]));
	  if(j==0) {
	    value += getISBeta() * s[3];
		if (mISVort) value += 0.5 * (s[5] + ROOT3*s[6]) * s[3] - 0.5*s[9]*s[2];
		return value;
	  } 
	  else if (j==1) {
		if (mISVort) value += 0.5 * (s[9]*s[1] + s[7]*s[3]);
		return value;
	  } 
	  else if (j==2) {
		value += getISBeta() * s[1];
		if (mISVort) value += -0.5 * ((s[5] + ROOT3*s[6])*s[1] + s[7]*s[2]);
		return value;	
	  }
	} //i7
	else if(i==8) {
	  value -= getISBeta()*(2*s[2]*s[3]) * s[j+1]/(s[0]*(1+s[0]));
	  if(j==0) {
		if (mISVort) value += 0.5 * (s[7]*s[3] + s[8]*s[2]);
		return value;
	  } 
	  else if (j==1) {
	    value += getISBeta() * s[3];
		if (mISVort) value += 0.5 * ( s[3]*(ROOT3*s[6] - s[5]) - s[8]*s[1]);
		return value;
	  } 
	  else if (j==2) {
		value += getISBeta() * s[2];
		if (mISVort) value += -0.5 * ( s[2]*(ROOT3*s[6] - s[5]) + s[7]*s[1]);
		return value;
	  }
	} //i8
	else if(i==9) {
	  return 0.;
	} //i9
  } // end IS EQs
  
  std::cout << "Unknown case of getM! (" << i << ',' << j << ")" << std::endl;
  return 0.;
}

void CCell::update() {
  fillEosVar();
  if (mViscNS) {
	s[5] = (-SV/alphaIS) * (dS[0][0] - dS[1][1]);	
	s[6] = (-SV/alphaIS) * (dS[0][0] + dS[1][1] - 2.*dS[2][2])/ROOT3;
	s[7] = (-SV/alphaIS) * (dS[0][1] + dS[1][0]);
	s[8] = (-SV/alphaIS) * (dS[0][2] + dS[2][0]);
	s[9] = (-SV/alphaIS) * (dS[1][2] + dS[2][1]);
  }

  for (int i=0;i<3;i++) 
	for (int j=0;j<3;j++) 
	  Pixy[i][j] = getPixy(i+1,j+1);

  uup = 0.;
  for (int i=0;i<3;i++) {
    up[i] = 0.;
    for (int j=0;j<3;j++) 
	  up[i] += s[j+1]*Pixy[j][i];	  
	uup += s[i+1] * up[i];
  }
  calcDeriv();
  fillM();
}

void CCell::fillEosVar(){ 
  cS2	= eos->getCs2( getE());
  P		= eos->getP( getE());
  S		= eos->getS( getE());
  T		= eos->getT(getE());
  SV	= eos->getSV(getE());
  BV	= eos->getBV(getE());
  
  tIS	= eos->getTIS(getE()); 
  tISB  = eos->getTISB(getE()); 
  
  alphaIS= eos->getISAlpha(getE()); 
  gammaIS= eos->getISGamma( getE());
  sigmaA = eos->getSigmaA( getE());
  sigmaB = eos->getSigmaB( getE());
  aIS	 = eos->getISA( getE());
  bIS	 = eos->getISB( getE());	
  dAlphaISDE = eos->getDISAlphaDE(getE());
  dGammaISDE = eos->getDISAlphaDE(getE());
  betaIS = SV/(alphaIS*tIS);
  
  if (mSVTrim) {
    SV  /= (1. + exp( (sqrt(x[1]*x[1]+x[2]*x[2]) - 10.)/0.6));
	tIS	/= (1. + exp( (sqrt(x[1]*x[1]+x[2]*x[2]) - 10.)/0.6));  
  }
}

void CCell::update(int ii) {
  uup  = (s[0]*s[0] - 1.) * s[10];
  uup += (s[1]*s[1] - s[2]*s[2]) * s[5];
  uup += (ROOT3/3.) * (s[1]*s[1] + s[2]*s[2] - 2.*s[3]*s[3]) * s[6];
  uup += s[1]*s[2]*s[7];
  uup += s[1]*s[3]*s[8];
  uup += s[2]*s[3]*s[9];

  up[0] = s[1] * ( s[10] + (ROOT3/3.) * s[6]);
  up[0] += s[2] * s[7];
  up[0] += s[3] * s[8];

  up[1] = s[1] * s[7];
  up[1] += s[2] * ( s[10] + (ROOT3/3.) * s[6]);
  up[1] += s[3] * s[9];

  up[2] = s[1] * s[8];
  up[2] += s[2] * s[9];
  up[2] += s[3] * ( s[10] + 2.*(ROOT3/3.)*s[6]);

  for (int i=0;i<3;i++) 
	for (int j=0;j<3;j++) Pixy[i][j] = getPixy(i+1,j+1);

//  calcDeriv();
  for (int i=0;i<10;i++)
    for (int j=0;j<3;j++) dS[i][j] = 0.;
  fillM();
  
  fillEosVar();
}

void CCell::initNS() {
	update();
	
	s[5] = - (getSV()/getISAlpha()) * (dULocal[1][1] - dULocal[2][2]);
	s[6] = - getSV()/(ROOT3*getISAlpha()) * (dULocal[1][1] + dULocal[2][2] - 2.*dULocal[3][3]);
	s[7] = - (getSV()/getISAlpha()) * (dULocal[1][2] + dULocal[2][1]);
	s[8] = - (getSV()/getISAlpha()) * (dULocal[1][3] + dULocal[3][1]);
	s[9] = - (getSV()/getISAlpha()) * (dULocal[3][2] + dULocal[2][3]);	
}

// turn this cell off; if impossible return false
// FIXME :assumes octant=true
void CCell::deactivate() {
	active = false;
	for (int i=0;i<3;i++) 
		if (neighbors[i][0]->x[i+1] == 0.)
			neighbors[i][0]->deactivate();
}

double CCell::getTxy(int x, int y) {
#ifdef DEBUG
  if (x<1 || x>3 || y<1 || y>3) {
	std::cout << "Unknown case of CCell::getTxy(" << x << ',' << y << "):" << std::endl;
	return 0.;
  }
#endif
  // off-diagonals
  if (x>y) {
    if (y==1) {
      if (x==2) return s[7]*getISAlpha();
      else return s[8]*getISAlpha();
    } else return s[9]*getISAlpha();
  } 
  else if(y>x) {
    if (x==1) {
      if (y==2) return s[7]*getISAlpha();
      else return s[8]*getISAlpha();
    } else return s[9]*getISAlpha();
  } 

  // diagonals
  if (x==1) return (s[10]*getISGamma() + getP() + s[5]*getISAlpha() + s[6]*getISAlpha()/ROOT3);
  else if (x==2) return (s[10]*getISGamma() + getP() - s[5]*getISAlpha() + s[6]*getISAlpha()/ROOT3);
  else return (s[10]*getISGamma() + getP() - (2./ROOT3) * s[6]*getISAlpha());
}

double CCell::getPixy(int x, int y) {
#ifdef DEBUG
  if (x<1 || x>3 || y<1 || y>3) {
	std::cout << "Unknwon case of CCell::getPixy(" << x << ',' << y << "):" << std::endl;
	return 0.;
  }
#endif
  // off-diagonals
  if (x>y) {
    if (y==1) {
      if (x==2) return s[7]*getISAlpha();
      else return s[8]*getISAlpha();
    } else return s[9]*getISAlpha();
  } 
  else if(y>x) {
    if (x==1) {
      if (y==2) return s[7]*getISAlpha();
      else return s[8]*getISAlpha();
    } else return s[9]*getISAlpha();
  } 

  // diagonals
  if (x==1) return (s[10]*getISGamma() + s[5]*getISAlpha() + s[6]*getISAlpha()/ROOT3);
  else if (x==2) return (s[10]*getISGamma() - s[5]*getISAlpha() + s[6]*getISAlpha()/ROOT3);
  else return (s[10]*getISGamma() - (2./ROOT3) * s[6]*getISAlpha());
}

double CCell::getTxyMesh(int x, int y) {
#ifdef DEBUG
  if (x<0 || x>3 || y<0 || y>3) {
	std::cout << "Unknwon case of CCell::getTxyMesh(" << x << ',' << y << "):" << std::endl;
	return 0.;
  }
#endif
  double value = 0.;
  if (x==y && y==0) {
    value = s[0]*s[0]*getE();
	for (int i=1;i<4;i++) 
	  for (int j=1;j<4;j++) 
	    value += s[i]*s[j]*getTxy(i,j);
	return value;
  } 
  else if (x==0) {
    value = s[0]*s[y]*getE();
	for (int i=1;i<4;i++) {
	  value += getTxy(i,y) * s[i];
	  for (int j=1;j<4;j++) value += (s[y]/(1.+s[0]))* s[i]*s[j]* getTxy(i,j);
	}
  } 
  else if (y==0) {
    value = s[0]*s[x]*getE();
	for (int i=1;i<4;i++) {
	  value += getTxy(i,x) * s[i];
	  for (int j=1;j<4;j++) value += (s[x]/(1.+s[0]))* s[i]*s[j]* getTxy(i,j);
	}
  } 
  else {
	value = s[x]*s[y]*getE() + getTxy(x,y);
	for (int i=1;i<4;i++) {
	  value += s[x]*s[i]/(1+s[0])*getTxy(i,y);
	  value += s[i]*s[y]/(1+s[0])*getTxy(x,i);
	  for (int j=1;j<4;j++) 
	    value += s[x]*s[y]* s[i]*s[j]/((1.+s[0])*(1.+s[0])) * getTxy(i,j);
	}
  }
}

double CCell::getPixyMesh(int x, int y) {
#ifdef DEBUG
  if (x<0 || x>3 || y<0 || y>3) {
	std::cout << "Unknwon case of CCell::getTxyMesh(" << x << ',' << y << "):" << std::endl;
	return 0.;
  }
#endif
  double value = 0.;
  if (x==y && y==0) {
	for (int i=1;i<4;i++) 
	  for (int j=1;j<4;j++) 
	    value += s[i]*s[j]*getPixy(i,j);
	return value;
  } 
  else if (x==0) {
	for (int i=1;i<4;i++) {
	  value += getPixy(i,y) * s[i];
	  for (int j=1;j<4;j++) value += (s[y]/(1.+s[0]))* s[i]*s[j]* getPixy(i,j);
	}
  } 
  else if (y==0) {
	for (int i=1;i<4;i++) {
	  value += getPixy(i,x) * s[i];
	  for (int j=1;j<4;j++) value += (s[x]/(1.+s[0]))* s[i]*s[j]* getPixy(i,j);
	}
  } 
  else {
	value = getPixy(x,y);
	for (int i=1;i<4;i++) {
	  value += s[x]*s[i]/(1+s[0])*getPixy(i,y);
	  value += s[i]*s[y]/(1+s[0])*getPixy(x,i);
	  for (int j=1;j<4;j++) 
	    value += s[x]*s[y]* s[i]*s[j]/((1.+s[0])*(1.+s[0])) * getPixy(i,j);
	}
  }
}

double CCell::getTxyEos(int x, int y) {
#ifdef DEBUG
  if (x<1 || x>3 || y<1 || y>3) {
	std::cout << "Unknwon case of CCell::getTxyEos(" << x << ',' << y << "):" << std::endl;
	return 0.;
  }
#endif
  // off-diagonals
  if (x>y) {
    if (y==1) {
      if (x==2) return s[7]*getISAlphaEos();
      else return s[8]*getISAlphaEos();
    } else return s[9]*getISAlphaEos();
  } 
  else if(y>x) {
    if (x==1) {
      if (y==2) return s[7]*getISAlphaEos();
      else return s[8]*getISAlphaEos();
    } else return s[9]*getISAlphaEos();
  } 

  // diagonals
  if (x==1) return (s[10]*getISGammaEos() + getPEos() + (s[5] + s[6]/ROOT3)*getISAlphaEos());
  else if (x==2) return (s[10]*getISGammaEos() + getPEos() + (s[6]/ROOT3- s[5])*getISAlphaEos());
  else return (s[10]*getISGammaEos() + getPEos() - (2./ROOT3) * s[6]*getISAlphaEos());
}

double CCell::getPixyEos(int x, int y) {
#ifdef DEBUG
  if (x<1 || x>3 || y<1 || y>3) {
	std::cout << "Unknwon case of CCell::getPixy(" << x << ',' << y << "):" << std::endl;
	return 0.;
  }
#endif
  // off-diagonals
  if (x>y) {
    if (y==1) {
      if (x==2) return s[7]*getISAlphaEos();
      else return s[8]*getISAlphaEos();
    } else return s[9]*getISAlphaEos();
  } 
  else if(y>x) {
    if (x==1) {
      if (y==2) return s[7]*getISAlphaEos();
      else return s[8]*getISAlphaEos();
    } else return s[9]*getISAlphaEos();
  } 

  // diagonals
  if (x==1) return (s[10]*getISGammaEos() + (s[5] + s[6]/ROOT3)*getISAlphaEos());
  else if (x==2) return (s[10]*getISGammaEos() + (s[6]/ROOT3 - s[5])*getISAlphaEos());
  else return (s[10]*getISGammaEos() - (2./ROOT3)*s[6]*getISAlphaEos());
}

double CCell::getTxyMeshEos(int x, int y) {
#ifdef DEBUG
  if (x<0 || x>3 || y<0 || y>3) {
	std::cout << "Unknwon case of CCell::getTxyMesh(" << x << ',' << y << "):" << std::endl;
	return 0.;
  }
#endif
  double value = 0.;
  if (x==y && y==0) {
    value = s[0]*s[0]*getE();
	for (int i=1;i<4;i++) 
	  for (int j=1;j<4;j++) 
	    value += s[i]*s[j]*getTxyEos(i,j);
  } 
  else if (x==0) {
    value = s[0]*s[y]*getE();
	for (int i=1;i<4;i++) {
	  value += getTxyEos(i,y) * s[i];
	  for (int j=1;j<4;j++) value += (s[y]/(1.+s[0]))* s[i]*s[j]* getTxyEos(i,j);
	}
  } 
  else if (y==0) {
    value = s[0]*s[x]*getE();
	for (int i=1;i<4;i++) {
	  value += getTxyEos(i,x) * s[i];
	  for (int j=1;j<4;j++) value += (s[x]/(1.+s[0]))* s[i]*s[j]* getTxyEos(i,j);
	}
  } 
  else {
	value = getTxyEos(x,y) + s[x]*s[y]*getE();
	for (int i=1;i<4;i++) {
	  value += s[x]*s[i]*getTxyEos(i,y) / (1.+s[0]);
	  value += getTxyEos(x,i)*s[i]*s[y] / (1.+s[0]);
	  for (int j=1;j<4;j++) 
	    value += s[x]*s[y]* s[i]*s[j]*getTxyEos(i,j)/( (1.+s[0])*(1.+s[0]));
	}
  }
  return value;
}

double CCell::getPixyMeshEos(int x, int y) {
#ifdef DEBUG
  if (x<0 || x>3 || y<0 || y>3) {
	std::cout << "Unknwon case of CCell::getTxyMesh(" << x << ',' << y << "):" << std::endl;
	return 0.;
  }
#endif
  double value = 0.;
  if (x==y && y==0) {
	for (int i=1;i<4;i++) 
	  for (int j=1;j<4;j++) 
	    value += s[i]*s[j]*getPixyEos(i,j);
  } 
  else if (x==0) {
	for (int i=1;i<4;i++) {
	  value += getPixyEos(i,y) * s[i];
	  for (int j=1;j<4;j++) 
	    value += (s[y]/(1.+s[0]))* s[i]*s[j]* getPixyEos(i,j);
	}
  } 
  else if (y==0) {
	for (int i=1;i<4;i++) {
	  value += getPixyEos(i,x) * s[i];
	  for (int j=1;j<4;j++) 
	    value += (s[x]/(1.+s[0]))* s[i]*s[j]* getPixyEos(i,j);
	}
  } 
  else {
	value = getPixyEos(x,y);
	for (int i=1;i<4;i++) {
	  value += s[x]*s[i]*getPixyEos(i,y) / (1.+s[0]);
	  value += getPixyEos(x,i) * s[i]*s[y] / (1.+s[0]);
	  for (int j=1;j<4;j++) 
	    value += s[x]*s[y]* s[i]*s[j]*getPixyEos(i,j)/( (1.+s[0])*(1.+s[0]));
	}
  }
  return value;
}

double CCell::getPixyNS(int x,int y) {
#ifdef DEBUG
  if (x<1 || x>3 || y<1 || y>3) {
	std::cout << "Unknwon case of CCell::getPixyNS(" << x << ',' << y << "):" << std::endl;
	return 0.;
  }
#endif
	double value = dULocal[x-1][y-1] + dULocal[y-1][x-1];
	value += s[x]*dUdT[y-1] + s[y]*dUdT[x-1];
	for (int i=0;i<3;i++) {
	  value -= 2.*(s[x]*s[y]*s[i+1]/(s[0]*(1+s[0]))) * dUdT[i];
	  if (x==y) value -= (2./3.) * ( dS[i][i] + (s[i+1]/s[0])*dUdT[i] ) ;
	}
	return -getSV()*value;
}

double CCell::getPixyNSMesh(int x,int y) {
	#ifdef DEBUG
  if (x<0 || x>3 || y<0 || y>3) {
	std::cout << "Unknwon case of CCell::getPixyNSMesh(" << x << ',' << y << "):" << std::endl;
	return 0.;
  }
#endif
  double value = 0.;
  if (x==y && y==0) {
	for (int i=1;i<4;i++) 
	  for (int j=1;j<4;j++) 
	    value += s[i]*s[j]*getPixyNS(i,j);
	return value;
  } 
  else if (x==0) {
	for (int i=1;i<4;i++) {
	  value += getPixyNS(i,y) * s[i];
	  for (int j=1;j<4;j++) 
	    value += (s[y]/(1.+s[0]))* s[i]*s[j]* getPixyNS(i,j);
	}
  } 
  else if (y==0) {
	for (int i=1;i<4;i++) {
	  value += getPixyNS(i,x) * s[i];
	  for (int j=1;j<4;j++) 
	    value += (s[x]/(1.+s[0]))* s[i]*s[j]* getPixyNS(i,j);
	}
  } 
  else {
	value = getPixyNS(x,y);
	for (int i=1;i<4;i++) {
	  value += s[x]*s[i]/(1+s[0])*getPixyNS(i,y);
	  value += s[i]*s[y]/(1+s[0])*getPixyNS(x,i);
	  for (int j=1;j<4;j++) 
	    value += s[x]*s[y]* s[i]*s[j]/((1.+s[0])*(1.+s[0])) * getPixyNS(i,j);
	}
  }
  return value;
}

double CCell::getDTxy(int i) {
#ifdef DEBUG 
	if (i<1 || i>3) {
//		printf("Unknown case of Ccell::getDTxy()\n");
		std::cout << "Unknown case of CCell::getDTxy(" << i << "):" << std::endl;
		return 0.;
	}
#endif
	double value = getISGamma()*dS[9][i-1] + s[10]*dTISS[i-1] + getCs2()*dS[3][i-1];  // trace elements
	if (i==1) {
		value += getISAlpha() * (dS[4][0] + dS[5][0]/ROOT3 + dS[6][1] + dS[7][2]);
		value += (s[5] + s[6]/ROOT3)*dTISS[0] + s[7]*dTISS[1] + s[8]*dTISS[2];
		return value;
	} 
	else if (i==2) {
		value += getISAlpha() * (dS[6][0] + dS[5][1]/ROOT3 - dS[4][1] + dS[8][2]);
		value += s[7]*dTISS[0] + (s[6]/ROOT3 - s[5])*dTISS[1] + s[9]*dTISS[2];
		return value;
	} 
	else if (i==3) {
		value += getISAlpha() * (dS[7][0] + dS[8][1] - (2./ROOT3)*dS[5][2]);
		value += s[8]*dTISS[0] + s[9]*dTISS[1] - (2./ROOT3)*s[6]*dTISS[2];
		return value;
	} 
	else return 0.;
}

double CCell::getDiTxy(int i, int j, int k) {
#ifdef DEBUG
	if (i<1 || i>3 || j<1 || j>3 || k<1 || k>3) {
//		printf("Unknown case of Ccell::getDTxy(%d,%d,%d)\n",i,j,k);
		std::cout << "Unknown case of CCell::getDTxy(" << i << ',' << j << ',' << k << ',' << ')' << std::endl;
		return 0.;
	}
#endif
	double value = 0.;
	
	if (j==k) {
		value += getISGamma()*dS[9][i-1] + s[10]*dTISS[i-1] + getCs2()*dS[3][i-1];
		if (j==1) 
		  value += getISAlpha()*dS[4][i-1] + s[5]*dTISS[i-1] + (getISAlpha()/ROOT3) *dS[5][i-1] + s[6]*dTISS[i-1]/ROOT3;
		else if (j==2) 
		  value +=  getISAlpha()*dS[5][i-1]/ROOT3 + s[6]*dTISS[i-1]/ROOT3 - getISGamma()*dS[4][i-1] - s[5]*dTISS[i-1];
		else if (j==3) 
		  value -= getISAlpha()*(2./ROOT3)*dS[5][i-1] + s[6]*dTISS[i-1]*(2./ROOT3);
		return value;
	} 
	else { 
		if ((j==1 && k==2) || (j==2 && k==1)) return (getISAlpha()*dS[6][i-1] + s[7]*dTISS[i-1]);
		if ((j==1 && k==3) || (j==3 && k==1)) return (getISAlpha()*dS[7][i-1] + s[8]*dTISS[i-1]);
		if ((j==2 && k==3) || (j==3 && k==2)) return (getISAlpha()*dS[8][i-1] + s[9]*dTISS[i-1]);
	}
	return 0.; // can't....
}

double CCell::getS() {
//	return eos->getS( getE());
    double value = eos->getS( getE());
	if (mSVRatio != 0.) 
      for (int i=5;i<10;i++)  
		value -= 0.5*s[i]*s[i];
//	    value *= exp(- 0.5 * s[i] * s[i]);
	if ( mBVRatio != 0.) 
	  value -= 0.5*s[10]*s[10];
//	  value *= exp(- 0.5 * s[10] * s[10]);
	return value;
}

// returns the spatial parts ONLY of
// the local velocity gradients required for IS
double CCell::getSpLocDS(int j,int i) {
#ifdef DEBUG
	if (i<0 || i>3 || j<0 || j>3) {
		std::cout << "Unknown case of CCell::getSpLocDS(" << i << ',' << j << std::endl; 
		return 0.;
	}
#endif 
	double value=0.;
	
	if (j==0 && i==0) return 0.;
	else if (i==0) {
	  for (int m=1;m<4;m++) {
	    value += s[m]*dS[j-1][m-1];
		for (int n=1;n<4;n++)
		  value -= s[j]*s[m]*s[n]/(s[0]*(1+s[0])) * dS[m-1][n-1];
	  }
	}
	else if (j==0) {
	  for (int m=1;m<4;m++) {
	    value += s[m]*dS[i-1][m-1];
		for (int n=1;n<4;n++)
		  value -= s[i]*s[m]*s[n]/(s[0]*(1+s[0])) * dS[m-1][n-1];
	  }
	}
	else {
	  value = dS[j-1][i-1];
	  for (int k=1;k<4;k++) {
		value += s[i]/(1.+s[0]) * s[k] * dS[j-1][k-1];
		value -= s[j]*s[k]/(s[0]*(1+s[0])) * dS[k-1][i-1];
		for (int l=1;l<4;l++) 
			value -= s[i]*s[j]*s[k]*s[l]/(s[0]*(1.+s[0])*(1.+s[0])) * dS[l-1][k-1];
	  }
/*
	value -= s[0]/getTau() * s[i]*s[j]*s[3]*s[3]/(s[0]*(1+s[0])*(1+s[0]));
	if (i==3) 
	  value -= s[0]/getTau() *s[j]*s[3]/(s[0]*(1+s[0]));
	if (j==3) {
	  value += s[0]/getTau() * s[i]*s[3]/(1+s[0]);
	  if (i==3) value += s[0]/getTau();
	}
	return value;
*/
	}
	
	return value;
}

double CCell::getDiPixyLocal(int i,int j,int k) {
#ifdef DEBUG
  if (j<1 || j>3 || k<1 || k>3 || i<0 || i>2) {
	std::cout << "Unknwon case of CCell::getDiPixyLocal(" << i << ',' << j << ',' << k << "):" << std::endl;
	return 0.;
  }
#endif
/*
  double value = DiPixy[i][j][k];
  
  value -= s[j] * DiPixy[i][0][k];
  value -= s[k] * DiPixy[i][j][0];
  value += s[j] * s[k] * DiPixy[i][0][0];

  for (int m=1;m<4;m++) {
    value += s[j]*s[m]/(1+s[0]) * DiPixy[i][m][k];
	value += s[k]*s[m]/(1+s[0]) * DiPixy[i][j][m];
	value -= 2.*s[j]*s[k]*s[m]/((1+s[0])*(1+s[0])) * DiPixy[i][0][m];
	for (int n=1;n<4;n++)
	  value += s[j]*s[k]*s[m]*s[n]/((1+s[0])*(1+s[0])) * DiPixy[i][m][n];
  }
*/
  if (k>j) {
    if (j==1) {
      if (k==2) return dS[6][i]*getISAlpha() + s[7]*dTISS[i];
      else      return dS[7][i]*getISAlpha() + s[8]*dTISS[i];
    } else      return dS[8][i]*getISAlpha() + s[9]*dTISS[i];
  } 
  else if(j>k) {
    if (k==1) {
      if (j==2) return dS[6][i]*getISAlpha() + s[7]*dTISS[i];
      else      return dS[7][i]*getISAlpha() + s[8]*dTISS[i];
    } else      return dS[8][i]*getISAlpha() + s[9]*dTISS[i];
  } 

  // diagonals
  if (k==1)       return dS[9][i]*getISGamma() + dS[4][i]*getISAlpha() + dS[5][i]*getISAlpha()/ROOT3 + dTISS[i]*Pixy[0][0]/getISAlpha();
  else if (k==2) return dS[9][i]*getISGamma() - dS[4][i]*getISAlpha() + dS[5][i]*getISAlpha()/ROOT3 + dTISS[i]*Pixy[1][1]/getISAlpha();
  else           return dS[9][i]*getISGamma() - (2./ROOT3)*dS[5][i]*getISAlpha() + dTISS[i]*Pixy[2][2]/getISAlpha();
}

void CCell::setU(double x0, double y0, double z0) {
  s[1] = x0;
  s[2] = y0;
  s[3] = z0;
  s[0] = sqrt(1. + x0*x0 + y0*y0 + z0*z0);
}

void CCell::setS(int i, double v) {
#ifdef DEBUG
	if (i<0 || i>10) {
//		printf("Unknown case of setS(%d,%0.6g)!\n",i,v);
		std::cout << "Unknown case of setS("<< i << ',' << v << ")!" << std::endl;
		return;
	}
#endif
	if (i!=4) s[i]=v;
	else s[i] = log(v);
	if (i<4) s[0] = sqrt(1. + s[1]*s[1]+s[2]*s[2]+s[3]*s[3]);
}

void CCell::fillM() {
	// weird due to row re-ordering...
	for (int i=0;i<11;i++) M[3][i] = getM(0,i);

	for (int i=1;i<10;i++) 
		for (int j=0;j<11;j++) {
			if (i>3) M[i][j] = getM(i,j);
			else	 M[i-1][j] = getM(i,j);
		}
}

// Integrates this cell forward in time
// NOTE: update() must be called for all cells before forward
void CCell::forward(CCell* mCell) {
  // fill the mixed partials and M_ij
  update();

  if (ISSCALE == 'c') {
  // remove tau_ij dep from cons eqs                                                                             
  for (int k=0;k<4;k++) {
    for (int j=4;j<10;j++) {
      for (int i=0; i<j;i++) 
	    M[k][i] -= (M[k][j]/M[j][j])*M[j][i];
	  M[k][10] -= (M[k][j]/M[j][j])*M[j][10];
    }
  }

  // finish LT form
  for (int k=3;k>0;k--)
    for (int j=k-1;j>=0;j--) {
	  M[j][10] -= (M[j][k]/M[k][k])*M[k][10];
      for (int i=0;i<k;i++)  M[j][i] -= (M[j][k]/M[k][k])*M[k][i];
	  M[j][k] = 0.;
    }

  // diagonalize M_ii, calculate change to vector side (M_j_11)                                 
  for (int i=0;i<3;i++)
    for (int j=i+1;j<11;j++) 
      M[j][10] -= (M[j][i]/M[i][i])*M[i][10];
  }
  else {
  // remove tau_ij dep from cons eqs                                                                             
  for (int k=0;k<4;k++) 
    for (int j=4;j<10;j++) {
      for (int i=0; i<j;i++) M[k][i] -= (M[k][j]/M[j][j])*M[j][i];
	  M[k][10] -= (M[k][j]/M[j][j])*M[j][10];
    }

  // finish LT form
  for (int k=3;k>0;k--)
    for (int j=k-1;j>=0;j--) {
	  M[j][10] -= (M[j][k]/M[k][k])*M[k][10];
      for (int i=0;i<=k;i++)  M[j][i] -= (M[j][k]/M[k][k])*M[k][i];
    }

  // diagonal M_ii is now as it will be, calculate change to vector side (M_j_11)                                
  for (int i=0;i<=3;i++)
    for (int j=i+1;j<11;j++) 
      M[j][10] -= (M[j][i]/M[i][i])*M[i][10];
}

  for (int i=0;i<3;i++) {
	mCell->dUdT[i] = (M[i][10]/M[i][i]);
    if (mLinT)			mCell->setS(i+1, s[i+1] + dx[0] * M[i][10]/M[i][i]);
	else if (mLogT)		mCell->setS(i+1, s[i+1] + mT0*exp(x[0])*dx[0]*M[i][10]/M[i][i]);
	else if (mLogSinhT)	mCell->setS(i+1, s[i+1] + tanh(getTau())*dx[0]*M[i][10]/M[i][i]);
  }
  
  if (x[1] == 6. && x[2] == 1. && false) {
    double v1 = 0.;	double v2 = 0.;
	
	for (int i=0;i<3;i++) {
	  v1 += dULocal[i][i] + s[i+1]*dUdT[i];
	  for (int j=0;j<3;j++) 
	    v1 += s[i+1]*s[i+1]*s[j+1] * ( 1./(1.+s[0]) - 1./s[0]) * dUdT[j];
	  v2 += dS[i][i] + (s[i+1]/s[0])*dUdT[i];
	}
	printf("%0.6g %0.6g\n",v1,v2);
  }
  
  if (mLinT) {
	mCell->setS(4, exp(s[4] + exp(-s[4]) * dx[0] * M[3][10]/M[3][3]));
	for (int i=4;i<10;i++) 
	  mCell->setS(i+1, s[i+1] + dx[0] * M[i][10]/M[i][i]);
  }
  else if (mLogT) {
	mCell->setS(4, exp(s[4] + exp(-s[4]+x[0]) * (dx[0]*mT0) * M[3][10]/M[3][3]));
	for (int i=4;i<10;i++) 
	  mCell->setS(i+1, s[i+1] +  mT0*exp(x[0]) * dx[0] * M[i][10]/M[i][i]);	  
  }
  else if (mLogSinhT){
	mCell->setS(4, exp(s[4] + exp(-s[4]) * dx[0] * tanh(getTau()) * M[3][10]/M[3][3]));
	for (int i=4;i<10;i++)
	  mCell->setS(i+1, s[i+1] + tanh(getTau())*dx[0]*M[i][10]/M[i][i]);
  }
  
  if (mLinT) mCell->setTau( x[0] + dx[0]);
  else if (mLogT) mCell->setTau(mT0*exp(x[0] + dx[0]));
  else if (mLogSinhT) mCell->setTau(asinh(exp(x[0]+dx[0])));
}

// Integrates this cell forward in time
void CCell::forward(CCell* onCell, CCell* offCell) {
  // fill the mixed partials and M_ij
  offCell->update();

  if (ISSCALE == 'c') {
  // remove tau_ij dep from cons eqs                                                                             
  for (int k=0;k<4;k++) {
    for (int j=4;j<10;j++) {
      for (int i=0; i<j;i++) M[k][i] -= (M[k][j]/M[j][j])*M[j][i];
	  M[k][10] -= (M[k][j]/M[j][j])*M[j][10];
    }
  }

  // finish LT form
  for (int k=3;k>0;k--)
    for (int j=k-1;j>=0;j--) {
	  M[j][10] -= (M[j][k]/M[k][k])*M[k][10];
      for (int i=0;i<k;i++)  M[j][i] -= (M[j][k]/M[k][k])*M[k][i];
	  M[j][k] = 0.;
    }

  // diagonal M_ii is now as it will be, calculate change to vector side (M_j_11)                                
  for (int i=0;i<3;i++)
    for (int j=i+1;j<11;j++) 
      M[j][10] -= (M[j][i]/M[i][i])*M[i][10];
	  
  } 
  else {
  // remove tau_ij dep from cons eqs                                                                             
  for (int k=0;k<4;k++) 
    for (int j=4;j<10;j++) {
      for (int i=0; i<j;i++) M[k][i] -= (M[k][j]/M[j][j])*M[j][i];
	  M[k][10] -= (M[k][j]/M[j][j])*M[j][10];
    }

  // finish LT form
  for (int k=3;k>0;k--)
    for (int j=k-1;j>=0;j--) {
	  M[j][10] -= (M[j][k]/M[k][k])*M[k][10];
      for (int i=0;i<=k;i++)  M[j][i] -= (M[j][k]/M[k][k])*M[k][i];
    }

  // diagonal M_ii is now as it will be, calculate change to vector side (M_j_11)                                
  for (int i=0;i<=3;i++)
    for (int j=i+1;j<11;j++) 
      M[j][10] -= (M[j][i]/M[i][i])*M[i][10];
}
  
  // note - off-diagonal elements HAVE NOT been set to zero, but HAVE been taken in to account

  for (int i=0;i<3;i++) {
    onCell->dUdT[i] = M[i][10]/M[i][i];
    if (mLinT)			onCell->setS(i+1, s[i+1] + dx[0]*M[i][10]/M[i][i]);
	else if (mLogT)		onCell->setS(i+1, s[i+1] + mT0*exp(offCell->x[0]) * dx[0] * M[i][10]/M[i][i]);
	else if (mLogSinhT)	onCell->setS(i+1, s[i+1] + tanh(offCell->getTau()) * dx[0] * M[i][10]/M[i][i]);
  }
  
  if (mLinT) {
	for (int i=4;i<10;i++) 
	    onCell->setS(i+1, s[i+1] + dx[0]*M[i][10]/M[i][i]);
    onCell->setS(4, exp( s[4] + exp(-offCell->s[4]) * dx[0] * M[3][10]/M[3][3]));
	onCell->setTau( x[0] + dx[0]);
  }
  else if (mLogT) {
    for (int i=4;i<10;i++) 
	    onCell->setS(i+1, s[i+1] + mT0*exp(offCell->x[0]) * dx[0] * M[i][10]/M[i][i]);
    onCell->setS(4, exp(s[4] + exp(-offCell->s[4]+offCell->x[0]) * (dx[0]*mT0) * M[3][10]/M[3][3]));
    onCell->setTau( mT0 * exp( x[0] + dx[0]));
  } 
  else if (mLogSinhT) {
	onCell->setS(4, exp(s[4] + dx[0] * exp(-offCell->s[4])*tanh(offCell->getTau())*M[3][10]/M[3][3]));
	for (int i=4;i<10;i++)
	    onCell->setS(i+1, s[i+1] + dx[0] * tanh(offCell->getTau()) * M[i][10]/M[i][i]);
	onCell->setTau(asinh(exp(x[0] + dx[0])));
  }

if (x[1] == 0. && false) {
  double mTemp = 0., mTemp2 =0.;
  printf("\nTest for cell: %0.6g %0.6g %0.6g %0.6g\n",getTau(),x[1],x[2],x[3]);
  
  mTemp += offCell->s[0] * M[3][10]/M[3][3];
  printf("1: %0.6g\n",mTemp);
  
  for (int i=0;i<3;i++) 
	mTemp2 = (4./3.)*offCell->getE()*offCell->s[i+1]/offCell->s[0] * onCell->dUdT[i];
	
  mTemp += mTemp2;
  printf("2: %0.6g %0.6g\n",mTemp,mTemp2);
  
  for (int i=0;i<3;i++)
	mTemp2 = up[i] * onCell->dUdT[i];
	
  mTemp += mTemp2;
  printf("3: %0.6g %0.6g\n",mTemp,mTemp2);
  
  for (int i=0;i<3;i++)
	mTemp2 = - uup*offCell->s[i+1]*onCell->dUdT[i]/(offCell->s[0]*(1+offCell->s[0]));
	
  mTemp += mTemp2;
  printf("4: %0.6g %0.6g\n",mTemp,mTemp2);
  
  for (int i=0;i<3;i++)
	mTemp2 = offCell->s[i+1]*dS[3][i];
	
  mTemp += mTemp2;
  printf("5: %0.6g %0.6g\n",mTemp,mTemp2);

  for (int i=0;i<3;i++)
	mTemp2 = (4./3.)*offCell->getE() * dS[i][i];
	
  mTemp += mTemp2;
  printf("6: %0.6g %0.6g\n",mTemp,mTemp2);
  	
  for (int i=0;i<3;i++)
	for (int j=0;j<3;j++)
	  mTemp2 = Pixy[i][j] * dULocal[i+1][j+1];
  mTemp += mTemp2;
  printf("7: %0.6g %0.6g\n",mTemp,mTemp2);

  }
}

void CCell::forward(CCell* k0, CCell* k1, CCell* k2, CCell* k3){
  // since the matrix for RK4 has contributions from 4 different
  // time steps, we'll compile it here from the usual update()
  // of the four cells to be used.
  double mM[10][11];

  k0->update();
  for (int i=0;i<10;i++) {
	for (int j=0;j<10;j++) mM[i][j] = (M[i][j]/6.);
	if (mLogT)
	  if (i!=3) mM[i][10] = M[i][10]*exp(k0->x[0])/6.;
	  else mM[i][10] = M[i][10]*exp(-k0->s[4] + k0->x[0])/6.;
	else if (mLogSinhT)
	  if (i!=3) mM[i][10] = M[i][10]*tanh(k0->getTau())/6.;
	  else mM[i][10] = M[i][10]*exp(-k0->s[4])*tanh(k0->getTau())/6.;
	else if (mLinT)
	  if (i!=3) mM[i][10] = M[i][10]/6.;
	  else mM[i][10] = M[i][10]*exp(-k0->s[4])/6.;
  } 
  
  k3->update();
  for (int i=0;i<10;i++) {
	for (int j=0;j<10;j++) mM[i][j] += (M[i][j]/6.);
	if (mLogT)
	  if (i!=3) mM[i][10] += (M[i][10]*exp(k3->x[0])/6.);
	  else mM[i][10] += (M[i][10]*exp(-k3->s[4] + k3->x[0])/6.);
	else if (mLogSinhT)
	  if (i!=3) mM[i][10] += (M[i][10]*tanh(k3->getTau())/6.);
	  else mM[i][10] += (M[i][10]*exp(-k3->s[4])*tanh(k3->getTau())/6.);
	else if (mLinT)
	  if (i!=3) mM[i][10] = M[i][10]/6.;
	  else mM[i][10] = M[i][10]*exp(-k3->s[4])/6.;
  }

  k1->update();
  for (int i=0;i<10;i++) {
	for (int j=0;j<10;j++) mM[i][j] += M[i][j]/3.;
	if (mLogT)
	  if (i!=3)  mM[i][10] += M[i][10]*exp(k1->x[0])/3.;
	  else mM[i][10] += M[i][10]*exp(-k1->s[4] + k1->x[0])/3.;
	else if (mLogSinhT)
	  if (i!=3)  mM[i][10] += M[i][10]*tanh(k1->getTau())/3.;
	  else mM[i][10] += M[i][10]*exp(-k1->s[4])*tanh(k1->getTau())/3.;
    else if (mLinT)
	  if (i!=3) mM[i][10] = M[i][10]/6.;
	  else mM[i][10] = M[i][10]*exp(-k1->s[4])/6.;
  }

  k2->update();
  for (int i=0;i<10;i++) {
	for (int j=0;j<10;j++) mM[i][j] += M[i][j]/3.;
	if (mLogT)
	  if (i!=3) mM[i][10] += M[i][10]*exp( k2->x[0])/3.;
	  else mM[i][10] += M[i][10]*exp(-k2->s[4] + k2->x[0])/3.;
	else if (mLogSinhT)
	  if (i!=3) mM[i][10] += M[i][10]*tanh(k2->getTau())/3.;
	  else mM[i][10] += M[i][10]*exp(-k2->s[4])*tanh(k2->getTau())/3.;
	else if (mLinT)
	  if (i!=3) mM[i][10] = M[i][10]/6.;
	  else mM[i][10] = M[i][10]*exp(-k2->s[4])/3.;
  }
  
  // remove tau_ij dep from cons eqs                                                                             
  for (int k=0;k<4;k++) {
    for (int j=4;j<10;j++) {
      for (int i=0; i<j;i++) mM[k][i] -= (mM[k][j]/mM[j][j])*mM[j][i];
	  mM[k][10] -= (mM[k][j]/mM[j][j])*mM[j][10];
    }
  }
  // finish LT form                                                                                              
  for (int k=3;k>0;k--)
    for (int j=k-1;j>=0;j--) {
      for (int i=0;i<k;i++)  mM[j][i] -= (mM[j][k]/mM[k][k])*mM[k][i];
	  mM[j][10] -= (mM[j][k]/mM[k][k])*mM[k][10];
	  mM[j][k] = 0.;
    }
  // diagonal M_ii is now as it will be, calculate change to vector side (M_j_10)                                
  for (int i=0;i<=3;i++)
    for (int j=i+1;j<11;j++) 
      mM[j][10] -= (mM[j][i]/mM[i][i])*mM[i][10];

  double mV[3];
  for (int i=0;i<3;i++) 
    if (mLogT) mV[i] = mT0 * dx[0] * mM[i][10]/mM[i][i];
	else      mV[i] = dx[0] * mM[i][10]/mM[i][i];

  //Now move it forward
  double uv = s[1]*mV[0] + s[2]*mV[1] + s[3]*mV[2];
  for (int i=0;i<3;i++)
    s[i+1] =  mV[i] + s[i+1]*sqrt(mV[0]*mV[0]+mV[1]*mV[1]+mV[2]*mV[2]) + (uv/(1.+s[0])) * s[i+1];

  //Now move it forward
  if (mLogT) 
    for (int i=4;i<10;i++) 
      s[i+1] += mT0 * dx[0] * mM[i][10]/mM[i][i];
  else if (mLogSinhT || mLinT)
    for (int i=4;i<10;i++)
	  s[i+1] += (dx[0]) * mM[i][10]/mM[i][i];

  s[0] = sqrt(1. + s[1]*s[1] + s[2]*s[2] + s[3]*s[3]);
  x[0] += dx[0];
}

void CCell::print() { 
	std::cout << "Cell (" << getTau() << ','
			  << x[1] << ',' << x[2] << ',' << x[3] << "):" << std::endl;

	std::cout << "Cell (" << dx[0] << ','
			  << dx[1] << ',' << dx[2] << ',' << dx[3] << "):" << std::endl;

	for (int i=0;i<11;i++){
		if (i!=4) std::cout << s[i] << '\t';
		else std::cout << exp(s[i]) << '\t';
		if (i==3 || i==5) std::cout << std::endl << '\t';
	}
	std::cout << std::endl;
	
	std::cout << std::endl << std::endl;
	std::cout << "sv = " << getSV() << std::endl;
	std::cout << "TIS = " << getTIS() << std::endl;
	std::cout << "ISAlpha = " << getISAlpha() << std::endl;
	std::cout << "T = " << getT() << std::endl;
	std::cout << "S = " << getS() << std::endl;
	std::cout << "dTISS:  ";
	for (int i=0;i<3;i++) std::cout << dTISS[i] << "  ";
	
	std::cout << std::endl << "dAlphaISDE:  " << dAlphaISDE << std::endl;
	std::cout << "cs2 = " << getCs2() << std::endl;

	for (int i=0;i<3;i++) {
	  for (int j=0;j<10;j++) 
	    if (fabs(dS[j][i]) > 0.) 
		  std::cout << "dS[" << j << "][" << i << "] = " << dS[j][i] << std::endl;
	std::cout << std::endl;
	}
	
	std::cout << std::endl;
	std::cout << "dULocal:  " << std::endl;
	for (int i=0;i<4;i++) {
	  for (int j=0;j<4;j++) std::cout << dULocal[i][j] << '\t';
	  std::cout << std::endl;
	}
	
	std::cout << "dUdT "<< dUdT[0] << " " << dUdT[1] << " " << dUdT[2] << std::endl;

	std::cout << "Pixy:" << std::endl;
	for (int i=0;i<3;i++) {
	  for (int j=0;j<3;j++) std::cout << Pixy[i][j] << '\t';
	  std::cout << std::endl;
	}
	
	std::cout << "u's:  " << uup << " " << up[0] << " " << up[1] << " " << up[2] << std::endl;
	std::cout << "DPixy "<< DPixy[0] << " " << DPixy[1] << " " << DPixy[2] << std::endl;
	std::cout << "DiPixyLocal:" << std::endl;
	for (int i=0;i<3;i++) {
	  for (int j=0;j<3;j++) {
		for (int k=0;k<3;k++) 
		  std::cout << DiPixyLocal[i][j][k] << " ";
		std::cout << std::endl;
		}
 	  std::cout << std::endl << std::endl;
	}
	return;
	std::cout << "DiPixy:" << std::endl;
	for (int i=0;i<3;i++) {
	  for (int j=0;j<4;j++) {
		for (int k=0;k<4;k++) std::cout << DiPixy[i][j][k] << " ";
		std::cout << std::endl;
		}
 	  std::cout << std::endl << std::endl;
	}
	return;
}

void CCell::print(int ii) { 
	if (ii==0){
		std::cout << "Cell (" << getTau() << ','
				  << x[1] << ',' << x[2] << ',' << x[3] << "):" << std::endl;

		for (int i=0;i<11;i++){
			if (i!=4) std::cout << s[i] << '\t';
			else std::cout << exp(s[i]) << '\t';
			if (i==3 || i==5) std::cout << std::endl << '\t';
		}
	std::cout << std::endl << std::endl;
	std::cout << "sv = " << getSV() << std::endl;
	std::cout << "TIS = " << getTIS() << std::endl;
	} 
}

// prints the matrix inverted in forward()                                                                       
// should only be called for debugging purposes...                                                               
void CCell::printM(){

  std::cout << "Matrix for Cell(" << getTau() << ',' << x[1] << ',' << x[2] << ',' << x[3] << ")" << std::endl;
 if (mSVRatio != 0. || mBVRatio != 0.)
    for (int i=0;i<10;i++) {
      for (int j=0;j<11;j++) 
	    if (j>3 && j!=10) std::cout << M[i][j] << '\t';
		else std::cout << M[i][j] << '\t';
	  std::cout << std::endl << std::endl;
    } 
  else
    for (int i=0;i<4;i++) {
	  for (int j=0;j<4;j++) std::cout << M[i][j] << '\t';
	  std::cout << M[i][10] << std::endl << std::endl;
    }
}

void CCell::printMm(double mM[10][11]) {
  std::cout << "mMatrix for Cell(" << getTau() << ',' << x[1] << ',' << x[2] << ',' << x[3] << ")" << std::endl;
  
  if (mSVRatio != 0. || mBVRatio != 0.) 
  for (int i=0;i<10;i++) {
    for (int j=0;j<11;j++)
	  std::cout << mM[i][j] << '\t' << '\t';
	  std::cout << std::endl << std::endl;
  }
  else 
  for (int i=0;i<4;i++) {
	  for (int j=0;j<4;j++) std::cout << mM[i][j] << '\t';
	  std::cout << mM[i][10] << std::endl << std::endl;
    }
}

//statics
double CCell::dx[4];// = {DT,DX,DY,DN};
double CCell::uup = 0.;
double CCell::up[3] = {0.,0.,0.};
double CCell::dS[10][3];
double CCell::dULocal[4][4];
double CCell::Pixy[3][3];
double CCell::DiPixy[3][4][4];
double CCell::DiPixyLocal[3][3][3];
double CCell::DPixy[3];
double CCell::M[10][11];
CEos*  CCell::eos = NULL;
double CCell::dTISS[3];
double CCell::cS2, CCell::P, CCell::S, CCell::tIS, CCell::tISB;
double CCell::SV, CCell::BV, CCell::T, CCell::sigmaA, CCell::sigmaB;
double CCell::alphaIS, CCell::gammaIS, CCell::betaIS, CCell::aIS, CCell::bIS;
double CCell::dAlphaISDE, CCell::dGammaISDE;
bool CCell::mDebug, CCell::mSVTrim, CCell::mViscNS, CCell::mPureBjorken, CCell::mBjorken;
bool CCell::mLinT, CCell::mLogT, CCell::mLogSinhT, CCell::mISVort, CCell::mISMax;
double CCell::mT0, CCell::mSVRatio, CCell::mBVRatio, CCell::mISAMax, CCell::mISBMax, CCell::mInitNS;
parameterMap* CCell::pMap;
