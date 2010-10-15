#ifndef __CEOSIDEAL_h_INCLUDE__
#define __CEOSIDEAL_h_INCLUDE__

/*
 *  CEosIdeal.h
 *  Created by Joshua Vredevoogd on 3/4/09.
 */
#include "coral.h"

class CEos{
  private:
    double *temp;
    double *tA;
	double *ed;
	double *pr;
	double *sd;
	double *cs2;
	int aSize;
	int lastAccess;

	parameterMap* pMap;

	static bool mLatEos;
	static double mSVRatio, mBVRatio;

  public:
    CEos(parameterMap* pMap);
	~CEos();
	double getCs2(double e);
	double getP(double e);
	double getS(double e);
	double getTIS(double e);
	double getTISB(double e);
	double getSV(double e);
	double getBV(double e);
	double getTA(double e);
	double getT(double e);
	double getSigmaA(double e);
	double getSigmaB(double e);
	double getISAlpha(double e);
	double getISGamma(double e);
	inline double getISA(double e) {return 0.;}
	inline double getISB(double e) {return 0.;}
	double getDISAlphaDE(double e);
	double getDISGammaDE(double e);
	double getISBeta(double e);
	
	double getEGivenT(double T);
};
//#include "CEosIdeal.cpp"
#endif //__CEOSIDEAL_h_INCLUDE__