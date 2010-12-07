#ifndef __CCELL_h_INCLUDE__
#define __CCELL_h_INCLUDE__
/*
 *  Ccell.h
 *  v3
 *
 *  Created by Joshua Vredevoogd on 2/11/09.
 
 * Clarifications:
  (1) CCells do not use (Bjorken) tau, but log(tau) or log(sinh(tau))
  This is done for numerical smoothness. 
  The result is that while users may request time steps of 0.01 fm/c, 
    we may start out taking smaller steps and then increase.
  However, the inputs and outputs are corrected so the rest of the code is done in tau.
  (JAV 4-22-09)

  (2) For the same reasons we use log(energy density). 
  Similar i/o corrections have been implimented.
  (JAV 4-22-09)
 */
 
#include <cstdlib>
#include <cmath>
#include "hydroDef.h"
#include <cstdio>
#include "coral.h" 
#include "CEos.h"

class CCell{ 

private:
		// local variables
	double x[4];					// cell location
	CCell* neighbors[3][2];		// pointers to cell's neighbors
	double s[11];					// [(u_0,u_i); e; a_i; b]
	double dUdT[3];				// time derivative of u
	bool active;
  
		// calculation variables
		//  static double value;			// junk variable, used by a bunch of functions
	static double uup;            // u_i u_j T_ij
	static double up[3];          // u_i T_ij
	static double dS[10][3];		// mesh derivatives of s_i
	static double dULocal[4][4];	// local velocity gradients
	static double DiPixy[3][4][4];// d_i T_jk
	static double DiPixyLocal[3][3][3];// d_i T_jk
	static double DPixy[3];       // d_i T_ij
	static double Pixy[3][3];		//
	static double M[10][11];		// local matrix for calculating the time derivative
	static double dTISS[3];		// derivative of the Israel-Stewart Scaling factor
	static double dx[4];			// differential elements
  
		//working class pointers
	static CEos* eos;
	static parameterMap* pMap;
	
		//EOS variables - avoids multiple calls
	static double cS2, P, S, tIS, tISB, SV, BV, T, sigmaA, sigmaB;
	static double alphaIS, gammaIS, betaIS, aIS, bIS, dAlphaISDE, dGammaISDE;
	
	static bool mDebug, mSVTrim, mViscNS, mPureBjorken, mBjorken, mLinT, mLogT, mLogSinhT, mISVort, mISMax;
	static double mT0, mSVRatio, mBVRatio, mISAMax, mISBMax, mInitNS;
	
protected:
		// functions to fill calc variables using internal information
	double getDTxy(int);					// returns d_i T_ij
	double getDiTxy(int i,int j,int k);	// returns d_i T_jk
	void   fillM();
	double getM(int,int);					// returns entries to the augmented matrix
	void   calcDeriv();					// recalculates derivatives
	void	 fillEosVar();
	
public:
		// constructors
	CCell(parameterMap* pMap);
	CCell(double,double,double,double);	// set positions
	
		// destructor
	~CCell();
	
	void paramFill(parameterMap* pMap);
	
		// grabs:
		// cell position
	inline double getTau() {if (mLinT) return x[0]; else if (mLogT) return mT0*exp(x[0]); else if (mLogSinhT) return asinh(exp(x[0]));}
	inline double getEta() {return x[3];}
	inline double getZ() {return getTau()*sinh(x[3]);}
	inline double getTime() {return getTau()*cosh(x[3]);}
	inline double getUMesh() {return sqrt(1+pow(tanh(x[3]),2))*tanh(x[3]);}
	inline double getX() {return x[1];}
	inline double getY() {return x[2];}
	inline double getX(int i) {if (i==0) return getTau(); else if(i<4 && i>0) return x[i]; else return 0.;}
	inline double getDx(int i) {if (i<4 && i>=0) return dx[i]; else return 0.;}

		//relativistic velocities 
	inline double getU0() {return s[0];}
	inline double getUx() {return s[1];}
	inline double getUy() {return s[2];}
	inline double getUz() {return s[3];}
	inline double getU(int v) {return s[v];}
	inline double getVr() { return sqrt(s[1]*s[1]+s[2]*s[2])/s[0];}
	inline double getUzRest() { return sinh(x[3])*(1.+ (s[3]*s[3]/(s[0]+1.))) + cosh(x[3])*s[3];}
	inline double getGammaRest() {return s[0]*cosh(x[3])+s[3]*sinh(x[3]);}
	inline double getVzRest() { return getUzRest()/getGammaRest();}
	inline double getYL() { return (0.5) * log((1+getVzRest())/(1-getVzRest()));}
	inline double getGammaTrans(){ return sqrt( 1. + s[1]*s[1] + s[2]*s[2]);}
	inline double getDUDT(int v) { if (v>0) return dUdT[v-1]; else return -(dUdT[0]*s[1]+dUdT[1]*s[2]+dUdT[2]*s[3])/s[0];}
	
		//spatial components of stress energy tensor
	inline double getA(int v) {return s[4+v];}
	inline double getB() {return s[10];}
	inline double getE() {return exp(s[4]);}
	double getTxy(int,int);						// returns spatial part of tensor in fluid frame
	double getTxy(CCell*,int,int);
	double getPixy(int,int);
	double getTxyMesh(int,int);
	double getPixyMesh(int,int);
	double getTxyEos(int,int);
	double getPixyEos(int,int);
	double getTxyMeshEos(int,int);
	double getPixyMeshEos(int,int);
	double getPixyNS(int,int);
	double getPixyNSMesh(int,int);
	
		// direct eos calls
	inline double getPEos() {return eos->getP(exp(s[4]));}
	inline double getCs2Eos() {return eos->getCs2(exp(s[4]));}
	inline double getSEquilEos() {return eos->getS(exp(s[4]));}
	inline double getTISEos() {if (mSVRatio != 0.) return eos->getTIS(exp(s[4])); else return 0.;}
	inline double getTISBEos() {if (mBVRatio != 0.) return eos->getTISB(exp(s[4])); else return 0.;}
	inline double getTEos() {return eos->getT(exp(s[4]));}
	inline double getSigmaAEos() {return eos->getSigmaA(exp(s[4]));}  
	inline double getSigmaBEos() {return eos->getSigmaB(exp(s[4]));}  
	inline double getISAlphaEos(){return eos->getISAlpha(exp(s[4]));}
	inline double getISGammaEos(){return eos->getISGamma(exp(s[4]));}
	inline double getISBetaEos(){return eos->getISBeta(exp(s[4]));}
	inline double getISAEos()    {return eos->getISA(exp(s[4]));}
	inline double getISBEos()    {return eos->getISB(exp(s[4]));}
	inline double getDISAlphaDEEos(){return eos->getDISAlphaDE(exp(s[4]));}
	inline double getDISGammaDEEos(){return eos->getDISGammaDE(exp(s[4]));}
	inline double getSVEos() {if (mSVTrim) return eos->getSV(exp(s[4])) / (1. + exp( (sqrt(x[1]*x[1]+x[2]*x[2]) - 10.)/0.6)); else return eos->getSV(exp(s[4]));}
	inline double getBVEos() {return eos->getBV(exp(s[4]));}
	
		// avoids eos call
	inline double getP() {return P;}
	inline double getCs2() {return cS2;}
	double getS();
	double getSEos();
	inline double getTIS() {if (mSVRatio != 0.) return tIS; else return 0.;}
	inline double getTISB() {if (mBVRatio != 0.) return tISB; else return 0.;}
	inline double getT() {return T;}
	inline double getSigmaA() {return sigmaA;}  
	inline double getSigmaB() {return sigmaB;}  
	inline double getISAlpha(){return alphaIS;}
	inline double getISGamma(){return gammaIS;}
	inline double getISBeta(){return betaIS;}
	inline double getISA()    {return aIS;}
	inline double getISB()    {return bIS;}
	inline double getISAMax() {return mISAMax;}
	inline double getISBMax() {return mISBMax;}
	inline double getDISAlphaDE(){return dAlphaISDE;}
	inline double getDISGammaDE(){return dGammaISDE;}
	inline double getSV(){return SV;}
	inline double getBV(){return BV;}
	
		// as they are in the s vector
	inline double getDS(int i,int j) { return dS[i][j];}
	inline double getS(int i) {if (i!= 4) return s[i]; else return exp(s[4]);}
	double getSpLocDS(int,int);
	double getDiPixyLocal(int,int,int);
	
		// working classes
	inline CEos* getEos() {return eos;} 
	
		// Variable sets:
  
		// cell positions
	inline void setTau(double v) {if (mLinT) x[0]=v; else if (mLogT) x[0] = log(v/mT0); else if (mLogSinhT) x[0] = log(sinh(v));}
	inline void setEta(double v) {x[3] = v;}
	inline void setX(double v)   {x[1] = v;}
	inline void setY(double v)   {x[2] = v;}

	inline void setDt(double v) {dx[0] = v;}

		// relativistic velocities  (gamma v)
	void setU(double,double,double);
	inline void setUx(double v) {s[1] = v; s[0] = sqrt(1.+s[1]*s[1]+s[2]*s[2]+s[3]*s[3]);}
	inline void setUy(double v) {s[2] = v; s[0] = sqrt(1.+s[1]*s[1]+s[2]*s[2]+s[3]*s[3]);}
	inline void setUz(double v) {s[3] = v; s[0] = sqrt(1.+s[1]*s[1]+s[2]*s[2]+s[3]*s[3]);}
	
		// spatial components of stress energy tensor
	inline void setA(int v, double v1) {if (v>0 && v<6) s[4+v] = v1;}
	inline void setB(double v) {s[10] = v;}
	inline void setE(double v) {s[4] = log(v);}
	void setS(int i, double v);
	
		//working classes
	inline void setEos(CEos* p) {eos = p;}
	
		// set neighbors
	inline void setEtaNeighbors(CCell* p1, CCell* p2) {neighbors[2][0]=p1; neighbors[2][1]=p2;}
	inline void setXNeighbors(CCell* p1, CCell* p2) {neighbors[0][0]=p1; neighbors[0][1]=p2;}
	inline void setYNeighbors(CCell* p1, CCell* p2) {neighbors[1][0]=p1; neighbors[1][1]=p2;}
	
		// updates: up[i], uup
		// calls: calcDeriv(), FillM()
		// called by: forward()
	void update();
	void update(int);  // update an unconnected edge
	void update(CCell*);
  
	void initNS();
	
		// active?
	void deactivate();
 
		// integrates the cell forward in time
		// puts result into cell passed
	void forward(CCell*);
	void forward(CCell*, CCell*);
	void forward(CCell*, CCell*, CCell*, CCell*);
  
		// prints
	void print();			// ten local variables (no EOS dep stuff)
	void print(int);
	void printM();		// state of dt coeff matrix
	void printMm(double Mm[10][11]);

};
//#include "CCell.cpp"
#endif //__CCELL_h_INCLUDE__
