#ifndef __CMESH_h_INCLUDE__
#define __CMESH_h_INCLUDE__
/*  Cmesh.h
 *  isHydro3
 *  Created by Joshua Vredevoogd on 2/11/09.
 */

#include "CCell.h"
#include "hydroDef.h"
#include "CEos.h"
#include <cstdio>
#include <coral.h>

class CMesh {
  private:
	// the cells contained
	CCell* mCells[2*NSIZE+1][2*XSIZE+1][2*YSIZE+1];
	
	static parameterMap *pMap;
	//parameters from map
	static bool mOctant, mPureBjorken, mBjorken, mPrintMs, mSVTrimInit, mSVTrim, mSawtooth;
	static double mE0, mFOTemp, mDeadT, mT0, mDx, mDy, mDn;
	static int mNSize, mXSize, mYSize, mNSizeOrig, mXSizeOrig, mYSizeOrig;
	static int mPrintX, mPrintY, mPrintN;
	
	static double wnRho0, wnRAu, wnXi, wnSigma, wnA, wnB, wnK;
	static double mRn, mRt, mRa, mRx, mRy, mRSig;
	static double mWnBinRatio, mInitFlow, mInitNS;
	
	double sLoss, sLoss2, eLoss, epLoss;
	
	// to assign initial conditions to the cells
	void initialCondition(CCell*, int, int, int);  
	
	// for 2d wounded nucleon initial condition generation
	double wnE(double, double);
	double wnT(double, double);
	double wnRho(double, double, double);
	
  public:
    // con/de-structors
	CMesh(parameterMap*);
	CMesh(const char*);   //OLD_FILE
//	CMesh(string);
	CMesh(CMesh* mesh);
	~CMesh();
	
	void initNS();
	void addInitialFlow();
	
	// integrate this mesh forward 
	void forward(CMesh*, double);
	void forward(CMesh*, CMesh*, double);
	void forward(CMesh*, CMesh*, CMesh*, CMesh*, double);
	
	//puts average of the two meshs here
	void average(CMesh*,CMesh*);
	
	void smoothEdge(CCell*,CCell*,CCell*,CCell*);
	void flipEdge(int,CCell*,CCell*,CCell*);
	
	void calcLoss();
	
	void setParameterMap(parameterMap* pM)  {pMap = pM;}
	void setEos(CEos* p) {mCells[0][0][0]->setEos(p);}
	
	void printCell(int eta,int x,int y);
	void printCell(int eta, int x, int y, int opt);
	
	double getTau(){return mCells[1][1][1]->getTau();}
	void setTau(double);
	
	inline bool getActive(int eta, int x, int y) 
	  {if (mOctant) return mCells[eta+1][x+1][y+1]->getActive(); else return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getActive();}
	inline void update(int eta, int x, int y) 
	  {if (mOctant) mCells[eta+1][x+1][y+1]->update(); else mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->update();}
	inline double getS(int eta, int x, int y, int s)
	  {if (!mOctant) return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getS(s); else return mCells[eta+1][x+1][y+1]->getS(s);}
	inline double getDS(int eta, int x, int y, int m, int n) 
	  {if (mOctant) return mCells[eta+1][x+1][y+1]->getDS(m,n); else return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getDS(m,n);}
	inline double getE(int eta, int x, int y)
	  { return getS(eta,x,y,4);}
    inline double getT(int eta, int x, int y)
	  {if (!mOctant) return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getTEos(); else return mCells[eta+1][x+1][y+1]->getTEos();}
    inline double getP(int eta, int x, int y)
	  {if (!mOctant) return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getPEos(); else return mCells[eta+1][x+1][y+1]->getPEos();}
	inline double getS(int eta, int x, int y)
	  { if(!mOctant) return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getS(); else return mCells[eta+1][x+1][y+1]->getS();}
	inline double getISAlpha(int eta, int x, int y)
	  { if(!mOctant) return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getISAlphaEos(); else return mCells[eta+1][x+1][y+1]->getISAlphaEos();}
	inline double getISGamma(int eta, int x, int y)
	  { if(!mOctant) return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getISGammaEos(); else return mCells[eta+1][x+1][y+1]->getISGammaEos();}
	inline double getTIS(int eta, int x, int y)
	  { if(!mOctant) return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getTISEos(); else return mCells[eta+1][x+1][y+1]->getTISEos();}
	inline double getTISB(int eta, int x, int y)
	  { if(!mOctant) return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getTISBEos(); else return mCells[eta+1][x+1][y+1]->getTISBEos();}
	inline double getX(int eta, int x, int y,int i)
	  { if (mOctant) return mCells[eta+1][x+1][y+1]->getX(i); else return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getX(i);}
	inline double getTxy(int eta, int x, int y, int xx, int yy)
	  { if (mOctant) return mCells[eta+1][x+1][y+1]->getTxyEos(xx,yy); else return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getTxyEos(xx,yy);}
	inline double getTxyMesh(int eta, int x, int y, int xx, int yy)
	  { if (mOctant) return mCells[eta+1][x+1][y+1]->getTxyMeshEos(xx,yy); else return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getTxyMeshEos(xx,yy);}
    inline double getPixy(int eta, int x, int y, int xx, int yy)
	  { if (mOctant) return mCells[eta+1][x+1][y+1]->getPixyEos(xx,yy); else return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getPixyEos(xx,yy);}
	inline double getPixyMesh(int eta, int x, int y, int xx, int yy)
	  { if (mOctant) return mCells[eta+1][x+1][y+1]->getPixyMeshEos(xx,yy); else return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getPixyMeshEos(xx,yy);}
	double getPixyNS(int eta, int x, int y, int xx, int yy);
	double getPixyNSMesh(int eta, int x, int y, int xx, int yy);
	inline double getTxx(int eta, int x, int y) 
	  {return getTxy(eta,x,y,1,1);}
	inline double getTyy(int eta, int x, int y) 
	  {return getTxy(eta,x,y,2,2);}
	inline double getTzz(int eta, int x, int y) 
	  {return getTxy(eta,x,y,3,3);}
	inline double getSV (int eta, int x, int y) 
	  {if (mOctant) return mCells[eta+1][x+1][y+1]->getSVEos(); else return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getSVEos();}
	inline double getBV (int eta, int x, int y) 
	  {if (mOctant) return mCells[eta+1][x+1][y+1]->getBVEos(); else return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getBVEos();}	
	
	double getTZR(int eta, int x, int y);
	double getTZPhi(int eta, int x, int y);
	double getTRR(int eta, int x, int y);
	double getTRPhi(int eta, int x, int y);
	double getTPhiPhi(int eta, int x, int y);
	
	double getPiZR(int eta, int x, int y);
	double getPiZPhi(int eta, int x, int y);
	double getPiRR(int eta, int x, int y);
	double getPiRPhi(int eta, int x, int y);
	double getPiPhiPhi(int eta, int x, int y);
	
	inline double getYL( int eta, int x, int y) 
	  {if (mOctant) return mCells[eta+1][x+1][y+1]->getYL(); else return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getYL();}
	inline double getVr( int eta, int x , int y) 
	  {if (mOctant) return mCells[eta+1][x+1][y+1]->getVr(); else return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getVr();}	
	inline double getGammaTrans( int eta, int x , int y) 
	  {if (mOctant) return mCells[eta+1][x+1][y+1]->getGammaTrans(); else return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getGammaTrans();}
	inline double getDUDT(int eta, int x, int y, int v) 
	  {if (mOctant) return mCells[eta+1][x+1][y+1]->getDUDT(v); else return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getDUDT(v);}
	
	double getDULocal(int eta, int x, int y, int u, int v);
	
	inline int getNSize() {return mNSize;}
	inline int getXSize() {return mXSize;}
	inline int getYSize() {return mYSize;}
	inline int getNSizeOrig() {return mNSizeOrig;}
	inline int getXSizeOrig() {return mXSizeOrig;}
	inline int getYSizeOrig() {return mYSizeOrig;}
	inline void setNSize(int i) {mNSize = i;}
	inline void setXSize(int i) {mXSize = i;}
	inline void setYSize(int i) {mYSize = i;}

	inline void setActive(int eta, int x, int y, bool v) 
	  {if (mOctant) mCells[eta+1][x+1][y+1]->setActive(v); else mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->setActive(v);}

	double integralE(int);
	double integralE();
	double integralP(int,int);
	double integralP(int);
	double integralS(int eta);
	double integralS();
	double integralEp(int eta);
	double integralEp();
	double integralEx(int eta);
	double integralEx();
	double averageVr(int eta);
	double averageVr();
	
	double getdNd3p(double px, double py, double pz, double m);
	
	void getFOS(double x[XSIZE+YSIZE][2], double y[XSIZE+YSIZE][2], int &foSize);
	void getFOS(double x[XSIZE+YSIZE][2], double y[XSIZE+YSIZE][2], double y[XSIZE+YSIZE][2], double y[XSIZE+YSIZE][5], int &foSize);
	double getFOP();
	double getFOSSST();
	
	inline double getSLoss() {return sLoss;}
	inline double getSLoss2(){return sLoss2;}
	double getELoss(CMesh*);
	double getELoss(CMesh*,int);
	double getPLoss(CMesh*,int);
	double getPLoss(CMesh*,int,int);
	
	inline void deactivate(int eta, int x, int y) 
	  {if (mOctant) mCells[eta+1][x+1][y+1]->deactivate(); else mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->deactivate();}
	
	// remove the cells 
	void deaden();
	void deaden(CMesh*);
	bool detectCrash();
	
	void checkAzimuthalSymmetry();

};
//#include "CMesh.cpp"
#endif // __CMESH_h_INCLUDE__
