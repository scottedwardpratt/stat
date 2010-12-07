#ifndef __CHYDRO_h_INCLUDE__
#define __CHYDRO_h_INCLUDE__

#include <cstdlib>
#include <cstdio> 
#include <cmath>

#include <ctime>
#include <math.h>
#include <string.h>
#include <parametermap.h>

#include "hydroDef.h"
#include "CMesh.h"
#include "CEos.h"

class CHydro {
private:
	CMesh *onMesh, *offMesh, *tempMesh, *k1, *k2, *deadMesh;
	CEos*  mEos;

	FILE *fEX, *fEY, *fEZ, *fS, *fT, *fP;
	FILE *fEN, *fU0, *fUx, *fUy, *fUz, *fULz, *fUzz;
	FILE *fE0, *fTxx0, *fTzz0;
	FILE *fIE, *fIEN, *fIS, *fIS2, *fISN, *fIEP;
	FILE *fIPX, *fIPY, *fIPN, *fPXLoss, *fPYLoss, *fPNLoss;
	FILE *fA1x, *fA2x, *fA3x, *fA4x, *fA5x, *fBx;
	FILE *fA1y, *fA2y, *fA3y, *fA4y, *fA5y, *fBy;
	FILE *fA1z, *fA2z, *fA3z, *fA4z, *fA5z, *fBz;
	FILE *fTxx, *fTyy, *fTzz, *fTxy, *fTxz, *fTyz, *fTxxNS;
	FILE *fTrr, *fTpp, *fTpr;
	FILE *fSLoss, *fELoss, *fSLoss2;
	FILE *fUL, *fEL, *fA1L, *fA2L, *fBL;
	FILE *fTIS, *fTISB, *fDzz, *fT0;
	FILE *fEP, *fEP0, *fEX0;
	FILE *fVR0;
	FILE *fFOS, *fFOSL, *fFOSSST;
	FILE *fTQMX, *fTQMY, *fTQMR;
	FILE *fTxxL, *fTxxNSL,*fTxxML, *fTxxNSML;
	FILE *fFull, *fOscarFull, *fOscarHyper;
	FILE *fPiDNDY, *fPiDNDPt;

	static string mDataRoot, mOldFileName;
	static bool mDebug, mIoIntegrals, mIoSlices, mIoSpots, mIoTechqm, mIoOscarFull, mIoOscarHyper, mIoFull, mIoSpectra, mOldFile;
	static bool mSawtooth, mRK2, mRK4, mLinT, mLogT, mLogSinhT;
	static bool mBjorken, mPureBjorken, mHalving, mOctant;
	static int mTStep, mXSize, mYSize, mNSize, mIoSliceTStep, mIoTechqmTStep, mIoOscarTStep, mPrintStep, mPrintX, mPrintY, mPrintN;
	static double mInitNS, mInitFlow, mFoTemp, mT0, mDt, mE0;

	static double intS,intE,intP[3],sLoss,sLoss2,eLoss,intS0,intE0,intP0[3],pLoss[3],epLoss;

	void printE(CMesh* lMesh,int lT);
	void printMesh(CMesh* lMesh);
	void openFileIntegrals();
	void closeFileIntegrals();
	void printIntegrals(int t);
	void printSpots();
	void openFileSpots();
	void closeFileSpots();
	void openFile(int lT);
	void closeFile();
	void printTQM(CMesh* lMesh);
	void openOscarHyper();
	void openOscarFull();
	void printOscarHyper(CMesh* lMesh, int mT);
	void printOscarFull(CMesh* lMesh, int mT);
	void printDNs(CMesh* lMesh,int,int,int);

	void copyCellActive(CMesh*, CMesh*);

	void testFileOpen();
	void zeroPointers();
	
	inline int min(int a, int b) {if (a>b) return b; else return a;}

	time_t start, now;

	parameterMap* pMap;
	
public:
	CHydro(parameterMap*);
	int runHydro();
};

#endif __CHYDRO_h_INCLUDE__
