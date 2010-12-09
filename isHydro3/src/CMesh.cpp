/*  Cmesh.cpp
 *  isHydro3
 *  Created by Joshua Vredevoogd on 2/11/09.
 */

#include "CMesh.h"
#include <iostream>
#include <fstream>
#include "CEos.h"
#include "coral.h"

CMesh::CMesh(parameterMap* pM){
 pMap = pM;
 CCell* dummyCell = new CCell(pMap);
 delete dummyCell;

 mE0 = parameter::getD(*pMap,"HYDRO_E0",1.0);

 mOctant = parameter::getB(*pMap,"HYDRO_OCTANT",true);
 mPureBjorken = parameter::getB(*pMap,"HYDRO_PURE_BJORKEN",false);
 mBjorken = parameter::getB(*pMap,"HYDRO_BJORKEN",true);
 mPrintMs = parameter::getB(*pMap,"HYDRO_PRINTMS",false);
 mSVTrimInit = parameter::getB(*pMap,"HYDRO_SVTRIMINIT",false);
 mSVTrim = parameter::getB(*pMap,"HYDRO_SVTRIM",false);
 mFOTemp = parameter::getD(*pMap,"HYDRO_FOTEMP",0.13);
 mDeadT = parameter::getD(*pMap,"HYDRO_DEADT",0.1);

 mNSize = parameter::getI(*pMap,"HYDRO_NSIZE",20);
 mXSize = parameter::getI(*pMap,"HYDRO_XSIZE",60);
 mYSize = parameter::getI(*pMap,"HYDRO_YSIZE",60);
 mNSizeOrig = mNSize;
 mXSizeOrig = mXSize;
 mYSizeOrig = mYSize;
 
 mPrintN = parameter::getI(*pMap,"HYDRO_PRINTN",0);
 mPrintX = parameter::getI(*pMap,"HYDRO_PRINTX",0);
 mPrintY = parameter::getI(*pMap,"HYDRO_PRINTY",0);

 mT0 = parameter::getD(*pMap,"HYDRO_T0",1.0);
 mDn = parameter::getD(*pMap,"HYDRO_DN",0.1);
 mDx = parameter::getD(*pMap,"HYDRO_DX",0.1);
 mDy = parameter::getD(*pMap,"HYDRO_DY",0.1);

 wnRho0 = parameter::getD(*pMap,"GLAUBER_Rho0",0.16);
 wnRAu  = parameter::getD(*pMap,"GLAUBER_RAu",6.3);
 wnXi   = parameter::getD(*pMap,"GLAUBER_Xi",0.3);
 wnSigma= parameter::getD(*pMap,"GLAUBER_Sigma",0.4);
 wnA    = parameter::getD(*pMap,"GLAUBER_A",197.);
 wnB    = parameter::getD(*pMap,"GLAUBER_B",0.);
 wnKTau = parameter::getD(*pMap,"GLAUBER_K_TAU",4.3);

 mWnBinRatio = parameter::getD(*pMap,"GLAUBER_WNBIN_RATIO",1.0);
 mInitFlow = parameter::getD(*pMap,"HYDRO_INIT_FLOW",0.0);
 mInitNS = parameter::getD(*pMap,"HYDRO_INIT_NS",0.0);

 if (!mOctant) {
	// make cells
	for (int i=-mNSize;i<=mNSize;i++)
	  for (int j=-mXSize;j<=mXSize;j++)
		for (int k=-mXSize;k<=mYSize;k++) {
		  mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig] = new CCell(mT0, 
														mDn*(double)(i), 
														mDx*(double)(j), 
														mDy*(double)(k));
		  initialCondition( mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig],i,j,k);
		}

	// connect your cells and reduce boundary problem...
	for (int i=-mNSizeOrig;i<=mNSizeOrig;i++)
	  for (int j=-mXSizeOrig;j<=mXSizeOrig;j++)
		for (int k=-mYSizeOrig;k<=mYSizeOrig;k++) {
		  if(i==-mNSizeOrig)	
		    mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]->setEtaNeighbors( mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig],
																mCells[i+1+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]);
		  else if (i==mNSizeOrig) 
		    mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]->setEtaNeighbors( mCells[i-1+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig], 
																mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]);
		  else				 
		    mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]->setEtaNeighbors( mCells[i-1+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig], 
																mCells[i+1+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]);
		  
		  if(j==-mXSizeOrig)		 
		    mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]->setXNeighbors( mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig], 
															  mCells[i+mNSizeOrig][j+1+mXSizeOrig][k+mYSizeOrig]);
		  else if (j==mXSizeOrig)
		    mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]->setXNeighbors( mCells[i+mNSizeOrig][j-1+mXSizeOrig][k+mYSizeOrig], 
															  mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]);
		  else				 
			mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]->setXNeighbors( mCells[i+mNSizeOrig][j-1+mXSizeOrig][k+mYSizeOrig], 
															  mCells[i+mNSizeOrig][j+1+mXSizeOrig][k+mYSizeOrig]);
		  
		  if(k==-mYSizeOrig)		 
		    mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]->setYNeighbors( mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig], 
															  mCells[i+mNSizeOrig][j+mXSizeOrig][k+1+mYSizeOrig]);
		  else if (k==mYSizeOrig) 
		    mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]->setYNeighbors( mCells[i+mNSizeOrig][j+mXSizeOrig][k-1+mYSizeOrig], 
															  mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]);
		  else				 
		    mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]->setYNeighbors( mCells[i+mNSizeOrig][j+mXSizeOrig][k-1+mYSizeOrig], 
															  mCells[i+mNSizeOrig][j+mXSizeOrig][k+1+mYSizeOrig]);
		}
  } 
 else { //mOctant
	// make cells
	for (int i=-1;i<=mNSize;i++)
	  for (int j=-1;j<=mXSize;j++)
		for (int k=-1;k<=mYSize;k++) {
		  mCells[i+1][j+1][k+1] = new CCell(mT0, 
											mDn*(double)(i), 
											mDx*(double)(j),
											mDy*(double)(k));
		  initialCondition( mCells[i+1][j+1][k+1],i,j,k);
		}

	// connect your cells and (antiquated..) reduce boundary problem...
	for (int i=0;i<=mNSize;i++)
	  for (int j=0;j<=mXSize;j++)
		for (int k=0;k<=mYSize;k++) {
			if (i<mNSize)
				mCells[i+1][j+1][k+1]->setEtaNeighbors( mCells[i][j+1][k+1], 
													   mCells[i+2][j+1][k+1]);
			else 
				mCells[i+1][j+1][k+1]->setEtaNeighbors( mCells[i][j+1][k+1], 
													   mCells[i+1][j+1][k+1]);
			if (j<mXSize)
				mCells[i+1][j+1][k+1]->setXNeighbors( mCells[i+1][j][k+1], 
													   mCells[i+1][j+2][k+1]);
			else 
				mCells[i+1][j+1][k+1]->setXNeighbors( mCells[i+1][j][k+1], 
													   mCells[i+1][j+1][k+1]);
			if (k<mYSize)
				mCells[i+1][j+1][k+1]->setYNeighbors( mCells[i+1][j+1][k], 
													   mCells[i+1][j+1][k+2]);
			else 
				mCells[i+1][j+1][k+1]->setYNeighbors( mCells[i+1][j+1][k], 
													   mCells[i+1][j+1][k+1]);
		}
 }
	
	deaden();
}

CMesh::CMesh(const char *fName) {
//CMesh::CMesh(const string fName) {
	double mX[4];
	double mS[11];
	std::ifstream mFile(fName);

  if (mOctant) {
	// make cells
	mFile >> mX[0] >> mNSize >> mXSize >> mYSize;
	if (!mPureBjorken)
	for (int i=-1;i<=mNSizeOrig;i++)
	  for (int j=-1;j<=mXSizeOrig;j++)
		for (int k=-1;k<=mYSizeOrig;k++) {
		  for (int l=1;l<4;l++) 
			mFile >> mX[l];
		  for (int l=0;l<11;l++)
			mFile >> mS[l];
		  mCells[i+1][j+1][k+1] = new CCell(mX[0],mX[1],mX[2],mX[3]);
		  for (int l=1;l<11;l++) mCells[i+1][j+1][k+1]->setS(l,mS[l]);
		}
	else 
	for (int j=-1;j<=mXSizeOrig;j++)
		for (int k=-1;k<=mYSizeOrig;k++) {
		  for (int l=1;l<4;l++) 
			mFile >> mX[l];
		  for (int l=0;l<11;l++)
			mFile >> mS[l];
		  for (int i=-1;i<=mNSizeOrig;i++) {
		    mCells[i+1][j+1][k+1] = new CCell(mX[0],mX[1],mX[2],mX[3]);
		    for (int l=1;l<11;l++) mCells[i+1][j+1][k+1]->setS(l,mS[l]);
		  }
		}

	// connect your cells...
	if (!mPureBjorken)
	for (int i=0;i<mNSizeOrig;i++)
	  for (int j=0;j<mXSizeOrig;j++)
		for (int k=0;k<mYSizeOrig;k++) {
		    mCells[i+1][j+1][k+1]->setEtaNeighbors( mCells[i][j+1][k+1], 
																mCells[i+2][j+1][k+1]);
			mCells[i+1][j+1][k+1]->setXNeighbors( mCells[i+1][j][k+1], 
															  mCells[i+1][j+2][k+1]);
		    mCells[i+1][j+1][k+1]->setYNeighbors( mCells[i+1][j+1][k], 
															  mCells[i+1][j+1][k+2]);
		}
	else 
	for (int i=0;i<mNSizeOrig;i++)
	  for (int j=0;j<mXSizeOrig;j++)
		for (int k=0;k<mYSizeOrig;k++) {
			mCells[i+1][j+1][k+1]->setXNeighbors( mCells[i+1][j][k+1], 
															  mCells[i+1][j+2][k+1]);
		    mCells[i+1][j+1][k+1]->setYNeighbors( mCells[i+1][j+1][k], 
															  mCells[i+1][j+1][k+2]);
		}
  }
}

CMesh::CMesh(CMesh* mesh) {
  for (int i=-1;i<=mNSize;i++)
	  for (int j=-1;j<=mXSize;j++)
		for (int k=-1;k<=mYSize;k++)
			mCells[i+1][j+1][k+1] = new CCell(mesh->mCells[i+1][j+1][k+1]);
		  
		// connect your cells and (antiquated..) reduce boundary problem...
	for (int i=0;i<=mNSize;i++)
		for (int j=0;j<=mXSize;j++)
			for (int k=0;k<=mYSize;k++) {
				if (i<mNSize)
					mCells[i+1][j+1][k+1]->setEtaNeighbors( mCells[i][j+1][k+1], 
														   mCells[i+2][j+1][k+1]);
				else 
					mCells[i+1][j+1][k+1]->setEtaNeighbors( mCells[i][j+1][k+1], 
														   mCells[i+1][j+1][k+1]);
				if (j<mXSize)
					mCells[i+1][j+1][k+1]->setXNeighbors( mCells[i+1][j][k+1], 
														 mCells[i+1][j+2][k+1]);
				else 
					mCells[i+1][j+1][k+1]->setXNeighbors( mCells[i+1][j][k+1], 
														 mCells[i+1][j+1][k+1]);
				if (k<mYSize)
					mCells[i+1][j+1][k+1]->setYNeighbors( mCells[i+1][j+1][k], 
														 mCells[i+1][j+1][k+2]);
				else 
					mCells[i+1][j+1][k+1]->setYNeighbors( mCells[i+1][j+1][k], 
														 mCells[i+1][j+1][k+1]);
			}
}

// destructor
CMesh::~CMesh() {
	delete [] mCells;
}

// assigns the cells initial conditions
void CMesh::initialCondition(CCell* mCell, int i, int j, int k) {

	mCell->setU( 0., 0., 0.);

		//  if ((mDn*i) < 3.)
		//    mCell->setE( wnK * wnE(mDx*(double)j, mDy*(double)k));
		//  else 
		//    mCell->setE( wnK * wnE(mDx*(double)j, mDy*(double)k) * exp(-0.5*pow( (mDn*i - 3.)/ 0.4, 2.)));

  mCell->setE( wnKTau/getX(i,j,k,0)  * wnE(mDx*(double)j, mDy*(double)k) * exp( -0.5*pow( (mDn*i)/1.6, 2.)));

	for (int l=5;l<11;l++) mCell->setS(l,0.);
}

void CMesh::addInitialFlow() {

	if (!mOctant) 
		for (int i=-mNSize+1;i<mNSize;i++)
			for (int j=-mXSize+1;j<mXSize;j++)
				for (int k=-mXSize+1;k<mYSize;k++) {
		  mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]->update();
		  mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]->setS(1, - 0.5*getDS(i,j,k,0,3)/( getE(i,j,k) + getP(i,j,k)) * getTau());
		  mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]->setS(2, - 0.5*getDS(i,j,k,1,3)/( getE(i,j,k) + getP(i,j,k)) * getTau());
		}
	else 
		if (!mPureBjorken)
			for (int i=0;i<=mNSize;i++)
				for (int j=0;j<=mXSize;j++)
					for (int k=0;k<=mYSize;k++) {
						if (!getActive(i,j,k)) break;
						update(i,j,k);
						mCells[i+1][j+1][k+1]->setS(1, - mInitFlow*0.5*getDS(i,j,k,3,0)/( getE(i,j,k) + getP(i,j,k)) * getTau());
						mCells[i+1][j+1][k+1]->setS(2, - mInitFlow*0.5*getDS(i,j,k,3,1)/( getE(i,j,k) + getP(i,j,k)) * getTau());
					}
		else 
			for (int j=0;j<=mXSize;j++)
				for (int k=0;k<=mYSize;k++) {
					if (!getActive(0,j,k)) break;
					update(0,j,k);
					mCells[1][j+1][k+1]->setS(1, - mInitFlow*0.5*getDS(0,j,k,3,0)/( getE(0,j,k) + getP(0,j,k)) * getTau());
					mCells[1][j+1][k+1]->setS(2, - mInitFlow*0.5*getDS(0,j,k,3,1)/( getE(0,j,k) + getP(0,j,k)) * getTau());
				}
	flipEdge();
}

double CMesh::wnRho(double x, double y, double z) {
	return wnRho0 / ( 1. + exp( (sqrt(x*x+y*y+z*z) - wnRAu)/wnXi));
}

double CMesh::wnT(double x, double y) {
	double wnDz = wnXi/2.;
	double value = 0.; 
	double i = -wnRAu*3.;

	for (;i < wnRAu*3.; i+=3*wnDz) 
	  value += (3.*wnDz/8.) * (wnRho(x,y,i) + 3.*wnRho(x,y,i+wnDz)
							+ 3.*wnRho(x,y,i+2*wnDz) + wnRho(x,y,i+3*wnDz));

	return value;
}

double CMesh::wnE(double x, double y) {
	double value1 = wnT(x+wnB/2.,y);
	double value2 = wnT(x-wnB/2.,y);
	return mWnBinRatio*(  value1 * ( 1. - pow(( 1. - (wnSigma/wnA)*value2),wnA))
						+ value2 * ( 1. - pow(( 1. - (wnSigma/wnA)*value1),wnA)))
		   + (1. - mWnBinRatio) * wnSigma * value1 * value2;
}

// calculates how this cell should change
// and puts the answer into the provided mesh, same cell
void CMesh::forward(CMesh* mMesh, double mDt) {
  // setting dt anywhere changes it everywhere

  mCells[1][1][1]->setDt(mDt);
  sLoss = 0.; epLoss = 0.; eLoss = 0.;

  mMesh->setNSize(mNSize);
  mMesh->setXSize(mXSize);
  mMesh->setYSize(mYSize);

 if (mPrintMs) {
  std::cout << std::endl << "pre-half..." << std::endl;
  mCells[mPrintN+1][mPrintX+1][mPrintY+1]->update();
  mCells[mPrintN+1][mPrintX+1][mPrintY+1]->print();
  mCells[mPrintN+1][mPrintX+1][mPrintY+1]->printM();
 }
 
 if (mOctant) {
  if (!mPureBjorken) {
	  for (int i=0;i<=mNSize;i++)
		  for (int j=0;j<=mXSize;j++)
			  for (int k=0;k<=mYSize;k++) 
				  if (getActive(i,j,k))
					  mCells[i+1][j+1][k+1]->forward(mMesh->mCells[i+1][j+1][k+1]);
				  else break;
	
	  mMesh->flipEdge();
  } else { //mPureBjorken
	  for (int j=0;j<=mXSize;j++)
		  for (int k=0;k<=mYSize;k++) 
			  if (getActive(0,j,k))
				  mCells[1][j+1][k+1]->forward(mMesh->mCells[1][j+1][k+1]);
			  else break;
	  
	  mMesh->flipEdge();
  }
 }  // !mOctant
 else { // !mOctant
  for (int i=-mNSize;i<=mNSize;i++)
    for (int j=-mXSize;j<=mXSize;j++)
	  for (int k=-mYSize;k<=mYSize;k++) 
	    mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]->forward( mMesh->mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]);
 }
 
 if (mPrintMs) {
   std::cout << "post-half" << std::endl;
   mMesh->mCells[mPrintN+1][mPrintX+1][mPrintY+1]->update();
   mMesh->mCells[mPrintN+1][mPrintX+1][mPrintY+1]->print();
   mMesh->mCells[mPrintN+1][mPrintX+1][mPrintY+1]->printM();
 }
}

void CMesh::forward(CMesh* onMesh, CMesh* offMesh, double mDt) {
		// setting dt anywhere changes it everywhere
	mCells[0][0][0]->setDt(mDt);
	eLoss = 0.;  sLoss = 0.;  epLoss = 0.;
	
	onMesh->setNSize(mNSize);
	onMesh->setXSize(mXSize);
	onMesh->setYSize(mYSize);
    
	if (mPrintMs) {
		std::cout << "pre-whole" << std::endl;
		onMesh->mCells[mPrintN+1][mPrintX+1][mPrintY+1]->update();
		onMesh->mCells[mPrintN+1][mPrintX+1][mPrintY+1]->print();
		if (mSawtooth) offMesh->mCells[mPrintN+1][mPrintX+1][mPrintY+1]->printM();
	}
  
	if (mOctant) 
		if (!mPureBjorken) 
			for (int i=0;i<=mNSize;i++)
				for (int j=0;j<=mXSize;j++)
					for (int k=0;k<=mYSize;k++) 
						if (getActive(i,j,k))
							mCells[i+1][j+1][k+1]->forward( onMesh->mCells[i+1][j+1][k+1],
														   offMesh->mCells[i+1][j+1][k+1]); 
						else break;
		else 
			for (int j=0;j<mXSize;j++)
				for (int k=0;k<mYSize;k++) 
					if (getActive(0,j,k))
						mCells[1][j+1][k+1]->forward(onMesh->mCells[1][j+1][k+1],
													 offMesh->mCells[1][j+1][k+1]);
					else break;
	else 
		for (int i=-mNSize;i<=mNSize;i++)
			for (int j=-mXSize;j<=mXSize;j++)
				for (int k=-mYSize;k<=mYSize;k++) 
					mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]->forward( onMesh->mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig],
																			  offMesh->mCells[i+mNSizeOrig][j+mXSizeOrig][k+mYSizeOrig]);

	
	onMesh->flipEdge();
	
	if (mPrintMs) {
		std::cout << "post-whole(t=" << onMesh->getTau() << ')' << std::endl;
		onMesh->mCells[mPrintN+1][mPrintX+1][mPrintY+1]->update();
		onMesh->mCells[mPrintN+1][mPrintX+1][mPrintY+1]->print();
		onMesh->mCells[mPrintN+1][mPrintX+1][mPrintY+1]->printM();
	}
}

void CMesh::forward(CMesh* k0, CMesh* k1, CMesh* k2, CMesh* k3, double mDt) {
  mCells[1][1][1]->setDt(mDt);
  if (mOctant && !mPureBjorken) {  
	
  for (int i=0; i<mNSize; i++)
    for (int jj=0; jj<mXSize; jj++) 
      for (int ll=0; ll<mYSize; ll++) 
		mCells[i+1][jj+1][ll+1]->forward(k0->mCells[i+1][jj+1][ll+1],
									     k1->mCells[i+1][jj+1][ll+1],
									     k2->mCells[i+1][jj+1][ll+1],
									     k3->mCells[i+1][jj+1][ll+1]);

//  mCells[mNSize][mXSize][mYSize]->print(0);

//fix the boundary cells
// outer eta boundaries
  for (int i=0;i<mXSize;i++)
    for (int j=0;j<mYSize;j++) {
		smoothEdge(mCells[mNSize-2][i+1][j+1],
				   mCells[mNSize-1][i+1][j+1],
				   mCells[mNSize][i+1][j+1],
				   mCells[mNSize+1][i+1][j+1]);
		mCells[mNSize+1][i+1][j+1]->setTau(mCells[mNSize][i+1][j+1]->getTau());
	}	
// mCells[mNSize][mXSize][mYSize]->print(0);

//outer radial boundaries
  for (int i=0;i<=mNSize;i++) {
    for (int j=0;j<mYSize;j++) {
		smoothEdge(mCells[i+1][mXSize-2][j+1],
				   mCells[i+1][mXSize-1][j+1],
				   mCells[i+1][mXSize][j+1],
				   mCells[i+1][mXSize+1][j+1]);
		mCells[i+1][mXSize+1][j+1]->setTau(mCells[i+1][mXSize][j+1]->getTau());
	} for (int j=0;j<mXSize;j++) {			   
		smoothEdge(mCells[i+1][j+1][mYSize-2],
		   	       mCells[i+1][j+1][mYSize-1],
				   mCells[i+1][j+1][mYSize],
				   mCells[i+1][j+1][mYSize+1]);
		mCells[i+1][j+1][mYSize+1]->setTau(mCells[i+1][j+1][mYSize]->getTau());
	}
	smoothEdge(mCells[i+1][mXSize+1][mYSize-2],
			   mCells[i+1][mXSize+1][mYSize-1],
			   mCells[i+1][mXSize+1][mYSize],
			   mCells[i+1][mXSize+1][mYSize+1]);
	mCells[i+1][mXSize+1][mYSize+1]->setTau(mCells[i+1][mXSize+1][mYSize]->getTau());
  }
// mCells[mNSize][mXSize][mYSize-1]->print(0);

//inner eta bound
  for (int i=0;i<=mXSize;i++)
    for (int j=0;j<=mYSize;j++) {
	// inner
		mCells[0][i+1][j+1]->setU(mCells[2][i+1][j+1]->getUx(), 
								  mCells[2][i+1][j+1]->getUy(), 
								 -mCells[2][i+1][j+1]->getUz());
		for (int ll=4;ll<11;ll++)
		  mCells[0][i+1][j+1]->setS(ll, mCells[2][i+1][j+1]->getS(ll));
		mCells[0][i+1][j+1]->setTau(mCells[2][i+1][j+1]->getTau());
	}
// mCells[mNSize][mXSize][mYSize-2]->print(0);
 
//inner rad bound
  for (int i=-1;i<=mNSize;i++) {
    for (int j=0;j<=mYSize;j++) {
		mCells[i+1][0][j+1]->setU(-mCells[i+1][2][j+1]->getUx(), 
								   mCells[i+1][2][j+1]->getUy(), 
								   mCells[i+1][2][j+1]->getUz());
		for (int kk=4;kk<11;kk++)
		  mCells[i+1][0][j+1]->setS(kk,mCells[i+1][2][j+1]->getS(kk));
		mCells[i+1][0][j+1]->setTau(mCells[i+1][2][j+1]->getTau());
	} for (int j=0;j<=mXSize;j++) {	
		mCells[i+1][j+1][0]->setU(mCells[i+1][j+1][2]->getUx(), 
								 -mCells[i+1][j+1][2]->getUy(), 
								  mCells[i+1][j+1][2]->getUz());									  
		for(int kk=4;kk<11;kk++)
		  mCells[i+1][j+1][0]->setS(kk,mCells[i+1][j+1][2]->getS(kk));
		mCells[i+1][j+1][0]->setTau(mCells[i+1][j+1][2]->getTau());
	}
	mCells[i+1][0][0]->setU(-mCells[i+1][2][2]->getUx(), 
							-mCells[i+1][2][2]->getUy(), 
							 mCells[i+1][2][2]->getUz());
	for (int kk=4;kk<11;kk++)
	  mCells[i+1][0][0]->setS(kk,mCells[i+1][2][2]->getS(kk));
	mCells[i+1][0][0]->setTau(mCells[i+1][2][2]->getTau());
  }
// mCells[mNSize][mXSize][mYSize]->print(0);
 } 
  else if (mOctant && mPureBjorken) {  
	for (int jj=0; jj<mXSize; jj++) 
      for (int ll=0; ll<mYSize; ll++) 
		mCells[1][jj+1][ll+1]->forward(k0->mCells[1][jj+1][ll+1],
									   k1->mCells[1][jj+1][ll+1],
									   k2->mCells[1][jj+1][ll+1],
									   k3->mCells[1][jj+1][ll+1]);
 
      // outer boundaries
	for (int j=0;j<mYSize;j++) {
		smoothEdge(mCells[1][mXSize-2][j+1],
				   mCells[1][mXSize-1][j+1],
				   mCells[1][mXSize][j+1],
				   mCells[1][mXSize+1][j+1]);
		mCells[1][mXSize+1][j+1]->setTau(mCells[1][mXSize][j+1]->getTau());
		
		if (mBjorken) {
		  sLoss += mCells[1][mXSize+1][j+1]->getS() * getTau()
				* mCells[1][mXSize+1][j+1]->getS(1)*mDn*mDy;
		  eLoss += mCells[1][mXSize+1][j+1]->getE() * getTau()
				* mCells[1][mXSize+1][j+1]->getS(1)*mDn*mDy;
		} 
		else {
		  sLoss += mCells[1][mXSize+1][j+1]->getS() 
				* mCells[1][mXSize+1][j+1]->getS(1)*mDn*mDy;
		  eLoss += mCells[1][mXSize+1][j+1]->getE() 
				* mCells[1][mXSize+1][j+1]->getS(1)*mDn*mDy;
		}
	} 
	for (int j=0;j<mXSize;j++) {		
		smoothEdge(mCells[1][j+1][mYSize-2],
		   	       mCells[1][j+1][mYSize-1],
				   mCells[1][j+1][mYSize],
				   mCells[1][j+1][mYSize+1]);
		mCells[1][j+1][mYSize+1]->setTau(mCells[1][j+1][mYSize]->getTau());
		
		if (mBjorken) {
		  sLoss += mCells[1][j+1][mYSize+1]->getS() * getTau()
				* mCells[1][j+1][mYSize+1]->getS(2)*mDx*mDn;
		  eLoss += mCells[1][j+1][mYSize+1]->getE() * getTau()
				* mCells[1][j+1][mYSize+1]->getS(2)*mDx*mDn;
		} else {
		  sLoss += mCells[1][j+1][mYSize+1]->getS() * mCells[1][j+1][mYSize+1]->getS(2)*mDx*mDn;
		  eLoss += mCells[1][j+1][mYSize+1]->getE() * mCells[1][j+1][mYSize+1]->getS(2)*mDx*mDn;
		}
				
	}
 
 // inner boundaries
    for (int j=0;j<=mYSize;j++) {
	  flipEdge(1,mCells[1][2][j+1],
			     mCells[1][1][j+1],
			     mCells[1][0][j+1]);
		mCells[1][0][j+1]->setTau(mCells[1][2][j+1]->getTau());
	}
	for (int j=0;j<=mXSize;j++) {	
		flipEdge(2,mCells[1][j+1][2],
				   mCells[1][j+1][1],
				   mCells[1][j+1][0]);
		mCells[1][j+1][0]->setTau(mCells[1][j+1][2]->getTau());
	}
 
 } else {
     //noop!
   printf(" full sphere option not implimented for RK4....\n");
 }
}

void CMesh::average(CMesh *m1, CMesh *m2) {
  // protection for averaging with yourself and another mesh (since getTau takes the central cell time, which will be changed)
  double mTau = getTau();
  double mTau1 = m1->getTau();
  double mTau2 = m2->getTau();
  
  for (int i=-1;i<=mNSize;i++)
	for (int j=-1;j<=mXSize;j++)
      for (int k=-1;k<=mYSize;k++)
		for (int s=1;s<11;s++)
		  mCells[i+1][j+1][k+1]->setS(s, m1->getS(i,j,k,s) + ((mTau - mTau1)/(mTau2 - mTau1))
															 *(m2->getS(i,j,k,s) - m1->getS(i,j,k,s)));	
}

void CMesh::setTau(double t) {
	for (int i=-1;i<mNSize+1;i++)
      for (int j=-1;j<mXSize+1;j++)
	    for (int kk=-1;kk<mYSize+1;kk++) 
		  mCells[i+1][j+1][kk+1]->setTau(t);
}

double CMesh::getDULocal(int eta, int x, int y, int u, int v){
  if (v==1) 
	return (getS(eta,x+1,y,u) - getS(eta,x-1,y,u))/(getX(eta,x+1,y,v) - getX(eta,x-1,y,v));
  if (v==2) 
	return (getS(eta,x,y+1,u) - getS(eta,x,y-1,u))/(getX(eta,x,y+1,v) - getX(eta,x,y-1,v));
  if (v==3) 
	return getS(eta,x,y,0)/getTau() 
			+ (getS(eta+1,x,y,u) - getS(eta-1,x,y,u))/(getTau()*(getX(eta+1,x,y,v) - getX(eta-1,x,y,v)));
}

double CMesh::integralE(int eta) {
	double mTemp=0.;
 if (!mOctant) {
	for (int i=-mXSize;i<=mXSize;i++)
		for (int j=-mYSize;j<=mYSize;j++) 
		  mTemp += mDx*mDy*getE(eta,i,j);
 }
 else { // mOctant
	mTemp = (mDx*mDy/9.) *(getTxyMesh(eta,0,0,0,0) 
					   - getTxyMesh(eta,0,mYSize,0,0)
					   - getTxyMesh(eta,mXSize,0,0,0)
					   - 3.* getTxyMesh(eta,mXSize,mYSize,0,0));
	
	for (int i=1;i<mXSize;i+=2)   
	  mTemp += (mDx*mDy/9.) * ( 4.*getTxyMesh(eta,i,0,0,0) 
							+ 2.*getTxyMesh(eta,i+1,0,0,0)
						    + 4.*getTxyMesh(eta,i,mYSize,0,0) 
							+ 2.*getTxyMesh(eta,i+1,mYSize,0,0));
							
	for (int i=1;i<mYSize;i+=2)   
	  mTemp += (mDx*mDy/9.) * ( 4.*getTxyMesh(eta,0,i,0,0) 
							+ 2.*getTxyMesh(eta,0,i+1,0,0)
						    + 4.*getTxyMesh(eta,mXSize,i,0,0) 
							+ 2.*getTxyMesh(eta,mXSize,i+1,0,0));
	  
	for (int i=1;i<mXSize;i++) 
	  for (int j=1;j<mYSize;j++) 
	    if ( (i+j)%2 == 1) 
		  mTemp += (8.*mDx*mDy/9.) * getTxyMesh(eta,i,j,0,0);
		else if (i%2 == 1)  
		  mTemp +=(16.*mDx*mDy/9.) * getTxyMesh(eta,i,j,0,0);
		else  
		  mTemp += (4.*mDx*mDy/9.) * getTxyMesh(eta,i,j,0,0);
 }
 return mTemp;
}

double CMesh::integralE(){
  if (mPureBjorken) 
    if (mBjorken) return mDn*getTau()*integralE(0);
	else return mDn*integralE(0);
  
 double mTemp=0.;
 if (!mOctant) {
  if (mBjorken) 
    for (int i=-mNSize;i<=mNSize;i++) 
      mTemp += mDn * getTau() * cosh(mCells[i+mNSizeOrig][0][0]->getEta()) * integralE(i);
  else 
    for (int i=-mNSize;i<=mNSize;i++) 
      mTemp += mDn * integralE(i);
 } 
 else //mOctant
  if (mBjorken) {
    mTemp += 0.5* mDn * getTau() * (integralE(0) + integralE(mNSize));
    for (int i=1;i<mNSize;i++) 
	  mTemp += mDn * getTau() * integralE(i);
  }
  else  {
    mTemp += 0.5 * mDn * (integralE(0) + integralE(mNSize));
    for (int i=1;i<mNSize;i++) 
	  mTemp += mDn * integralE(i);
  }
 return mTemp;
}

double CMesh::integralP(int mI, int eta) {
	double mTemp=0.;
 if (!mOctant) {
	for (int i=-mXSize;i<=mXSize;i++)
		for (int j=-mYSize;j<=mYSize;j++) 
		  mTemp += mDx*mDy*getE(eta,i,j);
 }
 else { // mOctant
	mTemp = (mDx*mDy/9.) *(getTxyMesh(eta,0,0,0,mI) 
					   - getTxyMesh(eta,0,mYSize,0,mI)
					   - getTxyMesh(eta,mXSize,0,0,mI)
					   - 3.* getTxyMesh(eta,mXSize,mYSize,0,mI));
	
	for (int i=1;i<mXSize;i+=2)   
	  mTemp += (mDx*mDy/9.) * ( 4.*getTxyMesh(eta,i,0,0,mI) 
							+ 2.*getTxyMesh(eta,i+1,0,0,mI)
						    + 4.*getTxyMesh(eta,i,mYSize,0,mI) 
							+ 2.*getTxyMesh(eta,i+1,mYSize,0,mI));
							
	for (int i=1;i<mYSize;i+=2)   
	  mTemp += (mDx*mDy/9.) * ( 4.*getTxyMesh(eta,0,i,0,mI) 
							+ 2.*getTxyMesh(eta,0,i+1,0,mI)
						    + 4.*getTxyMesh(eta,mXSize,i,0,mI) 
							+ 2.*getTxyMesh(eta,mXSize,i+1,0,mI));
	  
	for (int i=1;i<mXSize;i++) 
	  for (int j=1;j<mYSize;j++) 
	    if ( (i+j)%2 == 1) 
		  mTemp += (8.*mDx*mDy/9.) * getTxyMesh(eta,i,j,0,mI);
		else if (i%2 == 1)  
		  mTemp +=(16.*mDx*mDy/9.) * getTxyMesh(eta,i,j,0,mI);
		else  
		  mTemp += (4.*mDx*mDy/9.) * getTxyMesh(eta,i,j,0,mI);
 }
 return mTemp;
}

double CMesh::integralP(int mI){
  if (mPureBjorken) 
    if (mBjorken) return mDn*getTau()*integralP(mI,0);
	else return mDn*integralP(mI,0);
  
 double mTemp=0.;
 if (!mOctant) {
  if (mBjorken) 
    for (int i=-mNSize;i<=mNSize;i++) 
      mTemp += mDn * getTau() * cosh(mCells[i+mNSizeOrig][0][0]->getEta()) * integralP(mI,i);
  else 
    for (int i=-mNSize;i<=mNSize;i++) 
      mTemp += mDn * integralP(mI,i);
 } 
 else //mOctant
  if (mBjorken) {
    mTemp += 0.5* mDn * getTau() * (integralP(mI,0) + integralP(mI,mNSize));
    for (int i=1;i<mNSize;i++) 
	  mTemp += mDn * getTau() * integralP(mI,i);
  }
  else  {
    mTemp += 0.5 * mDn * (integralP(mI,0) + integralP(mI,mNSize));
    for (int i=1;i<mNSize;i++) 
	  mTemp += mDn * integralP(mI,i);
  }
 return mTemp;
}

double CMesh::integralS(int eta) {
	double mTemp=0.;
 if (!mOctant) {
	for (int i=-mXSize;i<=mXSize;i++)
		for (int j=-mYSize;j<=mYSize;j++) 
			mTemp += mDx*mDy*getS(eta,i,j)*getS(eta,i,j,0);
 }
 else { // mOctant
	mTemp = (mDx*mDy/9.)*( getS(eta,0,0) * getS(eta,0,0,0)
					   - getS(eta,0,mYSize) * getS(eta,0,mYSize,0)
					   - getS(eta,mXSize,0) * getS(eta,mXSize,0,0)
					   - 3.*getS(eta,mXSize,mYSize) * getS(eta,mXSize,mYSize,0));
	
	for (int i=1;i<mXSize;i+=2)   
	  mTemp += (mDx*mDy/9.)*( 4.*getS(eta,i,0)*getS(eta,i,0,0) 
						  + 2.*getS(eta,i+1,0)*getS(eta,i+1,0,0)
						  + 4.*getS(eta,i,mYSize)*getS(eta,i,mYSize,0) 
						  + 2.*getS(eta,i+1,mYSize)*getS(eta,i+1,mYSize,0));
	for (int i=1;i<mYSize;i+=2)   
	  mTemp += (mDx*mDy/9.)*( 4.*getS(eta,0,i)*getS(eta,0,i,0) 
						  + 2.*getS(eta,0,i+1)*getS(eta,0,i+1,0)
						  + 4.*getS(eta,mXSize,i)*getS(eta,mXSize,i,0) 
						  + 2.*getS(eta,mXSize,i+1)*getS(eta,mXSize,i+1,0));
	  
	for (int i=1;i<mXSize;i++) 
	  for (int j=1;j<mYSize;j++) 
	    if ((i+j)%2 == 1) 
		  mTemp += (8.*mDx*mDy/9.) * getS(eta,i,j)*getS(eta,i,j,0);
		else if (i%2 == 1)  
		  mTemp +=(16.*mDx*mDy/9.) * getS(eta,i,j)*getS(eta,i,j,0);
		else  
		  mTemp += (4.*mDx*mDy/9.) * getS(eta,i,j)*getS(eta,i,j,0);
 }
 return mTemp;
}

double CMesh::integralS(){
  if (mPureBjorken && mBjorken) return mDn*getTau()*integralS(0);
  else if (mPureBjorken) return mDn*integralS(0);
 
  double mTemp=0.;
 if (!mOctant){
  for (int i=-mNSize;i<=mNSize;i++) 
    if (mBjorken) 
	  mTemp += mDn * getTau() * cosh(mCells[i+mNSizeOrig][0][0]->getEta()) * integralS(i);
	else  
	  mTemp += mDn * integralS(i);
  }
 else //mOctant
  if (mBjorken) {
    mTemp += 0.5* mDn * getTau() * (integralS(0) + integralS(mNSize));
    for (int i=1;i<mNSize;i++) 
	  mTemp += mDn * getTau() * integralS(i);
  }
  else  {
    mTemp += 0.5 * mDn * (integralS(0) + integralS(mNSize));
    for (int i=1;i<mNSize;i++) 
	  mTemp += mDn * integralS(i);
  }
 return mTemp;
}

double CMesh::integralEp(int eta) {
 double mTemp = 0., mTemp2 = 0.;
 if (!mOctant) {
	for (int i=-mXSize;i<=mXSize;i++)
		for (int j=-mYSize;j<=mYSize;j++) 
		  mTemp += mDx*mDy*(getTxyMesh(eta,i,j,2,2) - getTxyMesh(eta,i,j,1,1));
 }
 else { // mOctant
	mTemp = (mDx*mDy/9.) *( (getTxyMesh(eta,0,0,2,2) - getTxyMesh(eta,0,0,1,1))
					   -  (getTxyMesh(eta,0,mYSize,2,2) - getTxyMesh(eta,0,mYSize,1,1))
					   -  (getTxyMesh(eta,mXSize,0,2,2) - getTxyMesh(eta,mXSize,0,1,1))
					   -3*(getTxyMesh(eta,mXSize,mYSize,2,2) - getTxyMesh(eta,mXSize,mYSize,1,1)));
	
	for (int i=1;i<mXSize;i+=2)   
	  mTemp += (mDx*mDy/9.) * ( 4.*(getTxyMesh(eta,i,0,2,2) - getTxyMesh(eta,i,0,1,1) )
							+ 2.*(getTxyMesh(eta,i+1,0,2,2) - getTxyMesh(eta,i+1,0,1,1))
						    + 4.*(getTxyMesh(eta,i,mYSize,2,2) - getTxyMesh(eta,i,mYSize,1,1))
							+ 2.*(getTxyMesh(eta,i+1,mYSize,2,2) - getTxyMesh(eta,i+1,mYSize,1,1)));
							
	for (int i=1;i<mYSize;i+=2)   
	  mTemp += (mDx*mDy/9.) * ( 4.*(getTxyMesh(eta,0,i,2,2) - getTxyMesh(eta,0,i,1,1) )
							+ 2.*(getTxyMesh(eta,0,i+1,2,2) - getTxyMesh(eta,0,i+1,1,1))
						    + 4.*(getTxyMesh(eta,mXSize,i,2,2) - getTxyMesh(eta,mXSize,i,1,1) )
							+ 2.*(getTxyMesh(eta,mXSize,i+1,2,2) - getTxyMesh(eta,mXSize,i+1,1,1)));
	  
	for (int i=1;i<mXSize;i++) 
	  for (int j=1;j<mYSize;j++) 
	    if ( (i+j)%2 == 1) 
		  mTemp += (8.*mDx*mDy/9.) * (getTxyMesh(eta,i,j,2,2) - getTxyMesh(eta,i,j,1,1));
		else if (i%2 == 1)  
		  mTemp +=(16.*mDx*mDy/9.) * (getTxyMesh(eta,i,j,2,2) - getTxyMesh(eta,i,j,1,1));
		else  
		  mTemp += (4.*mDx*mDy/9.) * (getTxyMesh(eta,i,j,2,2) - getTxyMesh(eta,i,j,1,1));
  
	mTemp2 = (mDx*mDy/9.) *( (getTxyMesh(eta,0,0,2,2) + getTxyMesh(eta,0,0,1,1))
					   -  (getTxyMesh(eta,0,mYSize,2,2) + getTxyMesh(eta,0,mYSize,1,1))
					   -  (getTxyMesh(eta,mXSize,0,2,2) + getTxyMesh(eta,mXSize,0,1,1))
					   -3*(getTxyMesh(eta,mXSize,mYSize,2,2) + getTxyMesh(eta,mXSize,mYSize,1,1)));
	
	for (int i=1;i<mXSize;i+=2)   
	  mTemp2 += (mDx*mDy/9.) *( 4.*(getTxyMesh(eta,i,0,2,2) + getTxyMesh(eta,i,0,1,1) )
							+ 2.*(getTxyMesh(eta,i+1,0,2,2) + getTxyMesh(eta,i+1,0,1,1))
						    + 4.*(getTxyMesh(eta,i,mYSize,2,2) + getTxyMesh(eta,i,mYSize,1,1))
							+ 2.*(getTxyMesh(eta,i+1,mYSize,2,2) + getTxyMesh(eta,i+1,mYSize,1,1)));
							
	for (int i=1;i<mYSize;i+=2)   
	  mTemp2 += (mDx*mDy/9.) * ( 4.*(getTxyMesh(eta,0,i,2,2) + getTxyMesh(eta,0,i,1,1) )
							 + 2.*(getTxyMesh(eta,0,i+1,2,2) + getTxyMesh(eta,0,i+1,1,1))
						     + 4.*(getTxyMesh(eta,mXSize,i,2,2) + getTxyMesh(eta,mXSize,i,1,1) )
							 + 2.*(getTxyMesh(eta,mXSize,i+1,2,2) + getTxyMesh(eta,mXSize,i+1,1,1)));
	  
	for (int i=1;i<mXSize;i++) 
	  for (int j=1;j<mYSize;j++) 
	    if ( (i+j)%2 == 1) 
		  mTemp2 += (8.*mDx*mDy/9.) * (getTxyMesh(eta,i,j,2,2) + getTxyMesh(eta,i,j,1,1));
		else if (i%2 == 1)  
		  mTemp2 +=(16.*mDx*mDy/9.) * (getTxyMesh(eta,i,j,2,2) + getTxyMesh(eta,i,j,1,1));
		else  
		  mTemp2 += (4.*mDx*mDy/9.) * (getTxyMesh(eta,i,j,2,2) + getTxyMesh(eta,i,j,1,1));
 }
 
 return -mTemp/mTemp2;
}

double CMesh::integralEp(){
  if (mPureBjorken) 
    if (mBjorken) mDn*getTau()*integralEp(0);
	else return mDn*integralEp(0);
  
 double mTemp=0.;
 if (!mOctant) {
  if (mBjorken) 
    for (int i=-mNSize;i<=mNSize;i++) 
      mTemp += mDn * getTau() * cosh(mCells[i+mNSizeOrig][0][0]->getEta()) * integralEp(i);
  else 
    for (int i=-mNSize;i<=mNSize;i++) 
      mTemp += mDn * integralEp(i);
 } else { //mOctant
  if (mBjorken) 
    for (int i=0;i<=mNSize;i++) 
	  mTemp += mDn * getTau() * integralEp(i);
  else 
    for (int i=0;i<=mNSize;i++) 
	  mTemp += mDn * integralEp(i);
 }
 return mTemp;
}

double CMesh::integralEx(int eta) {
 double mTemp = 0., mTemp2 = 0.;
 if (!mOctant) {
 }
 else { // mOctant
  for (int i=0;i<mXSize;i++)
	for (int j=0;j<mYSize;j++) {
	  mTemp += mDx*mDy*( mDy*mDy*i*i - mDx*mDx*j*j) * getE(eta,i,j);
	  mTemp2+= mDx*mDy*( mDy*mDy*i*i + mDx*mDx*j*j) * getE(eta,i,j);
	}
 /*
	mTemp = (mDx*mDy/9.) *( 0.
					   -  mDy*mDy*mYSize*mYSize * getE(eta,0,mYSize)
					   +  mDx*mDx*mXSize*mXSize * getE(eta,mXSize,0)
					   -3*( mDy*mDy*mYSize*mYSize - mDx*mDx*mXSize*mXSize) getE(eta,mXSize,mYSize));
	
	for (int i=1;i<mXSize;i+=2)   
	  mTemp += (mDx*mDy/9.) * ( 4.* -mDx*mDx*i*i*getE(eta,i,0)
							+ 2.*(getTxyMesh(eta,i+1,0,2,2) - getTxyMesh(eta,i+1,0,1,1))
						    + 4.*(getTxyMesh(eta,i,mYSize,2,2) - getTxyMesh(eta,i,mYSize,1,1))
							+ 2.*(getTxyMesh(eta,i+1,mYSize,2,2) - getTxyMesh(eta,i+1,mYSize,1,1)));
							
	for (int i=1;i<mYSize;i+=2)   
	  mTemp += (mDx*mDy/9.) * ( 4.*(getTxyMesh(eta,0,i,2,2) - getTxyMesh(eta,0,i,1,1) )
							+ 2.*(getTxyMesh(eta,0,i+1,2,2) - getTxyMesh(eta,0,i+1,1,1))
						    + 4.*(getTxyMesh(eta,mXSize,i,2,2) - getTxyMesh(eta,mXSize,i,1,1) )
							+ 2.*(getTxyMesh(eta,mXSize,i+1,2,2) - getTxyMesh(eta,mXSize,i+1,1,1)));
	  
	for (int i=1;i<mXSize;i++) 
	  for (int j=1;j<mYSize;j++) 
	    if ( (i+j)%2 == 1) 
		  mTemp += (8.*mDx*mDy/9.) * (getTxyMesh(eta,i,j,2,2) - getTxyMesh(eta,i,j,1,1));
		else if (i%2 == 1)  
		  mTemp +=(16.*mDx*mDy/9.) * (getTxyMesh(eta,i,j,2,2) - getTxyMesh(eta,i,j,1,1));
		else  
		  mTemp += (4.*mDx*mDy/9.) * (getTxyMesh(eta,i,j,2,2) - getTxyMesh(eta,i,j,1,1));
  
	mTemp2 = (mDx*mDy/9.) *( (getTxyMesh(eta,0,0,2,2) + getTxyMesh(eta,0,0,1,1))
					   -  (getTxyMesh(eta,0,mYSize,2,2) + getTxyMesh(eta,0,mYSize,1,1))
					   -  (getTxyMesh(eta,mXSize,0,2,2) + getTxyMesh(eta,mXSize,0,1,1))
					   -3*(getTxyMesh(eta,mXSize,mYSize,2,2) + getTxyMesh(eta,mXSize,mYSize,1,1)));
	
	for (int i=1;i<mXSize;i+=2)   
	  mTemp2 += (mDx*mDy/9.) *( 4.*(getTxyMesh(eta,i,0,2,2) + getTxyMesh(eta,i,0,1,1) )
							+ 2.*(getTxyMesh(eta,i+1,0,2,2) + getTxyMesh(eta,i+1,0,1,1))
						    + 4.*(getTxyMesh(eta,i,mYSize,2,2) + getTxyMesh(eta,i,mYSize,1,1))
							+ 2.*(getTxyMesh(eta,i+1,mYSize,2,2) + getTxyMesh(eta,i+1,mYSize,1,1)));
							
	for (int i=1;i<mYSize;i+=2)   
	  mTemp2 += (mDx*mDy/9.) * ( 4.*(getTxyMesh(eta,0,i,2,2) + getTxyMesh(eta,0,i,1,1) )
							 + 2.*(getTxyMesh(eta,0,i+1,2,2) + getTxyMesh(eta,0,i+1,1,1))
						     + 4.*(getTxyMesh(eta,mXSize,i,2,2) + getTxyMesh(eta,mXSize,i,1,1) )
							 + 2.*(getTxyMesh(eta,mXSize,i+1,2,2) + getTxyMesh(eta,mXSize,i+1,1,1)));
	  
	for (int i=1;i<mXSize;i++) 
	  for (int j=1;j<mYSize;j++) 
	    if ( (i+j)%2 == 1) 
		  mTemp2 += (8.*mDx*mDy/9.) * (getTxyMesh(eta,i,j,2,2) + getTxyMesh(eta,i,j,1,1));
		else if (i%2 == 1)  
		  mTemp2 +=(16.*mDx*mDy/9.) * (getTxyMesh(eta,i,j,2,2) + getTxyMesh(eta,i,j,1,1));
		else  
		  mTemp2 += (4.*mDx*mDy/9.) * (getTxyMesh(eta,i,j,2,2) + getTxyMesh(eta,i,j,1,1));
*/
 }
 
 return mTemp/mTemp2;
}

double CMesh::integralEx(){
  if (mPureBjorken) 
    if (mBjorken) return mDn*getTau()*integralEx(0);
	else return mDn*integralEx(0);
  
 double mTemp=0.;
 if (!mOctant) {
  if (mBjorken) 
    for (int i=-mNSize;i<=mNSize;i++) 
      mTemp += mDn * getTau() * cosh(mCells[i+mNSizeOrig][0][0]->getEta()) * integralEx(i);
  else 
    for (int i=-mNSize;i<=mNSize;i++) 
      mTemp += mDn * integralEx(i);
 } else { //mOctant
  if (mBjorken) 
    for (int i=0;i<=mNSize;i++) 
	  mTemp += mDn * getTau() * integralEx(i);
  else 
    for (int i=0;i<=mNSize;i++) 
	  mTemp += mDn * integralEx(i);
 }
 return mTemp;
}

double CMesh::averageVr(int eta) {
/*
	double mTemp=0.; double mTemp1=0.;
	for (int i=0;i<mXSize && i<mXSizeOrig;i++)
	  for (int j=0;j<mYSize && j<mYSizeOrig;j++) 
	    if (getT(eta,i,j) > mFOTemp) {
	      mTemp += mDx*mDy*getE(eta,i,j)*getGammaTrans(eta,i,j)*getVr(eta,i,j);
		  mTemp1+= mDx*mDy*getE(eta,i,j)*getGammaTrans(eta,i,j);
	    }
*/
	double mTemp = (mDx*mDy/9.) *(getE(eta,0,0) * getGammaTrans(eta,0,0) * getVr(eta,0,0)
							  - getE(eta,0,mYSize) * getGammaTrans(eta,0,mYSize) * getVr(eta,0,mYSize)
							  - getE(eta,mXSize,0) * getGammaTrans(eta,mXSize,0) * getVr(eta,mXSize,0)
							  - 3.* getE(eta,mXSize,mYSize) * getGammaTrans(eta,mXSize,mYSize) * getVr(eta,mXSize,mYSize));
	
	for (int i=1;i<mXSize;i+=2)  
	  mTemp += (mDx*mDy/9.) * ( 4.*getE(eta,i,0) * getGammaTrans(eta,i,0) * getVr(eta,i,0)
							+ 2.*getE(eta,i+1,0) * getGammaTrans(eta,i+1,0) * getVr(eta,i+1,0)
						    + 4.*getE(eta,i,mYSize) * getGammaTrans(eta,i,mYSize) * getVr(eta,i,mYSize)
							+ 2.*getE(eta,i+1,mYSize) * getGammaTrans(eta,i+1,mYSize) * getVr(eta,i+1,mYSize));
							
	for (int i=1;i<mYSize;i+=2)  
	  mTemp += (mDx*mDy/9.) * ( 4.*getE(eta,0,i) * getGammaTrans(eta,0,i) * getVr(eta,0,i)
	  						+ 2.*getE(eta,0,i+1) * getGammaTrans(eta,0,i+1) * getVr(eta,0,i+1)
							+ 4.*getE(eta,mXSize,i) * getGammaTrans(eta,mXSize,i) * getVr(eta,mXSize,i)
							+ 2.*getE(eta,mXSize,i+1) * getGammaTrans(eta,mXSize,i+1) * getVr(eta,mXSize,i+1));
	  
	for (int i=1;i<mXSize;i++) 
	  for (int j=1;j<mYSize;j++) 
	    if ( (i+j)%2 == 1) 
		  mTemp += (8.*mDx*mDy/9.) * getE(eta,i,j) * getGammaTrans(eta,i,j) * getVr(eta,i,j);
		else if (i%2 == 1)  
		  mTemp +=(16.*mDx*mDy/9.) * getE(eta,i,j) * getGammaTrans(eta,i,j) * getVr(eta,i,j);
		else  
		  mTemp += (4.*mDx*mDy/9.) * getE(eta,i,j) * getGammaTrans(eta,i,j) * getVr(eta,i,j);

	double mTemp1 = (mDx*mDy/9.) *(getE(eta,0,0) * getGammaTrans(eta,0,0) 
							  - getE(eta,0,mYSize) * getGammaTrans(eta,0,mYSize) 
							  - getE(eta,mXSize,0) * getGammaTrans(eta,mXSize,0) 
							  - 3.* getE(eta,mXSize,mYSize) * getGammaTrans(eta,mXSize,mYSize) );
	
	for (int i=1;i<mXSize;i+=2)  
	  mTemp1 += (mDx*mDy/9.) * ( 4.*getE(eta,i,0) * getGammaTrans(eta,i,0) 
							+ 2.*getE(eta,i+1,0) * getGammaTrans(eta,i+1,0) 
						    + 4.*getE(eta,i,mYSize) * getGammaTrans(eta,i,mYSize) 
							+ 2.*getE(eta,i+1,mYSize) * getGammaTrans(eta,i+1,mYSize) );
							
	for (int i=1;i<mYSize;i+=2)  
	  mTemp1 += (mDx*mDy/9.) * ( 4.*getE(eta,0,i) * getGammaTrans(eta,0,i) 
	  						+ 2.*getE(eta,0,i+1) * getGammaTrans(eta,0,i+1)
							+ 4.*getE(eta,mXSize,i) * getGammaTrans(eta,mXSize,i) 
							+ 2.*getE(eta,mXSize,i+1) * getGammaTrans(eta,mXSize,i+1) );
	  
	for (int i=1;i<mXSize;i++) 
	  for (int j=1;j<mYSize;j++) 
	    if ( (i+j)%2 == 1) 
		  mTemp1 += (8.*mDx*mDy/9.) * getE(eta,i,j) * getGammaTrans(eta,i,j);
		else if (i%2 == 1)  
		  mTemp1 +=(16.*mDx*mDy/9.) * getE(eta,i,j) * getGammaTrans(eta,i,j);
		else  
		  mTemp1 += (4.*mDx*mDy/9.) * getE(eta,i,j) * getGammaTrans(eta,i,j);

	return (mTemp/mTemp1);
}

double CMesh::averageVr() {
  if (mPureBjorken) 
    if (mBjorken) return averageVr(0);
	else return averageVr(0);
  
 double mTemp=0.;
 if (!mOctant) {
  if (mBjorken) 
    for (int i=-mNSize;i<=mNSize;i++) 
      mTemp += mDn * getTau() * cosh(mCells[i+mNSizeOrig][0][0]->getEta()) * averageVr(i);
  else 
    for (int i=-mNSize;i<=mNSize;i++) 
      mTemp += mDn * averageVr(i);
 } else { //mOctant
  if (mBjorken) 
    for (int i=0;i<=mNSize;i++) 
	  mTemp += averageVr(i);  //mDn * getTau() * averageVr(i);
  else 
    for (int i=0;i<=mNSize;i++) 
	  mTemp += mDn * averageVr(i);
 }
 return mTemp;
}

double CMesh::getdNd3p(double px, double py, double Y, double m) {
  double mTemp = 0.;
  double pz,T,E,P,p0,pu,puu;
//  p0 = sqrt(m*m + px*py + py*py)*cosh(Y);
//  pz = p0 * tanh(Y);
  
//  printf("%g %g %g %g:\n\n",px,py,Y,m);
  
  for (int eta=0;eta<=mNSize;eta++) {
	p0 = sqrt(m*m + px*py + py*py) * cosh(Y - mDn*eta);
	pz = p0 * tanh(Y - mDn*eta);
	
//	printf("\n%d %g %g:\n\n",eta,p0,pz);
	
    for (int i=0;i<=mXSize;i++)
	  for (int j=0;j<=mYSize;j++) {
	    pu = p0*getS(eta,i,j,0) - px*getS(eta,i,j,1) - py*getS(eta,i,j,2) - pz*getS(eta,i,j,3);

		puu = getPixyMesh(eta,i,j,0,0)*p0*p0;
		puu -= getPixyMesh(eta,i,j,0,1)*p0*px + getPixyMesh(eta,i,j,0,2)*p0*pz + getPixyMesh(eta,i,j,0,3)*p0*pz;
		puu +=  getPixyMesh(eta,i,j,1,1)*px*px +getPixyMesh(eta,i,j,2,2)*py*py +getPixyMesh(eta,i,j,3,3)*pz*pz;
		puu += 2.*(getPixyMesh(eta,i,j,1,2)*px*py +getPixyMesh(eta,i,j,1,3)*px*pz +getPixyMesh(eta,i,j,2,3)*py*pz);
		
		T = getT(eta,i,j);
		E = getE(eta,i,j);
		P = getP(eta,i,j);
//		printf("%d %d:  %g %g %g\n",i,j,exp(-pu/T),0.5*puu/((E+P)*T*T),cosh(Y - mDn*eta));
	    mTemp += p0 * exp( - pu/T) * (1. + 0.5*puu/((E+P)*T*T));// * cosh(Y - mDn*eta);
	  }
	}

  return mTemp*mDx*mDy*mDn*getTau()*JIFPI*JIFPI*JIFPI*8.;
}

double CMesh::getPixyNS(int eta, int x, int y, int xx, int yy) {
  if (mOctant) {
    mCells[eta+1][x+1][y+1]->update();
    return mCells[eta+1][x+1][y+1]->getPixyNS(xx,yy); 
  }
  else {
    mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->update();
    return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getPixyNS(xx,yy);
  }
}

double CMesh::getPixyNSMesh(int eta, int x, int y, int xx, int yy) {
  if (mOctant) {
    mCells[eta+1][x+1][y+1]->update();
    return mCells[eta+1][x+1][y+1]->getPixyNSMesh(xx,yy); 
  }
  else {
    mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->update();
    return mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->getPixyNSMesh(xx,yy);
  }
}

// freezeout surface at eta=0
void CMesh::getFOS(double mFOS[XSIZE+YSIZE][2], double mFOSigma[XSIZE+YSIZE][2], int &foSize) {
 foSize = 0;
 int j = mYSize;  //keep track of where we are in j

 for (int i=0; i<=mXSize; i++) {
  double v;
  if (j+3 <= mYSize) j+=3;
  else j=mYSize;
  
  
   while (j>=0) 
    if (mCells[1][i+1][j]->getTEos() > mFOTemp) {
	  mFOS[foSize][1] =  mCells[1][i+1][j]->getY() + (mCells[1][i+1][j+1]->getY() - mCells[1][i+1][j]->getY())
										* ((mFOTemp - mCells[1][i+1][j]->getTEos())
										/(mCells[1][i+1][j+1]->getTEos() - mCells[1][i+1][j]->getTEos()));
	  mFOS[foSize][0] = mCells[1][i+1][j]->getX();
	  
	  if (i==1) {
	    mFOSigma[foSize-1][0] = 0.;
		mFOSigma[foSize-1][1] = 1.;
	  }
	  else if (i>1) {
	    v = (mFOS[foSize-2][1] - mFOS[foSize][1])/(mFOS[foSize][0] - mFOS[foSize-2][0]);
	    mFOSigma[foSize-1][0] = v/sqrt(1.+v*v);
		mFOSigma[foSize-1][1] = 1./sqrt(1.+v*v);
	  }
  	  foSize++;
	  break;
	}
	else j--;
   
   if (v>1.) break;
 }
 
 j--;
 
 while(j>=0) {
   for (int i=mXSize;i>=0;i--) 
    if (mCells[1][i][j+1]->getTEos() > mFOTemp) {
	  mFOS[foSize][0] =  mCells[1][i][j+1]->getX() + (mCells[1][i+1][j+1]->getX() - mCells[1][i][j+1]->getX())
										* ((mFOTemp - mCells[1][i][j+1]->getTEos())
										/ (mCells[1][i+1][j+1]->getTEos() - mCells[1][i][j+1]->getTEos()));
	  mFOS[foSize][1] = mCells[1][i][j+1]->getY();
	  
	  if (j==0) {
	    mFOSigma[foSize][0] = 1.;
		mFOSigma[foSize][1] = 0.;
	  }
	  
	  double v = (mFOS[foSize-2][0] - mFOS[foSize][0])/(mFOS[foSize][1] - mFOS[foSize-2][1]);
	  mFOSigma[foSize-1][1] = v/sqrt(1.+v*v);
	  mFOSigma[foSize-1][0] = 1./sqrt(1.+v*v);
	  
  	  foSize++;
	  break;
	}
 
   j--;
 }
 
}

// freezeout surface at eta=0
void CMesh::getFOS(double mFOS[XSIZE+YSIZE][2], double mFOSigma[XSIZE+YSIZE][2], 
				   double mFOVelo[XSIZE+YSIZE][2], double mFODiss[XSIZE+YSIZE][5], int &foSize) {
 foSize = 0;
 int j = mYSize;  //keep track of where we are in j

 for (int i=0; i<=mXSize; i++) {
  double v;
  if (j+3 <= mYSize) j+=3;
  else j=mYSize;
  
  
   while (j>=0) 
    if (mCells[1][i+1][j]->getTEos() > mFOTemp) {
	  mFOS[foSize][1] =  mCells[1][i+1][j]->getY() + (mCells[1][i+1][j+1]->getY() - mCells[1][i+1][j]->getY())
										* ((mFOTemp - mCells[1][i+1][j]->getTEos())
										/(mCells[1][i+1][j+1]->getTEos() - mCells[1][i+1][j]->getTEos()));
	  mFOS[foSize][0] = mCells[1][i+1][j]->getX();
	  
	  for (int k=1;k<3;k++) 
	    mFOVelo[foSize][k-1] = mCells[1][i+1][j]->getS(k) + (mCells[1][i+1][j+1]->getS(k) - mCells[1][i+1][j]->getS(k))
										* ((mFOTemp - mCells[1][i+1][j]->getTEos())
										/(mCells[1][i+1][j+1]->getTEos() - mCells[1][i+1][j]->getTEos()));
	  for (int k=5;k<10;k++)
	    mFODiss[foSize][k-5] = mCells[1][i+1][j]->getS(k) + (mCells[1][i+1][j+1]->getS(k) - mCells[1][i+1][j]->getS(k))
										* ((mFOTemp - mCells[1][i+1][j]->getTEos())
										/(mCells[1][i+1][j+1]->getTEos() - mCells[1][i+1][j]->getTEos()));
	  
	  if (i==1) {
	    mFOSigma[foSize-1][0] = 0.;
		mFOSigma[foSize-1][1] = 1.;
	  }
	  else if (i>1) {
	    v = (mFOS[foSize-2][1] - mFOS[foSize][1])/(mFOS[foSize][0] - mFOS[foSize-2][0]);
	    mFOSigma[foSize-1][0] = v/sqrt(1.+v*v);
		mFOSigma[foSize-1][1] = 1./sqrt(1.+v*v);
	  }
  	  foSize++;
	  break;
	}
	else j--;
   
   if (v>1.) break;
 }
 
 j--;
 
 while(j>=0) {
   for (int i=mXSize;i>=0;i--) 
    if (mCells[1][i][j+1]->getTEos() > mFOTemp) {
	  mFOS[foSize][0] =  mCells[1][i][j+1]->getX() + (mCells[1][i+1][j+1]->getX() - mCells[1][i][j+1]->getX())
										* ((mFOTemp - mCells[1][i][j+1]->getTEos())
										/ (mCells[1][i+1][j+1]->getTEos() - mCells[1][i][j+1]->getTEos()));
	  mFOS[foSize][1] = mCells[1][i][j+1]->getY();
	  
	  for (int k=1;k<3;k++)
	    mFOVelo[foSize][k-1] = mCells[1][i][j+1]->getS(k) + (mCells[1][i+1][j+1]->getS(k) - mCells[1][i][j+1]->getS(k))
										* ((mFOTemp - mCells[1][i][j+1]->getTEos())
										/ (mCells[1][i+1][j+1]->getTEos() - mCells[1][i][j+1]->getTEos()));
	  for (int k=5;k<10;k++)
	    mFODiss[foSize][k-5] = mCells[1][i][j+1]->getS(k) + (mCells[1][i+1][j+1]->getS(k) - mCells[1][i][j+1]->getS(k))
										* ((mFOTemp - mCells[1][i][j+1]->getTEos())
										/ (mCells[1][i+1][j+1]->getTEos() - mCells[1][i][j+1]->getTEos()));
	  
	  if (j==0) {
	    mFOSigma[foSize][0] = 1.;
		mFOSigma[foSize][1] = 0.;
	  }
	  
	  double v = (mFOS[foSize-2][0] - mFOS[foSize][0])/(mFOS[foSize][1] - mFOS[foSize-2][1]);
	  mFOSigma[foSize-1][1] = v/sqrt(1.+v*v);
	  mFOSigma[foSize-1][0] = 1./sqrt(1.+v*v);
	  
  	  foSize++;
	  break;
	}
 
   j--;
 }
 
}

double CMesh::getFOP() {
	for (int i=mXSize; i> 0; i--) 
	  if ( mCells[1][i][1]->getTEos() > mFOTemp) 
	    return mCells[1][i][1]->getX() + (mCells[1][i+1][1]->getX() - mCells[1][i][1]->getX())
										* (( mFOTemp - mCells[1][i][1]->getTEos() )
										/(mCells[1][i+1][1]->getTEos() - mCells[1][i][1]->getTEos()));
}

double CMesh::getFOSSST() {
	for (int i=mXSize; i> 0; i--) 
	  if ( mCells[1][i][1]->getTEos() > mFOTemp) {
		double v1 =  ( mCells[1][i][1]->getPixyMeshEos(1,1) + mCells[1][i][1]->getPixyMeshEos(2,2))
					/( mCells[1][i][1]->getE() + mCells[1][i][1]->getPEos());
		double v2 =  ( mCells[1][i+1][1]->getPixyMeshEos(1,1) + mCells[1][i+1][1]->getPixyMeshEos(2,2))
					/( mCells[1][i+1][1]->getE() + mCells[1][i+1][1]->getPEos());
	    return (1./JHBARC)*(v1 + (v2- v1) * ((mFOTemp - mCells[1][i][1]->getTEos()) /(mCells[1][i+1][1]->getTEos() - mCells[1][i][1]->getTEos())));
	  }
}

double CMesh::getELoss(CMesh* mMesh, int eta) {
	double mTemp=0.; 
 if (!mOctant) {
	for (int i=-mXSize;i<=mXSize;i++)
		for (int j=-mYSize;j<=mYSize;j++) 
			mTemp += mDx*mDy*getS(eta,i,j)*getS(eta,i,j,0);
 }
 else { // mOctant
 
   if (mBjorken) {
	mTemp = (mDx*mDy/9.)*(  getTxyMesh(eta,0,0,3,3) 
						- getTxyMesh(eta,0,mYSize,3,3) 
					    - getTxyMesh(eta,mXSize,0,3,3) 
						- 3.*getTxyMesh(eta,mXSize,mYSize,3,3));
	
	for (int i=1;i<mXSize;i+=2)   
	  mTemp += (mDx*mDy/9.)*( 4.*getTxyMesh(eta,i,0,3,3) 
						  + 2.*getTxyMesh(eta,i+1,0,3,3) 
						  + 4.*getTxyMesh(eta,i,mYSize,3,3) 
						  + 2.*getTxyMesh(eta,i+1,mYSize,3,3));
	for (int i=1;i<mYSize;i+=2)   
	  mTemp += (mDx*mDy/9.)*( 4.*getTxyMesh(eta,0,i,3,3) 
						  + 2.*getTxyMesh(eta,0,i+1,3,3) 
						  + 4.*getTxyMesh(eta,mXSize,i,3,3) 
						  + 2.*getTxyMesh(eta,mXSize,i+1,3,3));
	  
	for (int i=1;i<mXSize;i++) 
	  for (int j=1;j<mYSize;j++) 
	    if ((i+j)%2 == 1) 
		  mTemp += (8.*mDx*mDy/9.) * (getTxyMesh(eta,i,j,3,3));
		else if (i%2 == 1)  
		  mTemp +=(16.*mDx*mDy/9.) * (getTxyMesh(eta,i,j,3,3));
		else  
		  mTemp += (4.*mDx*mDy/9.) * (getTxyMesh(eta,i,j,3,3));

	  mTemp /= getTau();
   }
    //outer radial boundaries
   for (int j=0;j<=mYSize;j++) 
	 mTemp += mDy*mMesh->getTxyMesh(eta,mXSize,j,0,1);
   mTemp -= 0.5*mDy*( mMesh->getTxyMesh(eta,mXSize,0,0,1) + mMesh->getTxyMesh(eta,mXSize,mYSize,0,1));
   for (int j=0;j<=mXSize;j++) 
	 mTemp += mDx*mMesh->getTxyMesh(eta,j,mYSize,0,2);
   mTemp -= 0.5*mDx*(mMesh->getTxyMesh(eta,0,mYSize,0,2) + mMesh->getTxyMesh(eta,mXSize,mYSize,0,2));
   
  }
  
 return mTemp;
}

double CMesh::getELoss(CMesh* mMesh) {
  if (mPureBjorken && mBjorken) 
	return mDn*getTau()*getELoss(mMesh,0);
  else if (mPureBjorken) 
	return mDn*getELoss(mMesh,0);
 
  double mTemp=0.;
 if (!mOctant){
  for (int i=-mNSize;i<=mNSize;i++) 
    if (mBjorken) 
	  mTemp += mDn * getTau() * cosh(mCells[i+mNSizeOrig][0][0]->getEta()) * getELoss(mMesh,i);
	else  
	  mTemp += mDn * getELoss(mMesh,i);
  }
 else //mOctant
  if (mBjorken) {
    mTemp += 0.5* mDn * getTau() * (getELoss(mMesh,0) + getELoss(mMesh,mNSize));
    for (int i=1;i<mNSize;i++) 
	  mTemp += mDn * getTau() * getELoss(mMesh,i);
	  
	for (int i=0;i<mXSize;i++)
     for (int j=0;j<mYSize;j++)
	   mTemp += mDx*mDy*getTxyMesh(mNSize,i,j,0,3);
  }
  else  {
    mTemp += 0.5 * mDn * (getELoss(mMesh,0) + getELoss(mMesh,mNSize));
    for (int i=1;i<mNSize;i++) 
	  mTemp += mDn * getELoss(mMesh,i);
  }
 return mTemp;

}

double CMesh::getPLoss(CMesh* mMesh, int mI, int eta) {
	double mTemp=0.; 
 if (!mOctant) {
	for (int i=-mXSize;i<=mXSize;i++)
		for (int j=-mYSize;j<=mYSize;j++) 
			mTemp += mDx*mDy*getS(eta,i,j)*getS(eta,i,j,0);
 }
 else { // mOctant
 
   if ((mBjorken) && mI==3 ) {
	mTemp = (mDx*mDy/9.)*(  getTxyMesh(eta,0,0,0,mI) 
						- getTxyMesh(eta,0,mYSize,0,mI) 
					    - getTxyMesh(eta,mXSize,0,0,mI) 
						- 3.*getTxyMesh(eta,mXSize,mYSize,0,mI));
	
	for (int i=1;i<mXSize;i+=2)   
	  mTemp += (mDx*mDy/9.)*( 4.*getTxyMesh(eta,i,0,0,mI) 
						  + 2.*getTxyMesh(eta,i+1,0,0,mI) 
						  + 4.*getTxyMesh(eta,i,mYSize,0,mI) 
						  + 2.*getTxyMesh(eta,i+1,mYSize,0,mI));
	for (int i=1;i<mYSize;i+=2)   
	  mTemp += (mDx*mDy/9.)*( 4.*getTxyMesh(eta,0,i,0,mI) 
						  + 2.*getTxyMesh(eta,0,i+1,0,mI) 
						  + 4.*getTxyMesh(eta,mXSize,i,0,mI) 
						  + 2.*getTxyMesh(eta,mXSize,i+1,0,mI));
	  
	for (int i=1;i<mXSize;i++) 
	  for (int j=1;j<mYSize;j++) 
	    if ((i+j)%2 == 1) 
		  mTemp += (8.*mDx*mDy/9.) * (getTxyMesh(eta,i,j,0,mI));
		else if (i%2 == 1)  
		  mTemp +=(16.*mDx*mDy/9.) * (getTxyMesh(eta,i,j,0,mI));
		else  
		  mTemp += (4.*mDx*mDy/9.) * (getTxyMesh(eta,i,j,0,mI));

	  mTemp /= getTau();
   }
   
    //outer radial boundaries
   for (int j=0;j<=mYSize;j++) 
	 mTemp += mDy*mMesh->getTxyMesh(eta,mXSize,j,mI,1);
   mTemp -= 0.5*mDy*( mMesh->getTxyMesh(eta,mXSize,0,mI,1) + mMesh->getTxyMesh(eta,mXSize,mYSize,mI,1));
   for (int j=0;j<=mXSize;j++) 
	 mTemp += mDx*mMesh->getTxyMesh(eta,j,mYSize,mI,2);
   mTemp -= 0.5*mDx*(mMesh->getTxyMesh(eta,0,mYSize,mI,2) + mMesh->getTxyMesh(eta,mXSize,mYSize,mI,2));
   
  }
  
 return mTemp;
}

double CMesh::getPLoss(CMesh* mMesh, int mI) {
  if (mPureBjorken) 
    if(mBjorken) 
	  return mDn*getTau()*getPLoss(mMesh,mI,0);
    else 
	  return mDn*getPLoss(mMesh,mI,0);
 
 double mTemp=0.;
 if (!mOctant){
  for (int i=-mNSize;i<=mNSize;i++) 
    if (mBjorken) 
	  mTemp += mDn * getTau() * cosh(mCells[i+mNSizeOrig][0][0]->getEta()) * getPLoss(mMesh,mI,i);
	else  
	  mTemp += mDn * getPLoss(mMesh,mI,i);
  }
 else //mOctant
  if (mBjorken) {
    mTemp += 0.5* mDn * getTau() * (getPLoss(mMesh,mI,0) + getPLoss(mMesh,mI,mNSize));
    for (int i=1;i<mNSize;i++) 
	  mTemp += mDn * getTau() * getPLoss(mMesh,mI,i);
	  
	for (int i=0;i<mXSize;i++)
     for (int j=0;j<mYSize;j++)
	   mTemp += mDx*mDy*getTxyMesh(mNSize,i,j,mI,3);
  }
  else  {
    mTemp += 0.5 * mDn * (getPLoss(mMesh,mI,0) + getPLoss(mMesh,mI,mNSize));
    for (int i=1;i<mNSize;i++) 
	  mTemp += mDn * getPLoss(mMesh,mI,i);
  }
 return mTemp;

}

//quadratic extrapolation from three evenly spaced cells to their neighbor
//used to smoothen grid edges
void CMesh::smoothEdge(CCell* c1,CCell* c2, CCell* c3, CCell* cOut){

  // quadratic fits to velocity
  for (int i=1;i<4;i++)	 
	if (c3->getS(i) == 0.) cOut->setS(i,0.);
      else cOut->setS(i, c1->getS(i) - 3.*(c2->getS(i) - c3->getS(i)));
//  	  else cOut->setS(i, - c2->getS(i) + 2.*c3->getS(i));
  
  // quadratic fits to scaled momentum anisotropies (a_i/alpha)
  for (int i=5;i<11;i++) 
    cOut->setS(i, c1->getS(i) - 3.*c2->getS(i) + 3.*c3->getS(i));
//	cOut->setS(i, - c2->getS(i) + 2.*c3->getS(i));
	
  //logarithmic smoothing for e
  
  // forces exponential tail
  cOut->setE( c3->getE() * (c3->getE()/c2->getE()));
  
  // allows exponential and gaussian components -- *******seems unstable to tail curling*******
//  cOut->setE(c1->getE()*pow( (c3->getE()/c2->getE()) ,3.));

  // forces gaussian tail
  /*
  if (c3->getX() - c2->getX() != 0.)
    cOut->setE( c2->getE() * pow( (c3->getE()/c2->getE()), 
				4.*c3->getX()/(c3->getX() + c2->getX())));
  else if (c3->getY() - c2->getY() != 0.)
    cOut->setE( c2->getE() * pow( (c3->getE()/c2->getE()), 
				4.*c3->getY()/(c3->getY() + c2->getY())));
  else if (c3->getEta() - c2->getEta() != 0.)
    cOut->setE( c2->getE() * pow( (c3->getE()/c2->getE()), 
				4.*c3->getEta()/(c3->getEta() + c2->getEta())));
*/
////******** RELICS **********//////
  /*
    if (c3->getS(i) == 0. || c2->getS(i) == 0.) 
	  cOut->setS(i,0.);
    else if (fabs(c3->getS(i)) < 1E-6)
//	  cOut->setS(i, 2.*c3->getS(i) - c2->getS(i));
	  cOut->setS(i, c1->getS(i) - 3.*c2->getS(i) + 3.*c3->getS(i));
	else  
	  cOut->setS(i,c1->getS(i)*pow(c3->getS(i)/c2->getS(i),3.));
  /*
  for (int i=1;i<11;i++) 
  if (c3->getX(1) == 10.2 && c3->getX(2) == 10.2 && c3->getX(3) == 0. && fabs(c3->getS(i)) < 1E-6 && cOut->getS(i) != 0.)
      printf("%d:\n %0.9g %0.9g %0.9g %0.9g\n %0.9g %0.9g\n %0.9g %0.9g\n",
				i,
				c1->getS(i),c2->getS(i),c3->getS(i),cOut->getS(i),
				c1->getS(i)*pow(c3->getS(i)/c2->getS(i),3.),
				cOut->getS(i) - c1->getS(i)*pow(c3->getS(i)/c2->getS(i),3.),
				cOut->getS(i)*pow(c2->getS(i)/c3->getS(i),3.), 
				c1->getS(i) - cOut->getS(i)*pow(c2->getS(i)/c3->getS(i),3.));
  
  /*
//  for (int i=5;i<11;i++) if (c3->getS(i) == 0.) cOut->setS(i,0.);
//  for (int i=5;i<11;i++) if (c2->getS(i) == c3->getS(i)) cOut->setS(i,c2->getS(i));
  
  cOut->setS(5, c1->getS(5) - 3.*c2->getS(5) + 3.*c3->getS(5));
  
//  if (SVRATIO != 0. && mBjorken && fabs(c2->getS(6)) > 1E-10 && fabs(c3->getS(6)) > 1E-10) 
//    cOut->setS(6, c1->getS(6)*pow( (c3->getS(6)/c2->getS(6)),3.));
//  else 
    cOut->setS(6, c1->getS(6) - 3.*c2->getS(6) + 3.*c3->getS(6));
  
  for (int i=7;i<10;i++) {
    cOut->setS(i, c1->getS(i) - 3.*c2->getS(i) + 3.*c3->getS(i));
  }
	
//  if (BVRATIO != 0. && mBjorken && fabs(c2->getS(10)) > 1E-10 && fabs(c3->getS(10)) > 1E-10) 
//    cOut->setS(10, c1->getS(10)*pow( (c3->getS(10)/c2->getS(10)),3.));	
//  else 
    cOut->setS(10, c1->getS(10) - 3.*c2->getS(10) + 3.*c3->getS(10));
  
  */
}

void CMesh::flipEdge() {
	double mTau = getTau();
	
	if (!mPureBjorken) {
		for (int i=0;i<=mXSize;i++)
			for (int j=0;j<=mYSize;j++) {	
				flipEdge(3,mCells[2][i+1][j+1],
						 mCells[1][i+1][j+1],
						 mCells[0][i+1][j+1]);
				mCells[0][i+1][j+1]->setTau(mTau);
			}
			//inner rad bound
		for (int i=-1;i<=mNSize;i++) {
			for (int j=0;j<=mYSize;j++) {
				flipEdge(1,mCells[i+1][2][j+1],
						 mCells[i+1][1][j+1],
						 mCells[i+1][0][j+1]);
				mCells[i+1][0][j+1]->setTau(mTau);
			}
			for (int j=-1;j<=mXSize;j++) {	
				flipEdge(2,mCells[i+1][j+1][2],
						 mCells[i+1][j+1][1],
						 mCells[i+1][j+1][0]);
				mCells[i+1][j+1][0]->setTau(mTau);
			}
		}
	}
	else {
			// inner boundaries
		for (int j=0;j<=mYSize;j++) {
			flipEdge(1,mCells[1][2][j+1],
					 mCells[1][1][j+1],
					 mCells[1][0][j+1]);
			mCells[1][0][j+1]->setTau(mTau);
		}
		for (int j=-1;j<=mXSize;j++) {	
			flipEdge(2,mCells[1][j+1][2],
					 mCells[1][j+1][1],
					 mCells[1][j+1][0]);
			mCells[1][j+1][0]->setTau(mTau);
		}	
	}
}

void CMesh::flipEdge(int fE, CCell* c1, CCell* c2, CCell* cOut) {
/*
  if (INIT_EXP_X && fE == 1) {
    for (int i=1;i<4;i++) cOut->setS(i,c1->getS(i));
	cOut->setE( pow(c2->getE(),2)/c1->getE());
	if (c1->getS(5) != 0.) 
	  cOut->setS(5, c2->getS(5)*c2->getS(5)/c1->getS(5));
	else cOut->setS(5, 0.);
	if (c1->getS(6) != 0.) 
	  cOut->setS(6, c2->getS(6)*c2->getS(6)/c1->getS(6));
	else cOut->setS(6, 0.);
  }	
  else if (INIT_EXP_Y && fE == 2) {
    for (int i=1;i<4;i++) cOut->setS(i,c1->getS(i));
	cOut->setE( pow(c2->getE(),2)/c1->getE());
  } 
  else if (INIT_EXP_Z && fE == 3) {
	for (int i=1;i<4;i++) cOut->setS(i,c1->getS(i));
	cOut->setE( pow(c2->getE(),2)/c1->getE());
  }
  else {
*/
	for (int i=1;i<4;i++)	 { 
		if (i == fE) 
			cOut->setS(i, -c1->getS(i));
		else 
			cOut->setS(i,  c1->getS(i));
    }
    cOut->setE(c1->getE());
    cOut->setS(5, c1->getS(5));
    cOut->setS(6, c1->getS(6));
		//  }
	
	if (fE == 1) {
		cOut->setS(7,-c1->getS(7));
		cOut->setS(8,-c1->getS(8));
		cOut->setS(9, c1->getS(9));
	} 
	else if (fE == 2) {
		cOut->setS(7,-c1->getS(7));
		cOut->setS(8, c1->getS(8));
		cOut->setS(9,-c1->getS(9));
	} 
	else if (fE == 3) {
		cOut->setS(7, c1->getS(7));
		cOut->setS(8,-c1->getS(8));
		cOut->setS(9,-c1->getS(9));
	}
	
	cOut->setS(10,c1->getS(10));
	cOut->setActive(c1->getActive());
}

void CMesh::deaden() {
	for (int i=mNSize; i>=0; i--) 
		for (int j=mXSize; j>=0; j--) 
			for (int k=mYSize; k>0; k--) {
				 		  
				if (!getActive(i,j,k)) continue;
				if (getT(i,j,k) > mDeadT) break;
				else {
					deactivate(i,j,k);
					sLoss += mDn*mDx*mDy * getTau() * getS(i,j,k)*getS(i,j,k,0);					
					eLoss += mDn*mDx*mDy * getTau() * getTxyMesh(i,j,k,0,0);
				}
			}
	
	for (int i=mNSize; i>0; i--) 
		if (!getActive(i,0,0)) mNSize = i;
		else break;
	
	for (int i=mXSize; i>0; i--) 
		if (!getActive(0,i,0)) mXSize = i;
		else break;
	
	for (int i=mYSize; i>0; i--) 
		if (!getActive(0,0,i)) mYSize = i;
		else break;
}

bool CMesh::detectCrash(){
  for (int i=-1;i<=mNSize;i++)
    for (int j=-1;j<=mXSize;j++)
	  for (int k=-1;k<=mYSize;k++)
	    for (int l=1;l<11;l++)
		  if ( mCells[i+1][j+1][k+1]->getS(l) != mCells[i+1][j+1][k+1]->getS(l)) {
		    printf("\n\ncrash cell (t=%0.6g) : n=%d x=%d y=%d s(%d) = %0.6g\n",
					getTau(),i,j,k,l, mCells[i+1][j+1][k+1]->getS(l));
			mCells[i+1][j+1][k+1]->print();		
		    return false;
		  }

  return true;
}

void CMesh::initNS() {
	if (!mPureBjorken) 
		for (int i=0;i<=mNSize;i++)
			for (int j=0;j<=mXSize;j++)
				for (int k=0;k<=mYSize;k++) {
					if (!getActive(i,j,k)) break;
					mCells[i+1][j+1][k+1]->initNS();
				}
	else 
		for (int j=0;j<=mXSize;j++)
			for (int k=0;k<=mYSize;k++) {
				if (!getActive(0,j,k)) break;
				mCells[1][j+1][k+1]->initNS();				
			}
	
	flipEdge();	
}

void CMesh::printCell(int eta,int x,int y) {
 if (!mOctant)
	mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->print();
 else
	mCells[eta+1][x+1][y+1]->print();
}

void CMesh::printCell(int eta,int x,int y,int opt) {
 if (!mOctant)
	mCells[eta+mNSizeOrig][x+mXSizeOrig][y+mYSizeOrig]->print(opt);
 else
	mCells[eta+1][x+1][y+1]->print(opt);
}

double CMesh::getTZR(int eta, int x, int y) {
	return 0.;
}

double CMesh::getTZPhi(int eta, int x, int y) {
	return 0.;
}

double CMesh::getTRR(int eta, int x, int y) {
	if (x!=0 || y!=0)
 	  return (mCells[eta+1][x+1][y+1]->getTxyMeshEos(1,1)*mDx*mDx*(double)(x*x)
		  +   mCells[eta+1][x+1][y+1]->getTxyMeshEos(2,2)*mDy*mDy*(double)(y*y)
		  +  2.*mDx*mDy*((double)(x*y))*mCells[eta+1][x+1][y+1]->getTxyMeshEos(1,2))/(mDx*mDx*(double)(x*x)+mDy*mDy*(double)(y*y));
	else return 0.;
}

double CMesh::getTRPhi(int eta, int x, int y) {
	if (x!=0 || y!=0)
	  return (mDx*mDy*x*y*(mCells[eta+1][x+1][y+1]->getTxyMeshEos(1,1)
						-mCells[eta+1][x+1][y+1]->getTxyMeshEos(2,2))
		  +  (mDx*mDx*x*x - mDy*mDy*y*y)*mCells[eta+1][x+1][y+1]->getTxyMeshEos(1,2))/(mDx*mDx*x*x+mDy*mDy*y*y);
	else return 0.;
}

double CMesh::getTPhiPhi(int eta, int x, int y) {
	if (x!=0 || y!=0)
	  return (mDx*mDx*x*x*mCells[eta+1][x+1][y+1]->getTxyMeshEos(2,2)
		  -  mDy*mDy*y*y*mCells[eta+1][x+1][y+1]->getTxyMeshEos(1,1)
		  +  2.*mDx*mDy*x*y*mCells[eta+1][x+1][y+1]->getTxyMeshEos(1,2))/(mDx*mDx*x*x+mDy*mDy*y*y);
	else return 0.;
}

double CMesh::getPiZR(int eta, int x, int y) {
	return 0.;
}

double CMesh::getPiZPhi(int eta, int x, int y) {
	return 0.;
}

double CMesh::getPiRR(int eta, int x, int y) {
	if (x!=0 || y!=0)
 	  return (mCells[eta+1][x+1][y+1]->getPixyMeshEos(1,1)*mDx*mDx*(double)(x*x)
		  +  mCells[eta+1][x+1][y+1]->getPixyMeshEos(2,2)*mDy*mDy*(double)(y*y)
		  +  2.*mDx*mDy*((double)(x*y))*mCells[eta+1][x+1][y+1]->getPixyMeshEos(1,2))/(mDx*mDx*(double)(x*x)+mDy*mDy*(double)(y*y));
	else return 1.;
}

double CMesh::getPiRPhi(int eta, int x, int y) {
	if (x!=0 || y!=0)
	  return (mDx*mDy*x*y*(mCells[eta+1][x+1][y+1]->getPixyMeshEos(1,1)
						-mCells[eta+1][x+1][y+1]->getPixyMeshEos(2,2))
		  +  (mDx*mDx*x*x - mDy*mDy*y*y)*mCells[eta+1][x+1][y+1]->getPixyMeshEos(1,2))/(mDx*mDx*x*x+mDy*mDy*y*y);
	else return 0.;
}

double CMesh::getPiPhiPhi(int eta, int x, int y) {
	if (x!=0 || y!=0)
	  return (mDx*mDx*x*x*mCells[eta+1][x+1][y+1]->getPixyMeshEos(2,2)
		  -  mDy*mDy*y*y*mCells[eta+1][x+1][y+1]->getPixyMeshEos(1,1)
		  +  2.*mDx*mDy*x*y*mCells[eta+1][x+1][y+1]->getPixyMeshEos(1,2))/(mDx*mDx*x*x+mDy*mDy*y*y);
	else return 0.;
}

void CMesh::checkAzimuthalSymmetry() {
	double tolerance = 1E-10;
	
	for (int i=0;i<=mNSize;i++) 
		for (int j=0;j<=mXSize;j++)
			for (int k=j;k<=mYSize;k++) {
				if (!getActive(i,j,k)) break;
				 				if ( abs(getS(i,j,k,1) - getS(i,k,j,2)) > tolerance*abs(getS(i,j,k,1)))
									printf("velo sym prob at %d %d %d (%0.6g != %0.6g)\n",i,j,k,getS(i,j,k,1),getS(i,k,j,2));
				if ( abs(getS(i,j,k,2) - getS(i,k,j,1)) > tolerance*abs(getS(i,j,k,2)))
					printf("velo sym prob at %d %d %d (%0.6g != %0.6g)\n",i,j,k,getS(i,j,k,2),getS(i,k,j,1));
				if ( abs(getS(i,j,k,4) - getS(i,k,j,4)) > tolerance*abs(getS(i,j,k,4)))
					printf("energy sym prob at %d %d %d (%0.6g != %0.6g)\n",i,j,k,getS(i,j,k,4),getS(i,k,j,4));
				
				if ( abs(getS(i,j,k,5) + getS(i,k,j,5)) > tolerance*abs(getS(i,j,k,5)) && k!=j)
					printf("a1 sym prob at %d %d %d (%0.6g != %0.6g)\n",i,j,k,getS(i,j,k,5),-getS(i,k,j,5));
				 				
				if ( abs(getS(i,j,k,6) - getS(i,k,j,6)) > tolerance*abs(getS(i,j,k,6)))
					printf("a2 sym prob at %d %d %d (%0.6g != %0.6g)\n",i,j,k,getS(i,j,k,6),getS(i,k,j,6));
				if ( abs(getS(i,j,k,7) - getS(i,k,j,7)) > tolerance*abs(getS(i,j,k,7)))
					printf("a3 sym prob at %d %d %d (%0.6g != %0.6g)\n",i,j,k,getS(i,j,k,7),getS(i,k,j,7));
				if ( abs(getS(i,j,k,8) - getS(i,k,j,9)) > tolerance*abs(getS(i,j,k,8)))
					printf("a4/a5 sym prob at %d %d %d (%0.6g != %0.6g)\n",i,j,k,getS(i,j,k,8),getS(i,k,j,9));
			}
}

//statics 
parameterMap* CMesh::pMap;
bool CMesh::mOctant, CMesh::mPureBjorken, CMesh::mBjorken, CMesh::mPrintMs, CMesh::mSVTrimInit, CMesh::mSVTrim, CMesh::mSawtooth;
double CMesh::mFOTemp, CMesh::mDeadT, CMesh::mT0, CMesh::mDx, CMesh::mDy, CMesh::mDn;
int CMesh::mNSize, CMesh::mXSize, CMesh::mYSize, CMesh::mNSizeOrig, CMesh::mXSizeOrig, CMesh::mYSizeOrig;
int CMesh::mPrintN, CMesh::mPrintX, CMesh::mPrintY;
double CMesh::wnRho0, CMesh::wnRAu, CMesh::wnXi, CMesh::wnSigma, CMesh::wnA, CMesh::wnB, CMesh::wnKTau;
double CMesh::mE0, CMesh::mRn, CMesh::mRt, CMesh::mRa, CMesh::mRx, CMesh::mRy, CMesh::mRSig, CMesh::mWnBinRatio;
double CMesh::mInitFlow, CMesh::mInitNS;
