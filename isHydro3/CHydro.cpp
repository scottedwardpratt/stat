#include "CHydro.h"

CHydro::CHydro(parameterMap* pM) {
	pMap = pM;

	mDataRoot = parameter::getS(*pMap,"HYDRO_OUTPUT_DATAROOT",mDataRoot);
	
	int sLength = mDataRoot.size();
	if (mDataRoot[sLength-1] != '/') mDataRoot.append("/");
  
	mDebug = parameter::getB(*pMap,"HYDRO_DEBUG",false);
	mIoIntegrals = parameter::getB(*pMap,"HYDRO_IO_INTEGRALS", false);
	mIoSlices = parameter::getB(*pMap,"HYDRO_IO_SLICES",false);
	mIoSpots = parameter::getB(*pMap,"HYDRO_IO_SPOTS",false);
	mIoTechqm = parameter::getB(*pMap,"HYDRO_IO_TECHQM",false);
	mIoOscarFull = parameter::getB(*pMap,"HYDRO_IO_OSCARFULL",false);
	mIoOscarHyper = parameter::getB(*pMap,"HYDRO_IO_OSCARHYPER",false);
	mIoFull = parameter::getB(*pMap,"HYDRO_IO_FULL",false);
	mIoSpectra = parameter::getB(*pMap,"HYDRO_IO_SPECTRA",false);

	mSawtooth = parameter::getB(*pMap,"HYDRO_SAWTOOTH",false);
	mRK2  = parameter::getB(*pMap,"HYDRO_RK2",true);
	mRK4  = parameter::getB(*pMap,"HYDRO_Rk4",false);
	mLinT = parameter::getB(*pMap,"HYDRO_LINT",false);
	mLogT = parameter::getB(*pMap,"HYDRO_LOGT",false);
	mLogSinhT = parameter::getB(*pMap,"HYDRO_LOGSINHT",false);

	mTStep = parameter::getI(*pMap,"HYDRO_TSTEP",100);
	mXSize = parameter::getI(*pMap,"HYDRO_XSIZE",60);
	mYSize = parameter::getI(*pMap,"HYDRO_YSIZE",60);
	mNSize = parameter::getI(*pMap,"HYDRO_NSIZE",20);
	mIoSliceTStep = parameter::getI(*pMap,"HYDRO_IO_SLICE_TSTEP",1);
	mIoTechqmTStep = parameter::getI(*pMap,"HYDRO_IO_TECHQM_TSTEP",1);
	mIoOscarTStep = parameter::getI(*pMap,"HYDRO_IO_OSCAR_TSTEP",1);
	mPrintStep = parameter::getI(*pMap,"HYDRO_PRINTSTEP",(mTStep/5));
	mPrintX = parameter::getI(*pMap,"HYDRO_PRINTX",0);
	mPrintY = parameter::getI(*pMap,"HYDRO_PRINTY",0);
	mPrintN = parameter::getI(*pMap,"HYDRO_PRINTN",0);

	mBjorken = parameter::getB(*pMap,"HYDRO_BJORKEN",true);
	mPureBjorken = parameter::getB(*pMap,"HYDRO_PURE_BJORKEN",false);
	mHalving = parameter::getB(*pMap, "HYDRO_HALVING", false);
	mOctant = parameter::getB(*pMap, "HYDRO_OCTANT", false);

	mInitNS = parameter::getD(*pMap,"HYDRO_INIT_NS",0.);
	mInitFlow = parameter::getD(*pMap,"HYDRO_INIT_FLOW",0.);

	mOldFile = parameter::getB(*pMap,"HYDRO_OLD_FILE",false);
	mOldFileName = parameter::getS(*pMap,"HYDRO_OLD_FILENAME","");

	mFoTemp = parameter::getD(*pMap,"HYDRO_FOTEMP",0.13);
	mT0 = parameter::getD(*pMap,"HYDRO_T0",1.0);
	mDt = parameter::getD(*pMap,"HYDRO_DT",0.03);
	mE0 = parameter::getD(*pMap,"HYDRO_E0",1.0);
  
}

int CHydro::runHydro() {
	time(&start);

	if (mIoIntegrals)	openFileIntegrals();
	if (mIoSpots)		openFileSpots();
	if (mIoTechqm) {
		string temp = mDataRoot;
		temp.append("fTQMX.txt");
		fTQMX = fopen(temp.c_str(),"w");
		
		temp = mDataRoot;
		temp.append("fTQMY.txt");
		fTQMY = fopen(temp.c_str(),"w");
		
		temp = mDataRoot;
		temp.append("fTQMR.txt");
		fTQMR = fopen(temp.c_str(),"w");
	}
  
	if (mOldFile) {
		onMesh   = new CMesh(mOldFileName.c_str());
		offMesh  = new CMesh(onMesh);
		tempMesh = new CMesh(onMesh);
		if (mRK4) {k1 = new CMesh(onMesh); k2 = new CMesh(onMesh);}
	} 
	else {
		onMesh  = new CMesh(pMap); 
		offMesh = new CMesh(onMesh);
		tempMesh = new CMesh(onMesh);
		if (mRK4) {k1 = new CMesh(onMesh); k2 = new CMesh(onMesh);}
	}
	
	mEos = new CEos(pMap);
    onMesh->setEos(mEos);  // static, so just once
	onMesh->setParameterMap(pMap);
	onMesh->deaden();

	if (mInitNS != 0.)   onMesh->initNS();
	if (mInitFlow != 0.) onMesh->addInitialFlow();

	sLoss=0.; sLoss2=0.; eLoss=0.;
	pLoss[0] = 0.; pLoss[1] = 0.; pLoss[2] = 0.;
	epLoss=0.;
	
	if (mIoIntegrals) {
		intS = onMesh->integralS();
		intS0 = intS;
		intE = onMesh->integralE();
		intE0 = intE;
		for (int i=0;i<3;i++) intP0[i] = onMesh->integralP(i+1);
	  
		fprintf(fIPX,"%0.9g %0.9g\n", mT0, intP0[0]+pLoss[0]);
		fprintf(fIPY,"%0.9g %0.9g\n", mT0, intP0[1]+pLoss[2]);
		fprintf(fIPN,"%0.9g %0.9g\n", mT0, intP0[2]+pLoss[2]);
		fprintf(fPXLoss,"%0.9g 0.\n",mT0);
		fprintf(fPYLoss,"%0.9g 0.\n",mT0);
		fprintf(fPNLoss,"%0.9g 0.\n",mT0);
	}
	if (mIoSpots) {
		double mTau = onMesh->getTau();
		fprintf(fIE,"%0.9g %0.9g\n", mTau, intE/intE0 - 1.);
		if (mBjorken) fprintf(fE0,"%0.9g %0.9g\n", mTau, pow(mTau,4./3.));
		else fprintf(fE0,"%0.9g 1.\n", mTau);
		fprintf(fTxx0,"%0.6g %0.6g\n", mTau, onMesh->getTxx(0,0,0));
		fprintf(fTzz0,"%0.6g %0.6g\n", mTau, onMesh->getTzz(0,0,0));
		fprintf(fT0,"%0.6g %0.6g\n", mTau, onMesh->getT(0,0,0)/0.36);
		fprintf(fIS,"%0.9g %0.9g\n", mTau, 0.);//intS/intS0);
		fprintf(fUL,"%0.6g %0.6g\n", mTau, onMesh->getS(onMesh->getNSize()-1,onMesh->getXSize()-1,0,1));
		fprintf(fEL,"%0.6g %0.6g\n", mTau, onMesh->getE(onMesh->getNSize()-1,onMesh->getXSize()-1,0));
		fprintf(fA1L,"%0.6g %0.6g\n",mTau, onMesh->getS(onMesh->getNSize()-1,onMesh->getXSize()-1,0,5));
		fprintf(fA2L,"%0.6g %0.6g\n",mTau, onMesh->getS(onMesh->getNSize()-1,onMesh->getXSize()-1,0,6));
		fprintf(fBL,"%0.6g %0.6g\n", mTau, onMesh->getS(onMesh->getNSize()-1,onMesh->getXSize()-1,0,10));
		fprintf(fEP0,"%0.6g %0.6g\n", mTau - mT0, onMesh->integralEp(0));
		fprintf(fEX0,"%0.6g %0.6g\n", mTau - mT0, onMesh->integralEx(0));
		fprintf(fVR0,"%0.6g %0.6g\n", mTau - mT0, onMesh->averageVr(0));
		fprintf(fFOSL,"%0.6g %0.6g\n", onMesh->getFOP(), mTau);
		fprintf(fFOSSST,"%0.6g %0.6g\n", mTau, onMesh->getFOSSST());
		fprintf(fTxxL,"%0.6g %0.6g\n", mTau, onMesh->getPixy(0,mPrintX,0,1,1));
		fprintf(fTxxNSL,"%0.6g %0.6g\n", mTau, onMesh->getPixyNS(0,mPrintX,0,1,1));
		fprintf(fTxxML,"%0.6g %0.6g\n", mTau, onMesh->getPixyMesh(0,mPrintX,0,1,1));
		fprintf(fTxxNSML,"%0.6g %0.6g\n", mTau, onMesh->getPixyNSMesh(0,mPrintX,0,1,1));
		fprintf(fDzz,"%0.6g %0.6g\n",mTau, mTau*onMesh->getDULocal(0,0,0,3,3)-1.);
	}	
	if (mIoOscarFull) openOscarFull();
	if (mIoOscarHyper) openOscarHyper();

    time(&now);
	if (!mPureBjorken) 
		printf("Expected Run Time < %0.3g sec.  \n", (mXSize*mYSize*mNSize)/1000. + (mXSize*mYSize*mNSize*mTStep)/2.5E4);
	else 
		printf("Expected Run Time < %0.3g sec.  \n", (mXSize*mYSize*mTStep)/2.5E4);
	fflush(stdout);
    printf("Time Step 0 [%0.6g](elapsed time %0.6g sec) .",onMesh->getTau(),difftime(now,start)); fflush(stdout);
	
    for (int t=1;t<=mTStep;t++) {

		if (mSawtooth) onMesh->forward(tempMesh,offMesh,mDt);
		else if (mRK2) {
			onMesh->forward(offMesh,mDt/2.);
			onMesh->forward(tempMesh,offMesh,mDt);
		}
		else if (mRK4) {
			onMesh->forward(offMesh,mDt/2.); 
			onMesh->forward(k1,offMesh,mDt/2.);
			onMesh->forward(k2,k1,mDt/2.);
			k2->setTau(k1->getTau());     //k[2] mesh is actually at half step
			onMesh->forward(tempMesh,k2,mDt);
			onMesh->forward(offMesh,k1,k2,tempMesh,mDt);
		}

		if (!tempMesh->detectCrash()) {
			printf("\n\n********** WARNING -- isHydro3 resulted in crash at t=%d.************\n Dumping pre-crash Mesh...... \n",t);
			printE(onMesh,5);
			return 1;
		}

		if (mIoIntegrals) printIntegrals(t);
	
		if (tempMesh->getT(0,0,0) < mFoTemp) {
			if (mIoSlices || mIoFull) printE(tempMesh,5);
			t=mTStep+1;
		}
	  
		if (t%mPrintStep==0 && (mIoSlices || mIoFull)) printE(tempMesh,t/mPrintStep);
		else if ((t*50)%mTStep==0) {printf(".");  fflush(stdout);}

		if (mIoSpots && t%mIoSliceTStep == 0) 
			printSpots();

		if (mIoTechqm && t%mIoTechqmTStep == 0)
			printTQM(tempMesh);
	
		if (mIoOscarFull && t%mIoOscarTStep == 0)
			printOscarFull(tempMesh,t);
	
		if (mIoOscarHyper && t%mIoOscarTStep == 0)
			printOscarHyper(tempMesh,t);
	  
		if (mHalving)
			if (tempMesh->getTau() - onMesh->getTau() > 2.*mDt*mT0) {
				printf("!");
				mDt /= 2.;
			}

		tempMesh->deaden();
		if (mSawtooth) {
			deadMesh = onMesh;
			onMesh = offMesh;
			offMesh = tempMesh;
			tempMesh = deadMesh;
		}
		else if (mRK2) {
			deadMesh = onMesh;
			onMesh = tempMesh;
			tempMesh = deadMesh;
		}
	  // no reassignments for RK4
	}

	time(&now);
	printf("\nHydro Finished...... (Total Time: %0.6g sec)\n",difftime(now,start)); fflush(stdout);
  
	onMesh->~CMesh(); 
	offMesh->~CMesh(); 
	if (mRK4) {
		k1->~CMesh();
		k2->~CMesh();
	}
	
	if (mIoSpectra) {
		printf("Computing Spectra.......\n"); fflush(stdout);
		printDNs(tempMesh,mNSize,mXSize,mYSize);
	}
	
	tempMesh->~CMesh();
  
  return 0;
}

void CHydro::printE(CMesh* lMesh, int lT) { 
	openFile(lT);
	time(&now);
	printf("\nAt Time: [%0.6g] (elapsed time %0.6g sec) .",lMesh->getTau(),difftime(now,start)); fflush(stdout);

	int nPrintX = min(parameter::getI(*pMap,"HYDRO_PRINTX",0),lMesh->getXSize()-1);
	int nPrintY = min(parameter::getI(*pMap,"HYDRO_PRINTY",0),lMesh->getYSize()-1);
	int lPrintN = min(parameter::getI(*pMap,"HYDRO_PRINTN",0),lMesh->getNSize()-1);
	if (mPureBjorken) lPrintN=0;

    if (!mPureBjorken)
		for (int i=-1;i<=lMesh->getNSize();i++) {
			double localN = lMesh->getX(i,nPrintX,nPrintY,3);
			if (parameter::getB(*pMap,"HYDRO_BJORKEN",true))  
				fprintf(fEN,"%0.9g %0.9g\n", localN, pow(lMesh->getTau(),(4./3.)) * lMesh->getS(i,nPrintX,nPrintY,4));
			else			
				fprintf(fEN,"%0.9g %0.9g\n", localN, lMesh->getS(i,nPrintX,nPrintY,4));

			fprintf(fEP, "%0.6g %0.6g\n", lMesh->getTau(), lMesh->integralEp(i));	
			fprintf(fEZ, "%0.9g %0.9g\n", localN, lMesh->getS(i,nPrintX,nPrintY,4));
			fprintf(fUz,"%0.9g %0.9g\n", localN, lMesh->getS(i,nPrintX,nPrintY,3));
			fprintf(fUzz,"%0.9g %0.9g\n", localN, lMesh->getS(i,nPrintX,nPrintY,3));

			fprintf(fA1z,"%0.9g %0.9g\n", localN, lMesh->getS(i,nPrintX,nPrintY,5));
			fprintf(fA2z,"%0.9g %0.9g\n", localN, lMesh->getS(i,nPrintX,nPrintY,6));
			fprintf(fA3z,"%0.9g %0.9g\n", localN, lMesh->getS(i,nPrintX,nPrintY,7));
			fprintf(fA4z,"%0.9g %0.9g\n", localN, lMesh->getS(i,nPrintX,nPrintY,8));
			fprintf(fA5z,"%0.9g %0.9g\n", localN, lMesh->getS(i,nPrintX,nPrintY,9));
			fprintf(fBz, "%0.9g %0.9g\n", localN, lMesh->getS(i,nPrintX,nPrintY,10));

			fprintf(fT,"%0.9g %0.9g\n", localN, lMesh->getT(i,nPrintX,nPrintY));
			fprintf(fP,"%0.9g %0.9g\n", localN, lMesh->getP(i,nPrintX,nPrintY));
			fprintf(fTxx,"%0.9g %0.9g\n", localN,lMesh->getTxyMesh(i,nPrintX,nPrintY,1,1));
			fprintf(fTyy,"%0.9g %0.9g\n", localN,lMesh->getTxyMesh(i,nPrintX,nPrintY,2,2));
			fprintf(fTzz,"%0.9g %0.9g\n", localN,lMesh->getTxyMesh(i,nPrintX,nPrintY,3,3));
			fprintf(fTxy,"%0.9g %0.9g\n", localN,lMesh->getTxyMesh(i,nPrintX,nPrintY,1,2));
			fprintf(fTyz,"%0.9g %0.9g\n", localN,lMesh->getTxyMesh(i,nPrintX,nPrintY,2,3));
		}

	double mFOS[XSIZE+YSIZE][2];
	double mFOSigma[XSIZE+YSIZE][2];
	int fosSize;
	lMesh->getFOS(mFOS,mFOSigma,fosSize);
	
	for (int i=0;i<fosSize;i++)
	    fprintf(fFOS,"%0.6g %0.6g\n", mFOS[i][0], mFOS[i][1]);

	int lPrintY = nPrintY;
	for (int i=0;i<lMesh->getXSize();i++) {
		double localX = lMesh->getX(lPrintN,i,lPrintY,1);
		fprintf(fEX,"%0.9g %0.9g\n", localX, lMesh->getS(lPrintN,i,lPrintY,4));
		fprintf(fS,"%0.9g %0.9g\n", localX, lMesh->getS(lPrintN,i,lPrintY));

		fprintf(fU0,"%0.9g %0.9g\n", localX, lMesh->getS(lPrintN,i,lPrintY,0));
		fprintf(fUx,"%0.9g %0.9g\n", localX, lMesh->getS(lPrintN,i,lPrintY,1));

		fprintf(fA1x,"%0.9g %0.9g\n", localX, lMesh->getS(lPrintN,i,lPrintY,5));
		fprintf(fA2x,"%0.9g %0.9g\n", localX, lMesh->getS(lPrintN,i,lPrintY,6));
		fprintf(fA3x,"%0.9g %0.9g\n", localX, lMesh->getS(lPrintN,i,lPrintY,7));
		fprintf(fA4x,"%0.9g %0.9g\n", localX, lMesh->getS(lPrintN,i,lPrintY,8));
		fprintf(fA5x,"%0.9g %0.9g\n", localX, lMesh->getS(lPrintN,i,lPrintY,9));
		fprintf(fBx, "%0.9g %0.9g\n", localX, lMesh->getS(lPrintN,i,lPrintY,10));

		double P = lMesh->getP(lPrintN,i,lPrintY);
		fprintf(fT,"%0.9g %0.9g\n", localX, lMesh->getT(lPrintN,i,lPrintY));
		fprintf(fP,"%0.9g %0.9g\n", localX, P);

		fprintf(fTxx,"%0.9g %0.9g\n", localX, lMesh->getPixyMesh(lPrintN,i,lPrintY,1,1)/P);
		fprintf(fTrr,"%0.9g %0.9g\n", localX, lMesh->getTRR(lPrintN,i,lPrintY)/P);
		fprintf(fTpp,"%0.9g %0.9g\n", localX, lMesh->getTPhiPhi(lPrintN,i,lPrintY)/P);
		fprintf(fTpr,"%0.9g %0.9g\n", localX, lMesh->getTRPhi(lPrintN,i,lPrintY)/P);

		fprintf(fTxxNS,"%0.9g %0.9g\n",localX, lMesh->getPixyNS(lPrintN,i,lPrintY,1,1)/P);
		fprintf(fTyy,"%0.9g %0.9g\n",  localX, lMesh->getPixyMesh(lPrintN,i,lPrintY,2,2)/P);
		fprintf(fTzz,"%0.9g %0.9g\n",  localX, lMesh->getPixyMesh(lPrintN,i,lPrintY,3,3)/P);

		fprintf(fTIS, "%0.9g %0.9g\n", localX, lMesh->getTIS(lPrintN,i,lPrintY));
		fprintf(fTISB,"%0.9g %0.9g\n", localX, lMesh->getTISB(lPrintN,i,lPrintY));

		if (mPureBjorken) 
			fprintf(fUz,"%0.9g %0.9g\n", lMesh->getX(lPrintN,i,lPrintY,2), lMesh->getS(lPrintN,i,lPrintY,3));
	}
	
	int mPrintN = lPrintN, mPrintX = nPrintX;
	for (int i=0; i<lMesh->getYSize();i++) {
		double localY = lMesh->getX(mPrintN,mPrintX,i,2);
		fprintf(fEY, "%0.9g %0.9g\n", localY, lMesh->getS(mPrintN,mPrintX,i,4));
		
		fprintf(fUy, "%0.9g %0.9g\n", localY,lMesh->getS(mPrintN,mPrintX,i,2));
	  
		fprintf(fA1y,"%0.9g %0.9g\n", localY, lMesh->getS(mPrintN,mPrintX,i,5));
		fprintf(fA2y,"%0.9g %0.9g\n", localY, lMesh->getS(mPrintN,mPrintX,i,6));
		fprintf(fA3y,"%0.9g %0.9g\n", localY, lMesh->getS(mPrintN,mPrintX,i,7));
		fprintf(fA4y,"%0.9g %0.9g\n", localY, lMesh->getS(mPrintN,mPrintX,i,8));
		fprintf(fA5y,"%0.9g %0.9g\n", localY, lMesh->getS(mPrintN,mPrintX,i,9));
		fprintf(fBy, "%0.9g %0.9g\n", localY, lMesh->getS(mPrintN,mPrintX,i,10));
	}

	for (int i=0;i < min(lMesh->getYSize(),lMesh->getXSize());i++) {
		double P = lMesh->getP(lPrintN,i,i);
		double localR = sqrt( lMesh->getX(lPrintN,i,i,1)*lMesh->getX(lPrintN,i,i,1) + lMesh->getX(lPrintN,i,i,2)*lMesh->getX(lPrintN,i,i,2));
		
		fprintf(fTxy,"%0.9g %0.9g\n", localR, lMesh->getPixy(lPrintN,i,i,1,2)/P);
		fprintf(fTyz,"%0.9g %0.9g\n", localR, lMesh->getPixy(lPrintN,i,i,2,3)/P);
		fprintf(fU0 ,"%0.9g %0.9g\n", localR, lMesh->getS(lPrintN,i,i,0));
	  
		fprintf(fTrr,"%0.9g %0.9g\n", localR, lMesh->getTRR(lPrintN,i,i)/P);
		fprintf(fTpp,"%0.9g %0.9g\n", localR, lMesh->getTPhiPhi(lPrintN,i,i)/P);
		fprintf(fTpr,"%0.9g %0.9g\n", localR, lMesh->getTRPhi(lPrintN,i,i)/P);
	}

    if (!mOctant)
		for (int i=-lMesh->getNSize();i<=lMesh->getNSize();i++) {
			fprintf(fIEN,"%0.9g %0.9g\n", lMesh->getYL(i,nPrintX,nPrintY), lMesh->integralE(i));
			fprintf(fISN,"%0.9g %0.9g\n", lMesh->getYL(i,nPrintX,nPrintY), lMesh->integralS(i));  
		}
	else 
		for (int i=0;i<=lMesh->getNSize();i++) {
			fprintf(fIEN,"%0.9g %0.9g\n", lMesh->getYL(i,nPrintX,nPrintY), lMesh->integralE(i));
			fprintf(fISN,"%0.9g %0.9g\n", lMesh->getYL(i,nPrintX,nPrintY), lMesh->integralS(i));  
		}

	if (parameter::getB(*pMap,"HYDRO_IO_FULL",false)) printMesh(lMesh);
	closeFile();
}

void CHydro::openFile(int lT) {
	string temp;
	if      (lT == 1) {
		temp = mDataRoot;
		temp.append("fEx1.txt");
		fEX = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fEy1.txt");
		fEY = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fEz1.txt");
		fEZ = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fS1.txt");
		fS = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fT1.txt");
		fT = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fP1.txt");
		fP = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fEN1.txt");
		fEN = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fIEN1.txt");
		fIEN = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fISN1.txt");
		fISN = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fU1.txt");
		fU0 = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fUx1.txt");
		fUx = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fUy1.txt");
		fUy = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fUz1.txt");
		fUz = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fUzz1.txt");
		fUzz = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fULz1.txt");
		fULz = fopen(temp.c_str(),"w");
			  
	  
		temp = mDataRoot;
		temp.append("fA1x-1.txt");
		fA1x = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fA2x-1.txt");
		fA2x = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fA3x-1.txt");
		fA3x = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fA4x-1.txt");
		fA4x = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fA5x-1.txt");
		fA5x = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fB1x.txt");
		fBx = fopen(temp.c_str(),"w");
	  
	  
		temp = mDataRoot;
		temp.append("fA1y-1.txt");
		fA1y = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fA2y-1.txt");
		fA2y = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fA3y-1.txt");
		fA3y = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fA4y-1.txt");
		fA4y = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fA5y-1.txt");
		fA5y = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fB1y.txt");
		fBy = fopen(temp.c_str(),"w");
	  
	  
		temp = mDataRoot;
		temp.append("fA1z-1.txt");
		fA1z = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fA2z-1.txt");
		fA2z = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fA3z-1.txt");
		fA3z = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fA4z-1.txt");
		fA4z = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fA5z-1.txt");
		fA5z = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fB1z.txt");
		fBz = fopen(temp.c_str(),"w");
	  
	  
		temp = mDataRoot;
		temp.append("fTxx1.txt");
		fTxx = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fTxxNS1.txt");
		fTxxNS = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fTyy1.txt");
		fTyy = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fTzz1.txt");
		fTzz = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fTxy1.txt");
		fTxy = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fTyz1.txt");
		fTyz = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fTrr1.txt");
		fTrr = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fTpp1.txt");
		fTpp = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fTpr1.txt");
		fTpr = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fTIS1.txt");
		fTIS = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fTISB1.txt");
		fTISB = fopen(temp.c_str(),"w");

		temp = mDataRoot;
		temp.append("fFull1.txt");
		fFull = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fEP1.txt");
		fEP = fopen(temp.c_str(),"w");
	  
		temp = mDataRoot;
		temp.append("fFOS1.txt");
		fFOS = fopen(temp.c_str(),"w");
	}  	
	else if (lT == 2) { 
	  temp = mDataRoot;
	  temp.append("fEx2.txt");
	  fEX = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fEy2.txt");
	  fEY = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fEz2.txt");
	  fEZ = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fS2.txt");
  	  fS = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fT2.txt");
	  fT = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fP2.txt");
	  fP = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fEN2.txt");
	  fEN = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fIEN2.txt");
	  fIEN = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fISN2.txt");
	  fISN = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fU2.txt");
	  fU0 = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fUx2.txt");
	  fUx = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fUy2.txt");
	  fUy = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fUz2.txt");
	  fUz = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fUzz2.txt");
	  fUzz = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fULz2.txt");
	  fULz = fopen(temp.c_str(),"w");
	  
	  
	  temp = mDataRoot;
	  temp.append("fA1x-2.txt");
	  fA1x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA2x-2.txt");
	  fA2x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA3x-2.txt");
	  fA3x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA4x-2.txt");
      fA4x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA5x-2.txt");
	  fA5x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fB2x.txt");
	  fBx = fopen(temp.c_str(),"w");
	  
	  
	  temp = mDataRoot;
	  temp.append("fA1y-2.txt");
	  fA1y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA2y-2.txt");
	  fA2y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA3y-2.txt");
	  fA3y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA4y-2.txt");
      fA4y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA5y-2.txt");
	  fA5y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fB2y.txt");
	  fBy = fopen(temp.c_str(),"w");
	  
	  
	  temp = mDataRoot;
	  temp.append("fA1z-2.txt");
	  fA1z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA2z-2.txt");
	  fA2z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA3z-2.txt");
	  fA3z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA4z-2.txt");
      fA4z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA5z-2.txt");
	  fA5z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fB2z.txt");
	  fBz = fopen(temp.c_str(),"w");
	  
	  
	  temp = mDataRoot;
	  temp.append("fTxx2.txt");
	  fTxx = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTxxNS2.txt");
	  fTxxNS = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTyy2.txt");
	  fTyy = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTzz2.txt");
	  fTzz = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTxy2.txt");
	  fTxy = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTyz2.txt");
	  fTyz = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTrr2.txt");
	  fTrr = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTpp2.txt");
	  fTpp = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTpr2.txt");
	  fTpr = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTIS2.txt");
	  fTIS = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTISB2.txt");
	  fTISB = fopen(temp.c_str(),"w");

	  temp = mDataRoot;
	  temp.append("fFull2.txt");
	  fFull = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fEP2.txt");
	  fEP = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fFOS2.txt");
	  fFOS = fopen(temp.c_str(),"w");
	}
	else if (lT == 3) {
	  temp = mDataRoot;
	  temp.append("fEx3.txt");
	  fEX = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fEy3.txt");
	  fEY = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fEz3.txt");
	  fEZ = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fS3.txt");
  	  fS = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fT3.txt");
	  fT = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fP3.txt");
	  fP = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fEN3.txt");
	  fEN = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fIEN3.txt");
	  fIEN = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fISN3.txt");
	  fISN = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fU3.txt");
	  fU0 = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fUx3.txt");
	  fUx = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fUy3.txt");
	  fUy = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fUz3.txt");
	  fUz = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fUzz3.txt");
	  fUzz = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fULz3.txt");
	  fULz = fopen(temp.c_str(),"w");
	  
	  
	  temp = mDataRoot;
	  temp.append("fA1x-3.txt");
	  fA1x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA2x-3.txt");
	  fA2x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA3x-3.txt");
	  fA3x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA4x-3.txt");
      fA4x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA5x-3.txt");
	  fA5x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fB3x.txt");
	  fBx = fopen(temp.c_str(),"w");
	  
	  
	  temp = mDataRoot;
	  temp.append("fA1y-3.txt");
	  fA1y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA2y-3.txt");
	  fA2y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA3y-3.txt");
	  fA3y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA4y-3.txt");
      fA4y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA5y-3.txt");
	  fA5y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fB3y.txt");
	  fBy = fopen(temp.c_str(),"w");
	  
	  
	  temp = mDataRoot;
	  temp.append("fA1z-3.txt");
	  fA1z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA2z-3.txt");
	  fA2z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA3z-3.txt");
	  fA3z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA4z-3.txt");
      fA4z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA5z-3.txt");
	  fA5z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fB3z.txt");
	  fBz = fopen(temp.c_str(),"w");
	  
	  
	  temp = mDataRoot;
	  temp.append("fTxx3.txt");
	  fTxx = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTxxNS3.txt");
	  fTxxNS = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTyy3.txt");
	  fTyy = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTzz3.txt");
	  fTzz = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTxy3.txt");
	  fTxy = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTyz3.txt");
	  fTyz = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTrr3.txt");
	  fTrr = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTpp3.txt");
	  fTpp = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTpr3.txt");
	  fTpr = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTIS3.txt");
	  fTIS = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTISB3.txt");
	  fTISB = fopen(temp.c_str(),"w");

	  temp = mDataRoot;
	  temp.append("fFull3.txt");
	  fFull = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fEP3.txt");
	  fEP = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fFOS3.txt");
	  fFOS = fopen(temp.c_str(),"w");
	}
	else if (lT == 4) {
	  temp = mDataRoot;
	  temp.append("fEx4.txt");
	  fEX = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fEy4.txt");
	  fEY = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fEz4.txt");
	  fEZ = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fS4.txt");
  	  fS = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fT4.txt");
	  fT = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fP4.txt");
	  fP = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fEN4.txt");
	  fEN = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fIEN4.txt");
	  fIEN = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fISN4.txt");
	  fISN = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fU4.txt");
	  fU0 = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fUx4.txt");
	  fUx = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fUy4.txt");
	  fUy = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fUz4.txt");
	  fUz = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fUzz4.txt");
	  fUzz = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fULz4.txt");
	  fULz = fopen(temp.c_str(),"w");
	  
	  
	  temp = mDataRoot;
	  temp.append("fA1x-4.txt");
	  fA1x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA2x-4.txt");
	  fA2x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA3x-4.txt");
	  fA3x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA4x-4.txt");
      fA4x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA5x-4.txt");
	  fA5x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fB4x.txt");
	  fBx = fopen(temp.c_str(),"w");
	  
	  
	  temp = mDataRoot;
	  temp.append("fA1y-4.txt");
	  fA1y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA2y-4.txt");
	  fA2y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA3y-4.txt");
	  fA3y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA4y-4.txt");
      fA4y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA5y-4.txt");
	  fA5y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fB4y.txt");
	  fBy = fopen(temp.c_str(),"w");
	  
	  
	  temp = mDataRoot;
	  temp.append("fA1z-4.txt");
	  fA1z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA2z-4.txt");
	  fA2z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA3z-4.txt");
	  fA3z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA4z-4.txt");
      fA4z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA5z-4.txt");
	  fA5z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fB4z.txt");
	  fBz = fopen(temp.c_str(),"w");
	  
	  
	  temp = mDataRoot;
	  temp.append("fTxx4.txt");
	  fTxx = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTxxNS4.txt");
	  fTxxNS = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTyy4.txt");
	  fTyy = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTzz4.txt");
	  fTzz = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTxy4.txt");
	  fTxy = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTyz4.txt");
	  fTyz = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTrr4.txt");
	  fTrr = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTpp4.txt");
	  fTpp = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTpr4.txt");
	  fTpr = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTIS4.txt");
	  fTIS = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTISB4.txt");
	  fTISB = fopen(temp.c_str(),"w");

	  temp = mDataRoot;
	  temp.append("fFull4.txt");
	  fFull = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fEP4.txt");
	  fEP = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fFOS4.txt");
	  fFOS = fopen(temp.c_str(),"w");
	}
	else if (lT == 5) {
	temp = mDataRoot;
	  temp.append("fEx5.txt");
	  fEX = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fEy5.txt");
	  fEY = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fEz5.txt");
	  fEZ = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fS5.txt");
  	  fS = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fT5.txt");
	  fT = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fP5.txt");
	  fP = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fEN5.txt");
	  fEN = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fIEN5.txt");
	  fIEN = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fISN5.txt");
	  fISN = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fU5.txt");
	  fU0 = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fUx5.txt");
	  fUx = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fUy5.txt");
	  fUy = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fUz5.txt");
	  fUz = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fUzz5.txt");
	  fUzz = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fULz5.txt");
	  fULz = fopen(temp.c_str(),"w");
	  
	  
	  temp = mDataRoot;
	  temp.append("fA1x-5.txt");
	  fA1x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA2x-5.txt");
	  fA2x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA3x-5.txt");
	  fA3x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA4x-5.txt");
      fA4x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA5x-5.txt");
	  fA5x = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fB5x.txt");
	  fBx = fopen(temp.c_str(),"w");
	  
	  
	  temp = mDataRoot;
	  temp.append("fA1y-5.txt");
	  fA1y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA2y-5.txt");
	  fA2y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA3y-5.txt");
	  fA3y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA4y-5.txt");
      fA4y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA5y-5.txt");
	  fA5y = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fB5y.txt");
	  fBy = fopen(temp.c_str(),"w");
	  
	  
	  temp = mDataRoot;
	  temp.append("fA1z-5.txt");
	  fA1z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA2z-5.txt");
	  fA2z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA3z-5.txt");
	  fA3z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA4z-5.txt");
      fA4z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fA5z-5.txt");
	  fA5z = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fB5z.txt");
	  fBz = fopen(temp.c_str(),"w");
	  
	  
	  temp = mDataRoot;
	  temp.append("fTxx5.txt");
	  fTxx = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTxxNS5.txt");
	  fTxxNS = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTyy5.txt");
	  fTyy = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTzz5.txt");
	  fTzz = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTxy5.txt");
	  fTxy = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTyz5.txt");
	  fTyz = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTrr5.txt");
	  fTrr = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTpp5.txt");
	  fTpp = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTpr5.txt");
	  fTpr = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTIS5.txt");
	  fTIS = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fTISB5.txt");
	  fTISB = fopen(temp.c_str(),"w");

	  temp = mDataRoot;
	  temp.append("fFull5.txt");
	  fFull = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fEP5.txt");
	  fEP = fopen(temp.c_str(),"w");
	  
	  temp = mDataRoot;
	  temp.append("fFOS5.txt");
	  fFOS = fopen(temp.c_str(),"w");
	}
}

void CHydro::openFileIntegrals() {
    string temp = mDataRoot;
	temp.append("fIE.txt");
    fIE = fopen( temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fIPX.txt");
	fIPX = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fIPY.txt");
	fIPY = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fIPN.txt");
	fIPN = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fPXLoss.txt");
	fPXLoss = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fPYLoss.txt");
	fPYLoss = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fPNLoss.txt");
	fPNLoss = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fIS.txt");
	fIS = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fIS2.txt");
    fIS2 = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fSLoss.txt");
    fSLoss = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fSLoss2.txt");
    fSLoss2 = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fELoss.txt");
    fELoss = fopen(temp.c_str(),"w");
}

void CHydro::printIntegrals(int t) {
	intS = tempMesh->integralS();  
	intE = tempMesh->integralE();
	for (int i=0;i<3;i++) intP[i] = tempMesh->integralP(i+1);
	if (mRK2) 
		if (mLogSinhT) {
			eLoss += mDt * tanh(offMesh->getTau()) * offMesh->getELoss(tempMesh);
			sLoss += mDt * tanh(offMesh->getTau()) * offMesh->getSLoss();
			for (int i=0;i<3;i++) 
				pLoss[i] += mDt * tanh(offMesh->getTau()) * offMesh->getPLoss(tempMesh,i+1);
		} 
		else if (mLogT) {
			eLoss += mDt * offMesh->getTau() * offMesh->getELoss(tempMesh);
			sLoss += mDt * offMesh->getTau() * offMesh->getSLoss();
			for (int i=0;i<3;i++) 
				pLoss[i] += mDt * offMesh->getTau() * offMesh->getPLoss(tempMesh,i+1);
		}
		else {
			eLoss += (tempMesh->getTau() - onMesh->getTau()) * onMesh->getELoss(tempMesh);
			sLoss += (tempMesh->getTau() - onMesh->getTau()) * onMesh->getSLoss();
			for (int i=0;i<3;i++) 
				pLoss[i] += mDt * onMesh->getPLoss(tempMesh,i+1);
		}

	if (t%mIoSliceTStep == 0) {
		double mTau = tempMesh->getTau();
		fprintf(fIE,"%0.9g %0.9g\n", mTau, (intE+eLoss)/intE0 - 1.);
		if (intP0[0] != 0.) fprintf(fIPX,"%0.9g %0.9g\n", mTau, (intP[0]+pLoss[0])/intP0[0] - 1.);
		else				fprintf(fIPX,"%0.9g %0.9g\n", mTau,  intP[0]+pLoss[0]);
		if(intP0[1] != 0.)  fprintf(fIPY,"%0.9g %0.9g\n", mTau, (intP[1]+pLoss[1])/intP0[1] - 1.);
		else				fprintf(fIPY,"%0.9g %0.9g\n", mTau,  intP[1]+pLoss[1]);
		if (intP0[2] != 0.) fprintf(fIPN,"%0.9g %0.9g\n", mTau, (intP[2]+pLoss[2])/intP0[2] - 1.);
		else				fprintf(fIPN,"%0.9g %0.9g\n", mTau,  intP[2]+pLoss[2]);
		fprintf(fIS,"%0.9g %0.9g\n", mTau, (intS+sLoss)/intS0 - 1.);
		fprintf(fSLoss,"%0.6g %0.6g\n", mTau, sLoss/intS0);
		fprintf(fELoss,"%0.6g %0.6g\n", mTau, eLoss/intE0);
		if (intP0[0] != 0.) fprintf(fPXLoss,"%0.9g %0.9g\n", mTau, pLoss[0]/intP0[0]);
		else				fprintf(fPXLoss,"%0.9g %0.9g\n", mTau, pLoss[0]);
		if (intP0[1] != 0.) fprintf(fPYLoss,"%0.9g %0.9g\n", mTau, pLoss[1]/intP0[1]);
		else				fprintf(fPYLoss,"%0.9g %0.9g\n", mTau, pLoss[1]);
		if (intP0[2] != 0.) fprintf(fPNLoss,"%0.9g %0.9g\n", mTau, pLoss[2]/intP0[2]);
		else				fprintf(fPNLoss,"%0.9g %0.9g\n", mTau, pLoss[2]);
	}
}

void CHydro::printSpots() {

	double mTau = tempMesh->getTau();
	
	if (mBjorken)
		fprintf(fE0,"%0.9g %0.9g\n", mTau, tempMesh->getE(0,0,0)/mE0 *pow(tempMesh->getTau(),4./3.));
	else fprintf(fE0,"%0.9g %0.9g\n", mTau, tempMesh->getE(0,0,0)/mE0);
		fprintf(fT0,"%0.6g %0.6g\n", mTau, tempMesh->getT(0,0,0)/0.36);
	
	fprintf(fTxx0,"%0.6g %0.6g\n", mTau, tempMesh->getTxx(0,0,0));
	fprintf(fTzz0,"%0.6g %0.6g\n", mTau, tempMesh->getTzz(0,0,0));
	if (mPureBjorken) {
		fprintf(fUL,"%0.6g %0.6g\n", mTau,  tempMesh->getS(0,tempMesh->getXSize()-1,0,1));
		fprintf(fEL,"%0.6g %0.6g\n", mTau,  tempMesh->getE(0,tempMesh->getXSize()-1,0));
		fprintf(fA1L,"%0.6g %0.6g\n", mTau, tempMesh->getS(0,tempMesh->getXSize()-1,0,5));
		fprintf(fA2L,"%0.6g %0.6g\n", mTau, tempMesh->getS(0,tempMesh->getXSize()-1,0,6));
		fprintf(fBL,"%0.6g %0.6g\n", mTau,  tempMesh->getS(0,tempMesh->getXSize()-1,0,10));
	} 
	else {
		fprintf(fUL,"%0.6g %0.6g\n", mTau,  tempMesh->getS(tempMesh->getNSize()-1,0,0,3));
		fprintf(fEL,"%0.6g %0.6g\n", mTau,  tempMesh->getE(0,0,0)/mE0);
		fprintf(fA1L,"%0.6g %0.6g\n", mTau, tempMesh->getS(tempMesh->getNSize()-1,0,0,5));
		fprintf(fA2L,"%0.6g %0.6g\n", mTau, tempMesh->getS(tempMesh->getNSize()-1,0,0,6));
		fprintf(fBL,"%0.6g %0.6g\n", mTau,  tempMesh->getS(tempMesh->getNSize()-1,0,0,10));
	}
	
	fprintf(fEP0,"%0.6g %0.6g\n", mTau - mT0, tempMesh->integralEp(0));
	fprintf(fEX0,"%0.6g %0.6g\n", mTau - mT0, tempMesh->integralEx(0));
	fprintf(fVR0,"%0.6g %0.6g\n", mTau - mT0, tempMesh->averageVr(0));
	fprintf(fFOSL,"%0.6g %0.6g\n", tempMesh->getFOP(), mTau);
	fprintf(fFOSSST,"%0.6g %0.6g\n", mTau, tempMesh->getFOSSST());
	fprintf(fTxxL,"%0.6g %0.6g\n", mTau, tempMesh->getPixy(0,mPrintX,0,1,1));
	fprintf(fTxxNSL,"%0.6g %0.6g\n", mTau, tempMesh->getPixyNS(0,mPrintX,0,1,1));
	fprintf(fTxxML,"%0.6g %0.6g\n", mTau, tempMesh->getPixyMesh(0,0,mPrintX,1,1));
	fprintf(fTxxNSML,"%0.6g %0.6g\n", mTau, tempMesh->getPixyNSMesh(0,mPrintX,0,1,1));
	fprintf(fDzz,"%0.6g %0.6g\n",mTau, mTau*tempMesh->getDULocal(0,0,0,3,3)-1.);
}

void CHydro::openFileSpots() {
	string temp = mDataRoot;
	temp.append("fE0.txt");
    fE0 = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fT0.txt");
    fT0  = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fTxx0.txt");
    fTxx0 = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fTzz0.txt");
    fTzz0 = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fUL.txt");
	fUL = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fEL.txt");
	fEL = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fA1L.txt");
	fA1L = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fA2L.txt");
	fA2L = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fBL.txt");
	fBL  = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fDzz.txt");
	fDzz = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fEp0.txt");
	fEP0 = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fEx0.txt");
	fEX0 = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fVr0.txt");
	fVR0 = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fFOSL.txt");
	fFOSL = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fFOSSST.txt");
	fFOSSST = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fTxxL.txt");
	fTxxL = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fTxxNSL.txt");
	fTxxNSL = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fTxxML.txt");
	fTxxML = fopen(temp.c_str(),"w");
	
	temp = mDataRoot;
	temp.append("fTxxNSML.txt");
	fTxxNSML = fopen(temp.c_str(),"w");
}

void CHydro::closeFile() {
  fclose(fEX); fclose(fEY); fclose(fEZ); fclose(fS); fclose(fT); fclose(fP);
  fclose(fEN); fclose(fU0); fclose(fUx); fclose(fUy); fclose(fUz); fclose(fULz); fclose(fUzz);
  fclose(fE0); fclose(fTxx0); fclose(fTzz0);
  fclose(fFull);
  fclose(fA1x); fclose(fA2x); fclose(fA3x); fclose(fA4x); fclose(fA5x); fclose(fBx);
  fclose(fA1y); fclose(fA2y); fclose(fA3y); fclose(fA4y); fclose(fA5y); fclose(fBy);
  fclose(fA1z); fclose(fA2z); fclose(fA3z); fclose(fA4z); fclose(fA5z); fclose(fBz);
  fclose(fTxx); fclose(fTyy); fclose(fTzz); fclose(fTxy); fclose(fTxz); fclose(fTyz); fclose(fTxxNS);
  fclose(fTrr); fclose(fTpp); fclose(fTpr);
  fclose(fTIS); fclose(fTISB);
  fclose(fFOS);
}

void CHydro::printMesh(CMesh* lMesh) {
	int mN = lMesh->getNSize();
	int mX = lMesh->getXSize();
	int mY = lMesh->getYSize();	
	fprintf(fFull,"%0.6g %d %d %d \n",lMesh->getTau(),mN,mX,mY);
	if (!mPureBjorken)
	  for (int i=-1;i<=mN;i++) 
	    for (int j=-1;j<=mX;j++) 
	      for (int k=-1;k<=mY;k++) {
		    for (int l=1;l<4;l++)  fprintf(fFull,"%0.6g ",lMesh->getX(i,j,k,l));
		    for (int l=0;l<11;l++) fprintf(fFull,"%0.9g ",lMesh->getS(i,j,k,l));
		    fprintf(fFull,"\n");
		  }
	else 
	  for (int j=-1;j<=mX;j++) 
	      for (int k=-1;k<=mY;k++) {
		    for (int l=1;l<4;l++)  fprintf(fFull,"%0.6g ",lMesh->getX(0,j,k,l));
		    for (int l=0;l<11;l++) fprintf(fFull,"%0.9g ",lMesh->getS(0,j,k,l));
		    fprintf(fFull,"\n");
		  }
}

void CHydro::printTQM(CMesh* lMesh) {
	double mZ = 0.;
	int mX = lMesh->getXSize();
	int mY = lMesh->getYSize();	
	for (int j=0;j<mX;j++) {
	  int k=0;
	  lMesh->update(0,j,k);
	  fprintf(fTQMX,"%0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g\n",
		  lMesh->getTau(), lMesh->getX(0,j,k,1), lMesh->getX(0,j,k,2), sqrt(lMesh->getX(0,j,k,1)*lMesh->getX(0,j,k,1)+lMesh->getX(0,j,k,2)*lMesh->getX(0,j,k,2)), //1-4
		  lMesh->getS(0,j,k,4), lMesh->getS(0,j,k,1)/lMesh->getS(0,j,k,0), lMesh->getS(0,j,k,2)/lMesh->getS(0,j,k,0),  //5-7
		  sqrt(lMesh->getS(0,j,k,1)*lMesh->getS(0,j,k,1) + lMesh->getS(0,j,k,2)*lMesh->getS(0,j,k,2))/lMesh->getS(0,j,k,0), //8
		  lMesh->getS(0,j,k,0), lMesh->getDUDT(0,j,k,0),lMesh->getTxyMesh(0,j,k,0,0), lMesh->getTxyMesh(0,j,k,0,1), //9-12
		  lMesh->getTxyMesh(0,j,k,0,2), lMesh->getP(0,j,k), lMesh->getBV(0,j,k), //13-15
		  (lMesh->getS(0,j,k,1)*lMesh->getDUDT(0,j,k,0) + lMesh->getS(0,j,k,2)*lMesh->getDUDT(0,j,k,1) + lMesh->getS(0,j,k,3)*lMesh->getDUDT(0,j,k,2))/lMesh->getS(0,j,k,0) 
		  + lMesh->getDS(0,j,k,0,0) + lMesh->getDS(0,j,k,1,1) + lMesh->getDS(0,j,k,2,2), //16
		  lMesh->getPixyMesh(0,j,k,0,0), lMesh->getPixyMesh(0,j,k,0,1), lMesh->getPixyMesh(0,j,k,0,2), lMesh->getPixyMesh(0,j,k,3,3), //17-20
		  lMesh->getPixyMesh(0,j,k,1,1), lMesh->getPixyMesh(0,j,k,2,2), lMesh->getPixyMesh(0,j,k,1,2), mZ, mZ, mZ, //21-26
		  lMesh->getPiZR(0,j,k), lMesh->getPiZPhi(0,j,k), lMesh->getPiRR(0,j,k), lMesh->getPiPhiPhi(0,j,k), lMesh->getPiRPhi(0,j,k), //27-31
		  lMesh->getT(0,j,k));
	}
	fprintf(fTQMX,"\n");
	
	for (int k=0;k<mY;k++) {
	  int j=0;
	  lMesh->update(0,j,k);
	  fprintf(fTQMY,"%0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g\n",
		  lMesh->getTau(), lMesh->getX(0,j,k,1), lMesh->getX(0,j,k,2), sqrt(lMesh->getX(0,j,k,1)*lMesh->getX(0,j,k,1)+lMesh->getX(0,j,k,2)*lMesh->getX(0,j,k,2)), //1-4
		  lMesh->getS(0,j,k,4), lMesh->getS(0,j,k,1)/lMesh->getS(0,j,k,0), lMesh->getS(0,j,k,2)/lMesh->getS(0,j,k,0),  //5-7
		  sqrt(lMesh->getS(0,j,k,1)*lMesh->getS(0,j,k,1) + lMesh->getS(0,j,k,2)*lMesh->getS(0,j,k,2))/lMesh->getS(0,j,k,0), //8
		  lMesh->getS(0,j,k,0), lMesh->getDUDT(0,j,k,0),lMesh->getTxyMesh(0,j,k,0,0), lMesh->getTxyMesh(0,j,k,0,1), //9-12
		  lMesh->getTxyMesh(0,j,k,0,2), lMesh->getP(0,j,k), lMesh->getBV(0,j,k), lMesh->getDUDT(0,j,k,0) + lMesh->getDS(0,j,k,0,0) + lMesh->getDS(0,j,k,1,1) + lMesh->getDS(0,j,k,2,2), //13-16
		  lMesh->getPixyMesh(0,j,k,0,0), lMesh->getPixyMesh(0,j,k,0,1), lMesh->getPixyMesh(0,j,k,0,2), lMesh->getPixyMesh(0,j,k,2,2), //17-20
		  lMesh->getPixyMesh(0,j,k,1,1), lMesh->getPixyMesh(0,j,k,2,2), lMesh->getPixyMesh(0,j,k,1,2), mZ, mZ, mZ, //21-26
		  lMesh->getPiZR(0,j,k), lMesh->getPiZPhi(0,j,k), lMesh->getPiRR(0,j,k), lMesh->getPiPhiPhi(0,j,k), lMesh->getPiRPhi(0,j,k), //27-31
		  lMesh->getT(0,j,k));
	}
	fprintf(fTQMY,"\n");
	
	for (int j=0;j<min(mX,mY);j++) {
	  int k=j;
	  lMesh->update(0,j,k);
	  fprintf(fTQMR,"%0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g %0.6g\n",
		  lMesh->getTau(), lMesh->getX(0,j,k,1), lMesh->getX(0,j,k,2), sqrt(lMesh->getX(0,j,k,1)*lMesh->getX(0,j,k,1)+lMesh->getX(0,j,k,2)*lMesh->getX(0,j,k,2)), //1-4
		  lMesh->getS(0,j,k,4), lMesh->getS(0,j,k,1)/lMesh->getS(0,j,k,0), lMesh->getS(0,j,k,2)/lMesh->getS(0,j,k,0),  //5-7
		  sqrt(lMesh->getS(0,j,k,1)*lMesh->getS(0,j,k,1) + lMesh->getS(0,j,k,2)*lMesh->getS(0,j,k,2))/lMesh->getS(0,j,k,0), //8
		  lMesh->getS(0,j,k,0), lMesh->getDUDT(0,j,k,0),lMesh->getTxyMesh(0,j,k,0,0), lMesh->getTxyMesh(0,j,k,0,1), //9-12
		  lMesh->getTxyMesh(0,j,k,0,2), lMesh->getP(0,j,k), lMesh->getBV(0,j,k), lMesh->getDUDT(0,j,k,0) + lMesh->getDS(0,j,k,0,0) + lMesh->getDS(0,j,k,1,1) + lMesh->getDS(0,j,k,2,2), //13-16
		  lMesh->getPixyMesh(0,j,k,0,0), lMesh->getPixyMesh(0,j,k,0,1), lMesh->getPixyMesh(0,j,k,0,2), lMesh->getPixyMesh(0,j,k,2,2), //17-20
		  lMesh->getPixyMesh(0,j,k,1,1), lMesh->getPixyMesh(0,j,k,2,2), lMesh->getPixyMesh(0,j,k,1,2), mZ, mZ, mZ, //21-26
		  lMesh->getPiZR(0,j,k), lMesh->getPiZPhi(0,j,k), lMesh->getPiRR(0,j,k), lMesh->getPiPhiPhi(0,j,k), lMesh->getPiRPhi(0,j,k), //27-31
		  lMesh->getT(0,j,k));
	}
	fprintf(fTQMR,"\n");
}

void CHydro::openOscarHyper() {
  // 1
//  fOscarHyper = fopen("fOscarHyper.OSCAR2008H","w");
	string fn = parameter::getS(*pMap,"HYDRO_IO_OSCARHYPER_FN","");
	fn = mDataRoot + fn;
	fOscarHyper = fopen( fn.c_str(),"w");
	fprintf(fOscarHyper,"OSCAR2008H  ");
	
  if ( parameter::getD(*pMap,"HYDRO_SVRATIO",0.0) > 0. || parameter::getD(*pMap,"HYDRO_BVRATIO",0.0) > 0.) 
    fprintf(fOscarHyper,"viscous     ");
  else 
    fprintf(fOscarHyper,"ideal       ");
  fprintf(fOscarHyper,"final_hs     \n");
  
  // 2
  fprintf(fOscarHyper,"INIT: ");
  int mA = parameter::getD(*pMap,"GLAUBER_A",197.); int mB = parameter::getD(*pMap,"GLAUBER_B",0.);
  fprintf(fOscarHyper," Glauber - rho = %g Gev/fm^3, xi = %g, sigma = %fb, A = %f, b = %f fm\n",
			parameter::getD(*pMap,"GLAUBER_Rho0",0.16)*1.,parameter::getD(*pMap,"GLAUBER_Xi",0.3)*1.,
			parameter::getD(*pMap,"GLAUBER_Sigma",0.4)*10.,parameter::getD(*pMap,"GLAUBER_A",197.),parameter::getD(*pMap,"GLAUBER_B",0.0));

  // 3
  fprintf(fOscarHyper,"EOS: ");
  if ( parameter::getB(*pMap,"EQOFST_LATEOS",true))
    fprintf(fOscarHyper," Lat-HRG interpolated, 0 Mev Shift, p4 action\n");
  else 
    fprintf(fOscarHyper," ideal gas of massless pions\n");
	
  // 4
  fprintf(fOscarHyper,"CHARGES: ");
  fprintf(fOscarHyper,"none\n");
  
  // 5
  double temp = parameter::getD(*pMap,"HYDRO_FOTEMP",0.13);
  fprintf(fOscarHyper,"HYPER:  hyper T=%g MeV isotherm \n",temp);
  
  // 6
  fprintf(fOscarHyper,"GEOM: ");
  if ( parameter::getB(*pMap,"HYDRO_BJORKEN",true)) 
    if (parameter::getB(*pMap,"HYDRO_PURE_BJORKEN",false)) 
	  fprintf(fOscarHyper,"scaling2d\n");
	else 
	  fprintf(fOscarHyper,"3d\n");
  else 
    if (parameter::getB(*pMap,"HYDRO_PURE_BJORKEN",false)) 
	  fprintf(fOscarHyper,"2d\n");
	else 
	  fprintf(fOscarHyper,"3d-cart");
 
  // 7
  fprintf(fOscarHyper,"GRID:  Euler\n");
  if (parameter::getB(*pMap,"HYDRO_OCTANT",true)) {
    fprintf(fOscarHyper,"%d %d %d %d 0 5 3\n",
			parameter::getI(*pMap,"HYDRO_TSTEP",100),parameter::getI(*pMap,"HYDRO_XSIZE",60)+1,
			parameter::getI(*pMap,"HYDRO_YSIZE",60)+1,parameter::getI(*pMap,"HYDRO_NSIZE",20)+1);
    fprintf(fOscarHyper,"%f ",parameter::getD(*pMap,"HYDRO_T0",1.0));
    if (parameter::getB(*pMap,"HYDRO_LINT",false)) 
	  fprintf(fOscarHyper,"%f ",parameter::getI(*pMap,"HYDRO_TSTEP",100)*parameter::getD(*pMap,"HYDRO_DT",0.1)+parameter::getD(*pMap,"HYDRO_T0",1.0));
    // else FIXME
  
    fprintf(fOscarHyper,"0. %f 0. %f 0. %f\n",parameter::getD(*pMap,"HYDRO_DX",0.1)*parameter::getI(*pMap,"HYDRO_XSIZE",60),
											  parameter::getD(*pMap,"HYDRO_DY",0.1)*parameter::getI(*pMap,"HYDRO_YSIZE",60),
											  parameter::getD(*pMap,"HYDRO_DN",0.1)*parameter::getI(*pMap,"HYDRO_NSIZE",20));
  }
  else {printf("!!!!!!!!!! OSCAR OUTPUT FAILING !!!!!!!!!!"); return;}
  
  // 8
  fprintf(fOscarHyper,"VISCOSITY:  ");
  if ( parameter::getD(*pMap,"HYDRO_SVRATIO",0.0) > 0.)
    if (parameter::getD(*pMap,"HYDRO_BVRATIO",0.0) > 0.) 
	  fprintf(fOscarHyper,"shear and bulk viscosity\n");
	else 
	  fprintf(fOscarHyper,"shear viscosity only\n");
  else 
    if (parameter::getD(*pMap,"HYDRO_BVRATIO",0.0) > 0.) 
	  fprintf(fOscarHyper,"bulk viscosity only\n");
	else 
	  fprintf(fOscarHyper,"none");
  
//  fprintf(fOscarHyper,"COMM:  Dissipative Quantities: eta, tau_pi, s\n");
  fprintf(fOscarHyper,"COMM:  Reflective symmetry in x,y, and eta assumed\n");
  
  fprintf(fOscarHyper,"END_OF_HEADER\n");
}

void CHydro::openOscarFull() {

  // 1
  fOscarFull = fopen("fOscarFull.OSCAR2008H","w");
  fprintf(fOscarFull,"OSCAR2008H  ");
  if ( parameter::getD(*pMap,"HYDRO_SVRATIO",0.0) > 0. || parameter::getD(*pMap,"HYDRO_BVRATIO",0.0) > 0.) 
    fprintf(fOscarFull,"viscous     ");
  else 
    fprintf(fOscarFull,"ideal       ");
  fprintf(fOscarFull,"history     \n");
  
  // 2
  fprintf(fOscarFull,"INIT: ");
  int mA = parameter::getD(*pMap,"GLAUBER_A",197.); int mB = parameter::getD(*pMap,"GLAUBER_B",197.);
  fprintf(fOscarFull," Glauber - rho = %g Gev/fm^3, xi = %g, sigma = %fb, A = %f, b = %f fm\n",
			parameter::getD(*pMap,"GLAUBER_Rho0",0.16)*1.,parameter::getD(*pMap,"GLAUBER_Xi",0.3)*1.,
			parameter::getD(*pMap,"GLAUBER_Sigma",0.4)*10.,parameter::getD(*pMap,"GLAUBER_A",197.),parameter::getD(*pMap,"GLAUBER_B",0.0));
  
  // 3
  fprintf(fOscarFull,"EOS: ");
  if (parameter::getB(*pMap,"EQOFST_LATEOS",true))
    fprintf(fOscarFull," Lat-HRG interpolated, 0 Mev Shift, p4 action\n");
  else 
    fprintf(fOscarFull," ideal gas of quarks and gluons\n");
	
  // 4
  fprintf(fOscarFull,"CHARGES: ");
  fprintf(fOscarFull,"none\n");
  	
  // 5
  fprintf(fOscarFull,"HYPER:  full evolution\n");

  // 6
  fprintf(fOscarFull,"GEOM: ");
  if (parameter::getB(*pMap,"HYDRO_BJORKEN",true)) 
    if (parameter::getB(*pMap,"HYDRO_PURE_BJORKEN",false)) 
	  fprintf(fOscarFull,"scaling2d\n");
	else 
	  fprintf(fOscarFull,"3d\n");
  else 
    if (parameter::getB(*pMap,"HYDRO_PURE_BJORKEN",false)) 
	  fprintf(fOscarFull,"2d\n");
	else fprintf(fOscarFull,"3d-cart");
 
  // 7
  fprintf(fOscarFull,"GRID:  Euler\n");
  if (parameter::getB(*pMap,"HYDRO_OCTANT",true)) {
	fprintf(fOscarFull,"%d %d %d %d 0 5 3\n",
			parameter::getI(*pMap,"HYDRO_TSTEP",100),parameter::getI(*pMap,"HYDRO_XSIZE",60)+1,
			parameter::getI(*pMap,"HYDRO_YSIZE",60)+1,parameter::getI(*pMap,"HYDRO_NSIZE",20)+1);
    fprintf(fOscarFull,"%f ",parameter::getD(*pMap,"HYDRO_T0",1.0));
    if (parameter::getB(*pMap,"HYDRO_LINT",false)) 
	  fprintf(fOscarFull,"%f ",parameter::getI(*pMap,"HYDRO_TSTEP",100)*parameter::getD(*pMap,"HYDRO_DT",0.1)+parameter::getD(*pMap,"HYDRO_T0",1.0));
    // else FIXME
  
	fprintf(fOscarFull,"0. %f 0. %f 0. %f\n",parameter::getD(*pMap,"HYDRO_DX",0.1)*parameter::getI(*pMap,"HYDRO_XSIZE",60),
											 parameter::getD(*pMap,"HYDRO_DY",0.1)*parameter::getI(*pMap,"HYDRO_YSIZE",60),
											 parameter::getD(*pMap,"HYDRO_DN",0.1)*parameter::getI(*pMap,"HYDRO_NSIZE",20));
  }
  else {printf("!!!!!!!!!! OSCAR OUTPUT FAILING !!!!!!!!!!"); return;}
  
  // 8
  fprintf(fOscarFull,"VISCOSITY:  ");
  if (parameter::getD(*pMap,"HYDRO_SVRATIO",0.0) > 0.)
    if (parameter::getD(*pMap,"HYDRO_BVRATIO",0.0) > 0.) 
	  fprintf(fOscarFull,"shear and bulk viscosity\n");
	else fprintf(fOscarFull,"shear viscosity only\n");
  else 
    if (parameter::getD(*pMap,"HYDRO_BVRATIO",0.0) > 0.) 
	  fprintf(fOscarFull,"bulk viscosity only\n");
	else fprintf(fOscarFull,"none");
  
  fprintf(fOscarFull,"COMM:  Dissipative Quantities: a_i (1-5, if SV>0), b (if BV>0)\n");
  fprintf(fOscarFull,"COMM:  Dissipative Qunatities: above will be rescaled to have units GeV/fm^3 (if necessary)\n");
  fprintf(fOscarFull,"COMM:  Transport Coefficients: eta, tau_pi, s\n");
  fprintf(fOscarFull,"COMM:  Reflective symmetry in x,y, and eta assumed\n");
  
  fprintf(fOscarFull,"END_OF_HEADER\n");
}

void CHydro::printOscarFull(CMesh* lMesh, int mT) {
  int mX = lMesh->getXSize();
  int mY = lMesh->getYSize();
  int mN = lMesh->getNSize();
  
  double mSVRatio = parameter::getD(*pMap,"HYDRO_SVRATIO",0.0);
  double mBVRatio = parameter::getD(*pMap,"HYDRO_BVRATIO",0.0);
  
  if (parameter::getB(*pMap,"HYDRO_PURE_BJORKEN",false))
    for (int i=0;i<mX;i++)
      for (int j=0;j<mY;j++) {
	    lMesh->update(0,i,j);
		fprintf(fOscarFull,"%d %d %d ",mT,i,j);
		fprintf(fOscarFull,"%f %f %f 1. ",lMesh->getE(0,i,j),lMesh->getP(0,i,j),lMesh->getT(0,i,j));
		
		fprintf(fOscarFull,"%f %f ",lMesh->getS(0,i,j,1)/lMesh->getS(0,i,j,0), lMesh->getS(0,i,j,2)/lMesh->getS(0,i,j,0));
		
		if (mSVRatio > 0.) {
		  double mAlpha = lMesh->getISAlpha(0,i,j);
		  for (int m=1;m<6;m++) 
		    fprintf(fOscarFull,"%f ",mAlpha*lMesh->getS(0,i,j,m+4));
		}
		if (mBVRatio > 0.) 
		  fprintf(fOscarFull,"%f ",lMesh->getS(0,i,j,10)*lMesh->getISGamma(0,i,j));
		
	    if (mSVRatio > 0.)  
		  fprintf(fOscarFull,"%f %f",lMesh->getSV(0,i,j),lMesh->getTIS(0,i,j));
		if (mBVRatio > 0.)
		  fprintf(fOscarFull,"%f %f",lMesh->getBV(0,i,j),lMesh->getTISB(0,i,j));
		
		fprintf(fOscarFull,"%f",lMesh->getS(0,i,j));
		fprintf(fOscarFull,"\n");
	  }
  else 
    for (int i=0;i<mX;i++)
      for (int j=0;j<mY;j++) 
	    for (int k=0;k<mN;k++) {
	    lMesh->update(k,i,j);
		fprintf(fOscarFull,"%d %d %d %d ",mT,i,j,k);
		fprintf(fOscarFull,"%g %g %g 1. ",lMesh->getE(k,i,j),lMesh->getP(k,i,j),lMesh->getT(k,i,j));
		
		fprintf(fOscarFull,"%g %g %g ",lMesh->getS(k,i,j,1)/lMesh->getS(k,i,j,0),
									   lMesh->getS(k,i,j,2)/lMesh->getS(k,i,j,0), 
									   lMesh->getS(k,i,j,3)/lMesh->getS(k,i,j,0));
		
		if (mSVRatio > 0.) {
		  double mAlpha = lMesh->getISAlpha(k,i,j);
		  for (int m=1;m<6;m++) 
		    fprintf(fOscarFull,"%g ",mAlpha*lMesh->getS(k,i,j,m+4));
		}
		if (mBVRatio > 0.) 
		  fprintf(fOscarFull,"%g ",lMesh->getS(k,i,j,10)*lMesh->getISGamma(k,i,j));
		
	    if (mSVRatio > 0.)  
		  fprintf(fOscarFull,"%g %g ",lMesh->getSV(k,i,j),lMesh->getTIS(k,i,j));
		if (mBVRatio > 0.)
		  fprintf(fOscarFull,"%g %g ",lMesh->getBV(k,i,j),lMesh->getTISB(k,i,j));
		
		fprintf(fOscarFull,"%g ",lMesh->getS(k,i,j));
		fprintf(fOscarFull,"\n");
		fflush(fOscarFull);
	  }
}

void CHydro::printOscarHyper(CMesh* lMesh, int mT) {
// it ix [iy iz] e p T R_qgp vx [vy y_L] [n(1)...n(C)] [mu(1)...mu(C)] dsig_t dsig_x [dsig_y dsig_eta] [Diss(1)...Diss(D)]  [Tr(1)...Tr(T)]
  if (mPureBjorken) {
  	double mFOS[XSIZE+YSIZE][2];
	double mFOSigma[XSIZE+YSIZE][2];
	double mFOVelo[XSIZE+YSIZE][2];
	double mFODiss[XSIZE+YSIZE][5];
	int fosSize;
	lMesh->getFOS(mFOS,mFOSigma,mFOVelo,mFODiss,fosSize);
	double temp = mFoTemp;
	double ed = mEos->getEGivenT(temp);

	fprintf(fOscarHyper,"time:  %g\n",lMesh->getTau());
	
	for (int i=0;i<fosSize;i++) {

//	  fprintf(fOscarHyper,"%g %g %g ", lMesh->getTau(), mFOS[i][0], mFOS[i][1]);
	  fprintf(fOscarHyper,"%g %g ", mFOS[i][0], mFOS[i][1]);
	  fprintf(fOscarHyper,"%g %g %g 1. ",ed,mEos->getP(ed),temp);
	  for (int k=0;k<2;k++)
	    fprintf(fOscarHyper,"%g ",mFOVelo[i][k]);
	  fprintf(fOscarHyper,"0. %g %g ", /*lMesh->getTau(),*/ mFOSigma[i][0], mFOSigma[i][1]);
	  for (int k=0;k<5;k++)
	    fprintf(fOscarHyper,"%g ",mFODiss[i][k]);
	  fprintf(fOscarHyper,"\n");
	}
  }
}

void CHydro::printDNs(CMesh* lMesh, int mNSize, int mXSize, int mYSize) {

/*
  double mDx = parameter::getD(*pMap,"DX");
  double mDy = parameter::getD(*pMap,"DY");
  double mDn = parameter::getD(*pMap,"DN");

  printf("wha??\n"); fflush(stdout);
  double Pi[mNSize][mXSize][mYSize][4][4];
  printf("wha??\n"); fflush(stdout);
  double T[mNSize][mXSize][mYSize];
  printf("wha??\n"); fflush(stdout);
  double E[mNSize][mXSize][mYSize];
  printf("wha??\n"); fflush(stdout);
  double P[mNSize][mXSize][mYSize];
/*
  double Pi[NSIZE][XSIZE][YSIZE][4][4];
  double T[NSIZE][XSIZE][YSIZE];
  double E[NSIZE][XSIZE][YSIZE];
  double P[NSIZE][XSIZE][YSIZE];

  printf("wha??\n"); fflush(stdout);
  for (int i=0;i<mNSize;i++)
    for (int j=0;j<mXSize;j++)
	  for (int k=0;k<mYSize;k++) {
	    for (int m=0;m<4;m++)
		  for (int n=0;n<4;n++) 
		    Pi[i][j][k][m][n] = lMesh->getPixyMesh(i,j,k,m,n);
		T[i][j][k] = lMesh->getT(i,j,k);
		E[i][j][k] = lMesh->getE(i,j,k);
		P[i][j][k] = lMesh->getP(i,j,k);
	  }
*/
  double mTemp;
//  double coeff = 72.*mDx*mDy*mDn*lMesh->getTau()*JIFPI*JIFPI*JIFPI;
  double dY=0.2;
  double dPT = 0.2;
  double lowPt = 0.;
  double m = 0.14;
  
  double dNd2PdY[21][11][11];
/*  
  double pu,puu;
  
  for (int i=0;i<=10;i++) {
    double px = dPT*i;
    for (int j=0;j<=10;j++) {
	  double py = dPT*j;
	  for (int k=0;k<=20;k++) {
	    double Y = k*dY;
		double mTemp = 0.;
        for (int eta=0;eta<=mNSize;eta++) {
		  double p0 = sqrt(m*m + px*py + py*py) * cosh(Y - mDn*eta);
		  double pz = p0 * tanh(Y - mDn*eta);
	
		  for (int x=0;x<=mXSize;x++)
			for (int y=0;y<=mYSize;y++) {
			  pu =  p0*lMesh->getS(eta,x,y,0) - px*lMesh->getS(eta,x,y,1) 
				  - py*lMesh->getS(eta,x,y,2) - pz*lMesh->getS(eta,x,y,3);

			  puu =      Pi[eta][x][y][0][0]*p0*p0;
			  puu -=     Pi[eta][x][y][0][1]*p0*px + Pi[eta][x][y][0][2]*p0*pz + Pi[eta][x][y][0][3]*p0*pz;
			  puu +=     Pi[eta][x][y][1][1]*px*px + Pi[eta][x][y][2][2]*py*py + Pi[eta][x][y][3][3]*pz*pz;
			  puu += 2.*(Pi[eta][x][y][1][2]*px*py + Pi[eta][x][y][1][3]*px*pz + Pi[eta][x][y][2][3]*py*pz);
		
			  mTemp += p0 * exp( - pu/T[eta][x][y]) 
						  * (1. + 0.5*puu/((E[eta][x][y]+P[eta][x][y])*T[eta][x][y]*T[eta][x][y]));
		    }
	    }
		dNd3p[k][i][j] = coeff * mTemp;
		printf("%0.6g %0.6g %0.6g %0.6g\n",px,py,Y,dNd3p[k][i][j]);
	  }
	  printf("\n");
	}
  }
*/

  printf("\n");
  for (int i=0;i<=10;i++)
    for (int j=0;j<=10;j++) {
	  for (int k=0;k<=20;k++) {
	    dNd2PdY[k][i][j] = 9.*lMesh->getdNd3p(lowPt+dPT*i, lowPt+dPT*j, dY*k, m);
		printf("%0.6g %0.6g %0.6g %0.6g\n",lowPt+dPT*i, lowPt+dPT*j, dY*k, dNd2PdY[k][i][j]);
	  }
	  printf("\n"); 
	}

  double dNdY[20];
  fPiDNDY = fopen("fPiDNDY.txt","w");

  for (int k=0;k<=20;k++) {
    dNdY[k] = 0.;
    for (int i=0;i<=10;i++)
      for (int j=0;j<=10;j++)
	    dNdY[k] += dPT*dPT*dNd2PdY[k][i][j];
	fprintf(fPiDNDY,"%0.6g %0.6g\n", dY*k, dNdY[k]);
  }
 
  for (int i=0;i<=10;i++) {
    
  }
  
 
// fPiDNDPt= fopen("fPiDNDPt.txt","w");

}

string CHydro::mDataRoot, CHydro::mOldFileName;
bool CHydro::mDebug, CHydro::mIoIntegrals, CHydro::mIoSlices, CHydro::mIoSpots, CHydro::mIoTechqm;
bool CHydro::mIoOscarFull, CHydro::mIoOscarHyper, CHydro::mIoFull, CHydro::mIoSpectra, CHydro::mOldFile;
bool CHydro::mSawtooth, CHydro::mRK2, CHydro::mRK4, CHydro::mLinT, CHydro::mLogT, CHydro::mLogSinhT;
bool CHydro::mBjorken, CHydro::mPureBjorken, CHydro::mHalving, CHydro::mOctant;
int CHydro::mTStep, CHydro::mXSize, CHydro::mYSize, CHydro::mNSize, CHydro::mIoSliceTStep, CHydro::mIoTechqmTStep;
int CHydro::mIoOscarTStep, CHydro::mPrintStep, CHydro::mPrintX, CHydro::mPrintY, CHydro::mPrintN;
double CHydro::mInitNS, CHydro::mInitFlow, CHydro::mFoTemp, CHydro::mT0, CHydro::mDt, CHydro::mE0;
double CHydro::intS, CHydro::intE, CHydro::intP[3], CHydro::sLoss, CHydro::sLoss2;
double CHydro::eLoss, CHydro::intS0, CHydro::intE0, CHydro::intP0[3], CHydro::pLoss[3], CHydro::epLoss;