#include "hydro.h"

int main (int argc, char * const argv[]) {
	time(&start);	

	pMap = new parameterMap();	
	parameter::ReadParsFromFile(*pMap ,argv[1]);

	for (int i=2; argv[i] != NULL; i++) parameter::ReadParsFromFile(*pMap, argv[i]);

//  parameter::PrintPars(*pMap);

	mDataRoot = parameter::getS(*pMap,"HYDRO_OUTPUT_DATAROOT","");

	int sLength = mDataRoot.size();
	if (mDataRoot[sLength-1] != '/') mDataRoot.append("/");
//  printf("%s\n",mDataRoot.c_str());

	bool mDebug = parameter::getB(*pMap,"HYDRO_DEBUG",false);
	bool mIoIntegrals, mIoSlices, mIoSpots, mIoTechqm, mIoOscarFull, mIoOscarHyper, mIoFull, mIoSpectra;
	mIoIntegrals = parameter::getB(*pMap,"HYDRO_IO_INTEGRALS", false);
	mIoSlices = parameter::getB(*pMap,"HYDRO_IO_SLICES",false);
	mIoSpots = parameter::getB(*pMap,"HYDRO_IO_SPOTS",false);
	mIoTechqm = parameter::getB(*pMap,"HYDRO_IO_TECHQM",false);
	mIoOscarFull = parameter::getB(*pMap,"HYDRO_IO_OSCARFULL",false);
	mIoOscarHyper = parameter::getB(*pMap,"HYDRO_IO_OSCARHYPER",false);
	mIoFull = parameter::getB(*pMap,"HYDRO_IO_FULL",false);
	mIoSpectra = parameter::getB(*pMap,"HYDRO_IO_SPECTRA",false);

	bool mSawtooth, mRK2, mRK4, mLinT, mLogT, mLogSinhT;
	mSawtooth = parameter::getB(*pMap,"HYDRO_SAWTOOTH",false);
	mRK2  = parameter::getB(*pMap,"HYDRO_RK2",true);
	mRK4  = parameter::getB(*pMap,"HYDRO_Rk4",false);
	mLinT = parameter::getB(*pMap,"HYDRO_LINT",false);
	mLogT = parameter::getB(*pMap,"HYDRO_LOGT",false);
	mLogSinhT = parameter::getB(*pMap,"HYDRO_LOGSINHT",false);

	int mTStep, mXSize, mYSize, mNSize, mIoSliceTStep, mIoTechqmTStep, mIoOscarTStep, mPrintStep, mPrintX, mPrintY, mPrintN;
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

	bool mBjorken, mPureBjorken, mHalving;
	mBjorken = parameter::getB(*pMap,"HYDRO_BJORKEN",true);
	mPureBjorken = parameter::getB(*pMap,"HYDRO_PURE_BJORKEN",false);
	mHalving = parameter::getB(*pMap, "HYDRO_HALVING", false);

	double mInitNS, mInitFlow;
	mInitNS = parameter::getD(*pMap,"HYDRO_INIT_NS",0.);
	mInitFlow = parameter::getD(*pMap,"HYDRO_INIT_FLOW",0.);

	bool mOldFile = parameter::getB(*pMap,"HYDRO_OLD_FILE",false);
	string mOldFileName = parameter::getS(*pMap,"HYDRO_OLD_FILENAME","");

	double mFoTemp, mT0, mDt, mE0;
	mFoTemp = parameter::getD(*pMap,"HYDRO_FOTEMP",0.13);
	mT0 = parameter::getD(*pMap,"HYDRO_T0",1.0);
	mDt = parameter::getD(*pMap,"HYDRO_DT",0.03);
	mE0 = parameter::getD(*pMap,"HYDRO_E0",1.0);

	if (mIoIntegrals) openFileIntegrals();
	if (mIoSpots) openFileSpots();
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

// midpoint method... 
	CMesh *onMesh, *offMesh, *tempMesh, *k1, *k2;
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
	}
	CMesh* deadMesh;// just the pointer, for switching... = new CMesh()

	CEos*  mEos = new CEos(pMap);
	onMesh->setEos(mEos);  // static, so just once
	onMesh->setParameterMap(pMap);
	onMesh->deaden();

	if (mInitNS != 0.)   onMesh->initNS();
	if (mInitFlow != 0.) onMesh->addInitialFlow();

	double intS,intE,intP[3];
	double sLoss=0.,sLoss2=0.,eLoss=0.;
	double intS0,intE0,intP0[3];
	double pLoss[3] = {0., 0., 0.};
	double epLoss=0.;

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
		if (mBjorken) 
			fprintf(fE0,"%0.9g %0.9g\n", mTau, pow(mTau,4./3.));
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
			t = mTStep;
			printE(onMesh,5);
			t++;
		}

		if (mIoIntegrals) {
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
		}

		if (tempMesh->getT(0,0,0) < mFoTemp) {
			if (mIoSlices || mIoFull) printE(tempMesh,5);
			t=mTStep+1;
		}

		if (t%mPrintStep==0 && (mIoSlices || mIoFull)) printE(tempMesh,t/mPrintStep);
		/*      if ( tempMesh->getTau() > 1. && onMesh->getTau() < 1.)		 printE(tempMesh,mPrintStep);
		else if ( tempMesh->getTau() > 1.6 && onMesh->getTau() < 1.6) printE(tempMesh,2*mPrintStep);
		else if ( tempMesh->getTau() > 2.0 && onMesh->getTau() < 2.)  printE(tempMesh,3*mPrintStep);
		else if ( tempMesh->getTau() > 5. && onMesh->getTau() < 5.)   printE(tempMesh,4*mPrintStep);
		else if ( t==TSTEP) printE(tempMesh,5*mPrintStep);
		/*	  if	  (t==18) printE(tempMesh,mPrintStep);
		else if (t==19) printE(tempMesh,2*mPrintStep);
		else if (t==20) printE(tempMesh,3*mPrintStep);
		else if (t==21) printE(tempMesh,4*mPrintStep);
		else if (t==22) printE(tempMesh,5*mPrintStep);*/
			else if ((t*50)%mTStep==0) {printf(".");  fflush(stdout);}

		if (t%mIoSliceTStep == 0) {
			double mTau = tempMesh->getTau();

			if (mIoIntegrals) {
				fprintf(fIE,"%0.9g %0.9g\n", mTau, (intE+eLoss)/intE0 - 1.);
				if (intP0[0] != 0.) fprintf(fIPX,"%0.9g %0.9g\n", mTau, (intP[0]+pLoss[0])/intP0[0] - 1.);
				else				fprintf(fIPX,"%0.9g %0.9g\n", mTau,  intP[0]+pLoss[0]);
				if (intP0[1] != 0.) fprintf(fIPY,"%0.9g %0.9g\n", mTau, (intP[1]+pLoss[1])/intP0[1] - 1.);
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
			if (mIoSpots) {
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
		}

		if (mIoTechqm && t%mIoTechqmTStep == 0)
			printTQM(tempMesh);

		if (mIoOscarFull && t%mIoOscarTStep == 0)
			printOscarFull(tempMesh,t);

		if (mIoOscarHyper && t%mIoOscarTStep == 0)
			printOscarHyper(tempMesh,t,mEos);

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
	k1->~CMesh();
	k2->~CMesh();;

	if (mIoSpectra) {
		printf("Computing Spectra.......\n"); fflush(stdout);
		printDNs(tempMesh,mNSize,mXSize,mYSize,pMap);
	}

	tempMesh->~CMesh();

//  system("cp def.h def.txt");

//  string temp = mDataRoot;
//  temp.append(argv[2]);

	return 0;
}
