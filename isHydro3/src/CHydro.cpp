#include "CHydro.h"

CHydro::CHydro(parameterMap* pM){
	pMap = pM;

	mDataRoot = parameter::getS(*pMap,"HYDRO_OUTPUT_DATAROOT","./");
	int sLength = mDataRoot.size();
	if (mDataRoot[sLength-1] != '/') mDataRoot.append("/");
  
	if (!checkDirectoryExistence(mDataRoot)){
		printf("\nHYDRO_OUTPUT_DATAROOT\n (((%s))) does not exist! Aborting.\n",
		mDataRoot.c_str());
		exit(1);
	}
	
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
	mRK4  = parameter::getB(*pMap,"HYDRO_RK4",false);
	mLinT = parameter::getB(*pMap,"HYDRO_LINT",false);
	mLogT = parameter::getB(*pMap,"HYDRO_LOGT",false);
	mLogSinhT = parameter::getB(*pMap,"HYDRO_LOGSINHT",true);

	mTStep = parameter::getI(*pMap,"HYDRO_TSTEP",-1);
	mXSize = parameter::getI(*pMap,"HYDRO_XSIZE",60);
	mYSize = parameter::getI(*pMap,"HYDRO_YSIZE",60);
	mNSize = parameter::getI(*pMap,"HYDRO_NSIZE",20);
	mPrintX = parameter::getI(*pMap,"HYDRO_PRINTX",0);
	mPrintY = parameter::getI(*pMap,"HYDRO_PRINTY",0);
	mPrintN = parameter::getI(*pMap,"HYDRO_PRINTN",0);

	mIoSliceTStep  = parameter::getI(*pMap,"HYDRO_IO_SLICES_TSTEP",10);
	mIoTechqmTStep = parameter::getI(*pMap,"HYDRO_IO_TECHQM_TSTEP",10);
	mIoOscarTStep  = parameter::getI(*pMap,"HYDRO_IO_OSCAR_TSTEP",10);
	mIoSpotsTStep  = parameter::getI(*pMap,"HYDRO_IO_SPOTS_TSTEP",1);
	mIoFullTStep   = parameter::getI(*pMap,"HYDRO_IO_FULL_TSTEP",100);
	
	mBjorken = parameter::getB(*pMap,"HYDRO_BJORKEN",true);
	mPureBjorken = parameter::getB(*pMap,"HYDRO_PURE_BJORKEN",false);
	mHalving = parameter::getB(*pMap, "HYDRO_HALVING", false);
	mOctant = parameter::getI(*pMap, "HYDRO_OCTANT", 3);

	mInitNS = parameter::getD(*pMap,"HYDRO_INIT_NS",0.);
	mInitNSTxxP = parameter::getD(*pMap,"HYDRO_INIT_NS_TXXP",1.);
	
	mInitFlow = parameter::getD(*pMap,"HYDRO_INIT_FLOW",0.);

	mOldFile = parameter::getB(*pMap,"HYDRO_OLD_FILE",false);
	mOldFileName = parameter::getS(*pMap,"HYDRO_OLD_FILENAME","");

	mFoTemp = parameter::getD(*pMap,"HYDRO_FOTEMP",0.165);
	mT0 = parameter::getD(*pMap,"HYDRO_T0",1.0);
	mDt = parameter::getD(*pMap,"HYDRO_DT",0.03);
	mE0 = parameter::getD(*pMap,"HYDRO_E0",1.0);

	
	nPrintSlice = 0; nPrintFull = 0;
	
	initializeHydro();
	
	//	zeroPointers();
}

CHydro::~CHydro(){
	if(onMesh!=NULL){
		delete onMesh;
		onMesh=NULL;
	}
	if(offMesh!=NULL){
		delete offMesh;
		offMesh=NULL;
	}
	if(tempMesh!=NULL){
		delete tempMesh;
		tempMesh=NULL;
	}
	delete mEos;
	
	if (mRK4) {
		delete k1;
		delete k2;
	}
}

void CHydro::zeroPointers(){
	fEX = NULL;
	fEN = NULL;
}

int CHydro::runHydro(){
	time(&start);
	time(&now);
	
	// time to complete estimate
	if (!mPureBjorken) 
		printf("Expected Run Time < %0.3g sec.  \n", 
	(mXSize*mYSize*mNSize*150.)/2.7E5);
	else 
		printf("Expected Run Time < %0.3g sec.  \n",
	(mXSize*mYSize*150)/5.E5);
	fflush(stdout);
	printf("Time Step 0 [%0.6g](elapsed time %0.6g sec) .",
	onMesh->getTau(),difftime(now,start)); fflush(stdout);
	
	nPrintUpdate=1;
	mTemp0 = tempMesh->getT(0,0,0);
	int t=1;
	
	if (mTemp0 < mFoTemp) 
		t = mTStep+1;
	printf("mTStep=%d, mTemp0=%g, mFoTemp=%g\n",mTStep,mTemp0,mFoTemp);
	for (;(t<=mTStep || mTStep==-1);t++) {
		onMesh->deaden();

		if (mSawtooth) 
			onMesh->forward(tempMesh,offMesh,mDt);
		else if (mRK2) {
			onMesh->copyActive(offMesh,tempMesh);

			onMesh->forward(offMesh,mDt/2.);

			/*
			CMesh* m = tempMesh;
			tempMesh = offMesh;
			printSlices();
			tempMesh = m;
			*/
			
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
			if (mIoSlices)
				printSlices();
			
			if (mIoSpots)
				closeFileSpots();
			
			if (mIoIntegrals)
				closeFileIntegrals();
			
			if (mIoOscarHyper)
				fclose(fOscarHyper);
			
			if (mIoOscarFull)
				fclose(fOscarFull);
			
			return 1;
		}

		if (tempMesh->getT(0,0,0) < mFoTemp) 
			if (!tempMesh->anyActiveCells())	
				break;

		if (mIoIntegrals) 
			printIntegrals();

		if (mIoSlices && t%mIoSliceTStep==0) 
			printSlices();

		/*
		if (tempMesh->getTau() > 1.6 && tempMesh->getTau() < 1.6 + mDt) {
		printf("\n1.6 print is number %d (%g)\n",nPrintSlice,tempMesh->getTau());
		printSlices();
		}
		*/

		if (mIoFull && t%mIoFullTStep==0) 
			printMesh();
		
		if (mIoSpots && t%mIoSpotsTStep == 0) 
			printSpots();

		if (mIoTechqm && t%mIoTechqmTStep == 0)
			printTQM();

		if (mIoOscarFull && t%mIoOscarTStep == 0)
			printOscarFull(t);

		if (mIoOscarHyper && t%mIoOscarTStep == 0){
			tempMesh->genFOS(oldMesh);
			oldMesh->copyMesh(tempMesh);
			
			if (!printOscarHyper(t)){
				printf("\nILLICIT FREEZEOUT SURFACE at t=%g\n",tempMesh->getTau());
				mIoOscarHyper=false;  // quit trying to print FOS
			}
		}
		
		if (mHalving)
		if (tempMesh->getTau() - onMesh->getTau() > 2.*mDt*mT0) {
			printf("!");
			mDt /= 2.;
		}

		if ( nPrintUpdate < 50.*(pow(tempMesh->getT(0,0,0),-2)-pow(mTemp0,-2))/(pow(mFoTemp,-2)-pow(mTemp0,-2))) {
			if (nPrintUpdate%10==0){
				time(&now);
				printf("\nAt Time %d: [%g] -- elapsed time %g sec .",t,tempMesh->getTau(),difftime(now,start)); 
				fflush(stdout);
			}
			else {
				printf(".");
				fflush(stdout);
			}
			nPrintUpdate++;
		}

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
	// if we quit early, undo mesh reassignments
	if (mRK2 && tempMesh->getT(0,0,0) > mFoTemp) {
		deadMesh = tempMesh;
		tempMesh = onMesh;
		onMesh = deadMesh;
	}

	if (mIoSlices)
		printSlices();
	
	if (mIoFull){
		printMesh();
#ifdef HYDRO_XDMF_OUTPUT
		StopXdmf();
#endif
	}
	
	if (mIoTechqm)
		printTQM();
	
	if (mIoSpots) {
		printSpots();
		closeFileSpots();
	}

	if (mIoIntegrals)
		closeFileIntegrals();

	if (mIoOscarHyper)	{
		fclose(fOscarHyper);
		if (tempMesh->getT(0,0,0) < mFoTemp)
			tempMesh->genFOS(oldMesh);
		fclose(fOSU);
		delete oldMesh;
		oldMesh=NULL;
	}

	if (mIoOscarFull) {
		printOscarFull(t);
		fclose(fOscarFull);
	}
	
	time(&now);
	printf("\nHydro Finished %d [%g]...... (Total Time: %0.3g sec)\n",
	t,tempMesh->getTau(),difftime(now,start)); 
	fflush(stdout);
	
	if (mIoSpectra) {
		time(&start);
		printf("\nCalculating Spectra (Exp Time %0.3g s)....",mXSize*mYSize*mNSize*1500./7.E6); fflush(stdout);
		
		string fn = "fdNdY.txt";
		fn = mDataRoot + fn;
		FILE *fdNdY = fopen( fn.c_str(),"w");
		
		const int nRap = 10;
		double dRap = 0.4;   
		double dNdY[nRap];
		
		tempMesh->calcdNdY(0.1396, 1.0, dNdY, dRap, nRap);
		
		for (int i=0;i<nRap;i++)
			fprintf(fdNdY,"%g %g\n",dRap*i,dNdY[i]);
		
		time(&now);
		printf("Finished! [%g s]\n",difftime(now,start));
	}
	
	tempMesh->writeIEPvEta();
	
#ifdef HYDRO_BOOST_THREADS
	delete mHelper;
#endif
	return 0;
}

void CHydro::initializeHydro(){
	// Open IO Files
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
  
	// make main mesh
	//printf("making main mesh...., mOctant=%d\n",mOctant);
	oldMesh=onMesh=offMesh=tempMesh=deadMesh=k1=k2=NULL;
	if (mOldFile) {
		if (mOctant == 3)
			onMesh = new octMesh(pMap,mOldFileName.c_str());
		else 
			onMesh = new fullMesh(pMap,mOldFileName.c_str());
	}
	else {
		if (mOctant == 3)
			onMesh = new octMesh(pMap); 
		else 
			onMesh = new fullMesh(pMap);
	}
	
	mEos = new CEos(pMap);
	onMesh->setEos(mEos);
	
	if (parameter::getB(*pMap,"HYDRO_KLN_INPUT",false))
		onMesh->fixEntropyMesh();

	//	printf("filling active cells...\n");
	onMesh->fillActiveCells();
	
	// alter initial conditions
	if (mInitFlow != 0.) 
		onMesh->addInitialFlow();
	
	if (mInitNS != 0. || mInitNSTxxP != 0.)   
		onMesh->initNS();
	
	// make additional meshes
	//	printf("making additional meshes...\n");
	if (mOctant == 3){
		offMesh  = new octMesh(onMesh);
		tempMesh = new octMesh(onMesh);
		if (mRK4) {
			k1 = new octMesh(onMesh); 
			k2 = new octMesh(onMesh);
		}
	}
	else {
		offMesh  = new fullMesh(onMesh);
		tempMesh = new fullMesh(onMesh);
		if (mRK4) {
			k1 = new fullMesh(onMesh); 
			k2 = new fullMesh(onMesh);
		}
	}
	
	//	printf("active matching between meshes...\n");
	offMesh->fillActiveCells();
	tempMesh->fillActiveCells();
	onMesh->copyActive(offMesh,tempMesh);
	
	if (mRK4) {
		k1->fillActiveCells();
		k2->fillActiveCells();
		onMesh->copyActive(k1,k2);
	}
	
#ifdef HYDRO_BOOST_THREADS
	offMesh->setEos(mEos);
	tempMesh->setEos(mEos);
	CCellHelper* mHelper = new CCellHelper();
	tempMesh->setHelper(mHelper);
	onMesh->setHelper(mHelper);
	offMesh->setHelper(mHelper);
#endif
	
	//printActive();
	
#ifdef HYDRO_XDMF_OUTPUT
	StartXdmf();
	WriteXdmf();
#endif
	
	sLoss=0.; 
	eLoss=0.;
	
	//writeIC();
	
	// first IO dump
	//printf("first i/o dump...\n");
	if (mIoIntegrals) {
		openFileIntegrals();
		onMesh->calcIntegrals();
		double intS = onMesh->integralS();
		intS0 = intS;
		double intE = onMesh->integralE();
		intE0 = intE;

		fprintf(fIS,"%0.9g 0.\n", mT0);
		fprintf(fIE,"%0.9g 0.\n", mT0);
		
		fprintf(fIEP,"0. %0.6g\n", onMesh->integralEp());
		fprintf(fIEX,"0. %0.6g\n", onMesh->integralEx());
		fprintf(fIVR,"0. %0.6g\n", onMesh->averageVr());
	}
	if (mIoSpots) {
		openFileSpots();
		printSpots();
	}
	if (mIoOscarFull) 
		openOscarFull();
	if (mIoOscarHyper) {
		openOscarHyper();
		if (mOctant == 3)
			oldMesh = new octMesh(tempMesh);
		else 
			oldMesh = new fullMesh(tempMesh);

		string osuFN("/osuFOS.dat");
		osuFN = mDataRoot + osuFN;
		fOSU = fopen(osuFN.c_str(),"w");
		oldMesh->setOsuFos(fOSU);
	}
	if (mIoSlices) printSlices();
	if (mIoFull && false)
		printMesh();
	
}

string CHydro::mDataRoot, CHydro::mOldFileName;
bool CHydro::mDebug, CHydro::mIoIntegrals, CHydro::mIoSlices;
bool CHydro::mSawtooth, CHydro::mRK2, CHydro::mRK4;
bool CHydro::mLinT, CHydro::mLogT, CHydro::mLogSinhT;
bool CHydro::mBjorken, CHydro::mPureBjorken, CHydro::mHalving;
int CHydro::mOctant;
int CHydro::mTStep, CHydro::mXSize, CHydro::mYSize, CHydro::mNSize;
double CHydro::mInitNS, CHydro::mInitNSTxxP, CHydro::mInitFlow, CHydro::mFoTemp;
double CHydro::mT0, CHydro::mDt, CHydro::mE0;
double CHydro::intS0, CHydro::intE0;
double CHydro::eLoss, CHydro::sLoss;