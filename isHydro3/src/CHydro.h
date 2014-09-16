#ifndef __CHYDRO_h_INCLUDE__
#define __CHYDRO_h_INCLUDE__

#include "hydroDef.h"
#include <cstdlib>
#include <cstdio> 
#include <cmath>

#include <ctime>
#include <math.h>
#include <string.h>
#include <sys/stat.h> 

#include "coralutils.h"
#include "CMesh.h"
#include "octMesh.h"
#include "fullMesh.h"
#include "CEos.h"

#ifdef HYDRO_XDMF_OUTPUT
#include <H5Cpp.h>
#include <Xdmf.h>
#endif

#ifdef HYDRO_XDMF_OUTPUT
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif
#endif

using namespace std;

class CHydro {
private:
	CMesh *onMesh, *offMesh, *tempMesh, *k1, *k2, *deadMesh,*oldMesh;
	CEos*  mEos;

	FILE *fEX, *fEY, *fEZ, *fS, *fT, *fP;
	FILE *fEN, *fU0, *fUxx, *fUxy, *fUyy, *fUyx, *fUz, *fULz, *fUzz;
	FILE *fE0, *fTxx0, *fTzz0;
	FILE *fIE, *fIS, *fIEP, *fIEX, *fIVR;
	FILE *fSLoss, *fELoss, *fSLoss2;
	FILE *fA1x, *fA2x, *fA3x, *fA4x, *fA5x, *fBx;
	FILE *fA1y, *fA2y, *fA3y, *fA4y, *fA5y, *fBy;
	FILE *fA1z, *fA2z, *fA3z, *fA4z, *fA5z, *fBz;
	FILE *fTxx, *fTyy, *fTzz, *fTxy, *fTxz, *fTyz, *fTxxNS;
	FILE *fTrr, *fTpp, *fTpr;
	FILE *fUL, *fEL, *fA1L, *fA2L, *fBL;
	FILE *fTIS, *fTISB, *fDzz, *fT0;
	FILE *fFOS, *fFOSL, *fFOSSST;
	FILE *fTQMX, *fTQMY, *fTQMR;
	FILE *fTxxL, *fTxxNSL,*fTxxML;
	FILE *fFull, *fOscarFull, *fOscarHyper;
	FILE *fPiDNDY, *fPiDNDPt;
	FILE *fOSU;

#ifdef HYDRO_XDMF_OUTPUT
	H5File *h5outfile, *h5infile;
	XdmfDOM *d;
	XdmfRoot *root;
	XdmfDomain *domain;
	XdmfGrid *tGrid;
	XdmfInt64 sShape[3], vShape[4], tShape[4];
#endif
	
	static string mDataRoot, mOldFileName;
	static bool mDebug, mIoIntegrals, mIoSlices, mIoSpots, mIoTechqm;
	static bool mIoOscarFull, mIoOscarHyper, mIoFull, mIoSpectra, mOldFile;
	static bool mSawtooth, mRK2, mRK4, mLinT, mLogT, mLogSinhT;
	static bool mBjorken, mPureBjorken, mHalving;
	static int mOctant;
	static int mTStep, mXSize, mYSize, mNSize, mPrintX, mPrintY, mPrintN;
	static int mIoSliceTStep, mIoTechqmTStep, mIoOscarTStep, mIoSpotsTStep, mIoFullTStep;
	static double mInitNS, mInitNSTxxP, mInitFlow, mFoTemp, mT0, mDt, mE0;

	static double intS0,intE0,sLoss,eLoss;

	double mTemp0;
	int nPrintSlice, nPrintFull, nPrintUpdate;
	
	void printSlices();
	void printMesh();
	void printIntegrals();
	void printSpots();
	void printTQM();
	bool printOscarHyper(int mT);
	void printOscarFull(int mT);
	void printDNs(CMesh* lMesh,int,int,int);
	
	void openFileIntegrals();
	void openFileSpots();
	void openSliceFiles();
	void openOscarHyper();
	void openOscarFull();
	
	void closeFileIntegrals();
	void closeFileSpots();
	void closeSliceFiles();

	void printActive();
	
	void testFileOpen();
	void zeroPointers();
	
	inline int min(int a, int b) {if (a>b) return b; else return a;}

	time_t start, now;

	parameterMap* pMap;
	
	bool checkDirectory(string fName);
	bool checkDirectoryExistence(string dirname);
	
#ifdef HYDRO_XDMF_OUTPUT
	void fillCellsH5(CCellH5& cH, CCell* p);
	int ReadDataH5();
	int WriteDataH5();
	
	void StartXdmf();
	void WriteXdmf();
	void StopXdmf();
	XdmfArray* xdmfScalar(XdmfAttribute& mAtt, const char name[], XdmfInt64 sShape[3]);
	void MoveH5();
	void ReadXdmf();
#endif
	
public:
	CHydro(parameterMap*);
	~CHydro();
	int runHydro();
	void spectraFromFile();
	void printVorticity();
	void initializeHydro();
	void writeIC(); // summarized initial conditions, and writes to file
	inline void setEos(CEos* eos_ptr) {mEos=eos_ptr;};
};

#endif //__CHYDRO_h_INCLUDE__
