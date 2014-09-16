#ifndef __CMESH_h_INCLUDE__
#define __CMESH_h_INCLUDE__
/*  Cmesh.h
 *  isHydro3
 *  Created by Joshua Vredevoogd on 2/11/09.
 */

#include "hydroDef.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <list>
#ifdef HYDRO_BOOST_THREADS
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/tss.hpp>
#endif
#include <string.h>
#include "coralutils.h"
#include "CCell.h"
#include "CEos.h"
#include "cornelius.h"

class CMesh {
	typedef double (CMesh::*meshMemFunc)(const double x);
	
protected:
	std::list<CCell*> activeCells;
	
	static double ****localE4;
	static double ***localE3;
	static double *dxCorn;
	
#ifdef HYDRO_BOOST_THREADS
	static CEos* eosVector[NTHREADS];
	static CCellHelper* helperVector[NTHREADS];
	std::list<CCell*>::iterator marks[NTHREADS][2];  // break points
#endif
	
	static parameterMap *pMap;
	static Cornelius mCornelius;	
	static string mDataRoot;
	
		//parameters from map
	static int mOctant;
	static bool mPureBjorken, mBjorken, mPrintMs, mSVTrimInit, mSVTrim, mSawtooth;
	static double mE0, mFOTemp, mDeadT, mT0, mDx, mDy, mDn;
	static int mNSize, mXSize, mYSize, mNSizeOrig, mXSizeOrig, mYSizeOrig;
	static int mPrintX, mPrintY, mPrintN;
	
	static bool bGSat,bBGK;
	static double wnRho0, wnRAu, wnXi, wnSigma, wnA, wnB, wnKTau, gSatN;
	static double mRn, mRt, mRa, mRx, mRy, mRSig;
	static double mWnBinRatio, mInitFlow, mInitNS, mInitNSTxxP;
	static double icNiceA, icNiceB;
	
	static double fracPPE, sigmaAtt,ppet;
	static double wnAttRatio, wnTx, wnTy;
	static double collRootS;
	
	static int activeX, activeY, activeN;
	
	static FILE *osuFOS;
	
	double intS, intE, intEp, intEx, intVr;
	double sLoss, eLoss;
	double sDead, eDead, epDead, exDead, vrDead;
	
	bool mFBose, mFFerm;
	
	void fillParams();
	
		// to assign initial conditions to the cells
	void initialCondition(CCell*);
	void initialCondition(CCell*,CCell*);
	
	virtual void klnInput()=0;
	
	virtual void initBGK()=0;
	
		//	double qromb(double func(const double), double a, double b);
		//	double trapzd(double func(const double), const double a, const double b, const int n);
	double qromb(meshMemFunc func, double a, double b);
	double trapzd(meshMemFunc func, const double a, const double b, const int n);
	void polint(std::vector<double> &xa, std::vector<double> &ya, const double x, double &y, double &dy);
	
		// for 2d Glauber-type initial condition generation
	double wnE(double, double);
	
	inline double wnRhoR(double r){return wnRho0/(1.+exp((r-wnRAu)/wnXi));}
	inline double wnRho(double x, double y, double z){return wnRhoR(sqrt(x*x+y*y+z*z));}
	inline double wnRhoIntegrand(double r){return r*r*wnRhoR(r);}
	inline double wnTIntegrand(double z){return wnRho(wnTx,wnTy,z);}
	double wnT(double x,double y); 
	void wnNormalize(); 
	
public:
	CCell **** mCells;
	static void setFOTemp(double TFO);
		
	CMesh() {}
	CMesh(parameterMap *pM);
	virtual ~CMesh();
	
	void initNS();
	void addInitialFlow();
	bool detectCrash();
	
		// integrate this mesh forward 
	void forward(CMesh*, double);
	
	void forward(CMesh*, CMesh*, double);
	void forward(CMesh*, CMesh*, CMesh*, CMesh*, double);
	
#ifdef HYDRO_BOOST_THREADS
	void tForward(CMesh* mMesh, int iThread);
	void tForward2(CMesh* mM1, CMesh* mM2, int iThread);
	void removeCell(list<CCell*>::iterator it);
#endif
		//puts average of the two meshs here
	void average(CMesh*,CMesh*);
	
	void smoothEdge(CCell*,CCell*,CCell*,CCell*);
	virtual void flipEdge() {/*noop*/}
	virtual void flipEdge(int,CCell*,CCell*,CCell*) {/*noop*/}
	
	inline void setParameterMap(parameterMap* pM)  {pMap = pM;}
	
#ifdef HYDRO_BOOST_THREADS
	void setEos(CEos* mEos);
	void setHelper(CCellHelper* mHelper);
#else
	inline void setEos(CEos* mEos) {getCell(0,0,0)->setEos(mEos);}
#endif
	inline void setEos(int eta, int x, int y, CEos* mEos) {getCell(eta,x,y)->setEos(mEos);}
	inline void setHelper(int eta, int x, int y, CCellHelper* mHelper) {getCell(eta,x,y)->setHelper(mHelper);}
	inline CEos* getEos(){return getCell(0,0,0)->getEos();}
	inline CCellHelper* getHelper(){return getCell(0,0,0)->getHelper();}
	
	inline void printCell(int eta,int x,int y) {return getCell(eta,x,y)->print();}
	inline void selfPrint(int eta,int x,int y) {return getCell(eta,x,y)->selfPrint();}
	
	void setTau(double mTau);
	inline double getTau(){return getCell(0,0,0)->getTau();}
	
	virtual CCell* getCell(int eta, int x, int y)=0;
	virtual void update(int eta, int x, int y)=0;
	virtual void selfUpdate(int eta, int x, int y)=0;
	inline bool getActive(int eta, int x, int y) {return getCell(eta,x,y)->getActive();}
	inline double getS(int eta, int x, int y, int s) {return getCell(eta,x,y)->getS(s);}
	inline double getU0(int eta, int x, int y) {return getCell(eta,x,y)->getU0();}
	inline double getDS(int eta, int x, int y, int m, int n) {return getCell(eta,x,y)->getDS(m,n);}
	inline double getE(int eta, int x, int y){return getCell(eta,x,y)->getE();}
    inline double getT(int eta, int x, int y){return getCell(eta,x,y)->getTCalc();}
    inline double getP(int eta, int x, int y){return getCell(eta,x,y)->getPCalc();}
	inline double getS(int eta, int x, int y){return getCell(eta,x,y)->getSCalc();}
	inline double getISAlpha(int eta, int x, int y){return getCell(eta,x,y)->getISAlphaCalc();}
	inline double getISGamma(int eta, int x, int y){return getCell(eta,x,y)->getISGammaCalc();}
	inline double getTIS(int eta, int x, int y){return getCell(eta,x,y)->getTISCalc();}
	inline double getTISB(int eta, int x, int y){return getCell(eta,x,y)->getTISBCalc();}
	inline double getX(int eta, int x, int y,int i){return getCell(eta,x,y)->getX(i);}
	inline double getTxy(int eta, int x, int y, int xx, int yy){return getCell(eta,x,y)->getTxyCalc(xx,yy);}
    inline double getPixy(int eta, int x, int y, int xx, int yy){return getCell(eta,x,y)->getPixyCalc(xx,yy);}
	inline double getPixyNS(int eta, int x, int y, int xx, int yy) {return getCell(eta,x,y)->getPixyNSLocal(xx,yy);}
	inline double getPixyNSMesh(int eta, int x, int y, int xx, int yy){return getCell(eta,x,y)->getPixyNS(xx,yy);}
	inline double getTxx(int eta, int x, int y) {return getTxy(eta,x,y,1,1);}
	inline double getTyy(int eta, int x, int y) {return getTxy(eta,x,y,2,2);}
	inline double getTzz(int eta, int x, int y) {return getTxy(eta,x,y,3,3);}
	inline double getSV (int eta, int x, int y) {return getCell(eta,x,y)->getSVCalc();}
	inline double getBV (int eta, int x, int y) {return getCell(eta,x,y)->getBVCalc();}
	
	inline double getX(int s) {return mCells[activeN][activeX][activeY]->getX(s);}
	inline double getS(int s) {return mCells[activeN][activeX][activeY]->getS(s);}
	inline double getDS(int m, int n) {return mCells[activeN][activeX][activeY]->getDS(m,n);}
	inline double getU0(){return mCells[activeN][activeX][activeY]->getU0();}
	inline double getE () {return mCells[activeN][activeX][activeY]->getE();}
	inline double getP () {return mCells[activeN][activeX][activeY]->getP();}
	inline double getS () {return mCells[activeN][activeX][activeY]->getS();}
	inline double getT () {return mCells[activeN][activeX][activeY]->getT();}
	inline double getSV () {return mCells[activeN][activeX][activeY]->getSV();}
	inline double getTIS () {return mCells[activeN][activeX][activeY]->getTIS();}
	inline double getISAlpha() {return mCells[activeN][activeX][activeY]->getISAlpha();}
	inline double getBV () {return mCells[activeN][activeX][activeY]->getBV();}
	inline double getTISB () {return mCells[activeN][activeX][activeY]->getTISB();}	
	inline double getISGamma() {return mCells[activeN][activeX][activeY]->getISGamma();}
	inline double getB () {return mCells[activeN][activeX][activeY]->getB();}
	inline void printActive() {mCells[activeN][activeX][activeY]->selfPrint();}
	inline double getTxy(int xx, int yy) {return mCells[activeN][activeX][activeY]->getTxy(xx,yy);}
	inline double getPixy(int xx, int yy) {return mCells[activeN][activeX][activeY]->getPixy(xx,yy);}
	inline double getPixyNS(int xx, int yy) {return mCells[activeN][activeX][activeY]->getPixyNS(xx,yy);}
	
	inline double getTZR() {return mCells[activeN][activeX][activeY]->getTZR();}
	inline double getTZPhi() {return mCells[activeN][activeX][activeY]->getTZPhi();}
	inline double getTRR() {return mCells[activeN][activeX][activeY]->getTRR();}
	inline double getTRPhi() {return mCells[activeN][activeX][activeY]->getTRPhi();}
	inline double getTPhiPhi() {return mCells[activeN][activeX][activeY]->getTPhiPhi();}
	
	inline double getPiZR() {return mCells[activeN][activeX][activeY]->getPiZR();}
	inline double getPiZPhi() {return mCells[activeN][activeX][activeY]->getPiZPhi();}
	inline double getPiRR() {return mCells[activeN][activeX][activeY]->getPiRR();}
	inline double getPiRPhi() {return mCells[activeN][activeX][activeY]->getPiRPhi();}
	inline double getPiPhiPhi() {return mCells[activeN][activeX][activeY]->getPiPhiPhi();}
	
	inline double getTZR(int eta, int x, int y){return getCell(eta,x,y)->getTZR();}
	inline double getTZPhi(int eta, int x, int y){return getCell(eta,x,y)->getTZPhi();};
	inline double getTRR(int eta, int x, int y){return getCell(eta,x,y)->getTRR();};
	inline double getTRPhi(int eta, int x, int y){return getCell(eta,x,y)->getTRPhi();};
	inline double getTPhiPhi(int eta, int x, int y){return getCell(eta,x,y)->getTPhiPhi();};
	
	inline double getPiZR(int eta, int x, int y){return getCell(eta,x,y)->getPiZR();};
	inline double getPiZPhi(int eta, int x, int y){return getCell(eta,x,y)->getPiZPhi();}
	inline double getPiRR(int eta, int x, int y){return getCell(eta,x,y)->getPiRR();}
	inline double getPiRPhi(int eta, int x, int y){return getCell(eta,x,y)->getPiRPhi();}
	inline double getPiPhiPhi(int eta, int x, int y){return getCell(eta,x,y)->getPiPhiPhi();}
	
	inline double getYL( int eta, int x, int y){return getCell(eta,x,y)->getYL();};
	inline double getVr( int eta, int x , int y){return getCell(eta,x,y)->getVr();} 
	inline double getGammaTrans( int eta, int x , int y){return getCell(eta,x,y)->getGammaTrans();}
	inline double getDUDT(int eta, int x, int y, int v){return getCell(eta,x,y)->getDUDT(v);}
	inline double getVorticity(int eta, int x, int y, int v, int w){return getCell(eta,x,y)->getVorticity(v,w);}
	
	inline double getYL() {return mCells[activeN][activeX][activeY]->getYL();}
	inline double getVr() {return mCells[activeN][activeX][activeY]->getVr();}
	inline double getGammaTrans() {return mCells[activeN][activeX][activeY]->getGammaTrans();}
	inline double getDUDT(int v) {return mCells[activeN][activeX][activeY]->getDUDT(v);}
	inline double getVorticity(int v, int w) {return mCells[activeN][activeX][activeY]->getVorticity(v,w); }
	
	double getDULocal(int eta, int x, int y, int u, int v);
	
	CCell* getNeighbor(int eta, int x, int y, int m, int n) {return getCell(eta,x,y)->getNeighbor(m,n);}
	
	inline int getNSize() {return mNSize;};
	inline int getXSize() {return mXSize;};
	inline int getYSize() {return mYSize;};
	inline int getNSizeOrig() {return mNSizeOrig;};
	inline int getXSizeOrig() {return mXSizeOrig;};
	inline int getYSizeOrig() {return mYSizeOrig;};
	
	inline void setNSize(int i) {mNSize = i;};
	inline void setXSize(int i) {mXSize = i;};
	inline void setYSize(int i) {mYSize = i;};
	
	void calcIntegrals();
	void calcLossIntegrals();
	
	virtual void writeIEPvEta()=0;
	
	inline double integralE() {return intE+eDead;};
	inline double integralS() {return intS+sDead;};
	inline double integralEp(){return intEp+epDead;};
	inline double integralEx(){return intEx+exDead;};
	inline double averageVr() {return intVr+vrDead;};
	
	virtual double getdNd3p(double px, double py, double pz, double m)=0;
	
	virtual void calcdNdY(double m, double d, double *output, double dRap, int nRap)=0;
	virtual double getF( int x, int y, int n, double pt, double phi, double rap, double m)=0;
	double getF(double beta, double dF=0.);
	
	virtual void getFOS(double mFOS[XSIZE+YSIZE][2], double mFOSigma[XSIZE+YSIZE][2], int &foSize)=0;
	virtual bool getFOS(double mFOS[XSIZE+YSIZE][2], double mFOSigma[XSIZE+YSIZE][2], 
						double mFOVelo[XSIZE+YSIZE][2], double mFODiss[XSIZE+YSIZE][6], int &foSize)=0; 
	virtual bool getFOS(int eta, double mFOS[XSIZE+YSIZE][2], double mFOSigma[XSIZE+YSIZE][2], 
						double mFOVelo[XSIZE+YSIZE][2], double mFODiss[XSIZE+YSIZE][6], int &foSize)=0;

	virtual bool goodFOS(double mFOS[XSIZE+YSIZE][2], int foSize)=0;
	
	double getFOSX();
	double getFOSY();
	double getFOSSST();
	double getFOSV();
	
	virtual void copyMesh(CMesh* m)=0;
	virtual int genFOS(CMesh* m)=0;
	bool containsSurf(double &fosE, double localE []);
		//	bool containsSurf(double &fosE, double localE[2][2][2][2]);
	bool containsSurf(double &fosE, double**** localE);
	
	inline double getELoss()  {return eLoss;}
	inline double getSLoss()  {return sLoss;}
	
		// remove the cells 
	virtual void deaden()=0;
	virtual void fillActiveCells()=0;
	
	inline void setActive(int eta, int x, int y, bool v){return getCell(eta,x,y)->setActive(v);}
	inline CCell* deactivate(int eta, int x, int y) {return getCell(eta,x,y)->deactivate();}
	
	bool anyActiveCells();
	void copyActive(CMesh* mMesh);
	void copyActive(CMesh* mM1, CMesh* mM2);
	void cleanActiveCells();
	
	virtual void checkAzimuthalSymmetry() {}
	virtual void fixEntropyMesh()=0;
	
	void cubeInterp(double mS[11], double midX[], CCell* cube[], double mDt);
	void hcubeInterp(double mS[11], double foPoint[4], CCell* cube[2][2][2][2], double mDt);
	
	int cubical2p1(double e0, double eCubic[8], double surfvec[3], 
				   int* nsurfptr, double Vmid[3], double dt, double dx, 
				   double dy, int* Nerrptr);
	
	inline void setOsuFos(FILE *f){osuFOS=f;}
	inline void closeFOS(){fclose(osuFOS);}
};

#endif // __CMESH_h_INCLUDE__