#ifndef __CCELL_h_INCLUDE__
#define __CCELL_h_INCLUDE__
/*
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

#include "hydroDef.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <coralutils.h>
#ifdef HYDRO_BOOST_THREADS
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#endif
#include "CEos.h"

	// Container for variables that we'll calculate
	// every time we update() a cell
	// Usually CCell will have one static pointer to a CCH
	// but if multithreaded each cell has a pointer
class CCellHelper {
private:
		// Viscous Shear Tensor
	double Pixy[4][4];	
	
		// Dynamic Bulk Viscosity
	double BPi;
	
		// derivatives of s_i
	double dS[11][3];	
	
		// spatial portion of derivatives including affine connections
	double dU[4][4];
	double DPixy[4][4][4];
	
		// local matrix for calculating the time derivative
	double M[11][12];		
		// derivative of the Israel-Stewart Scaling factor
	double dTISS[3];	
		// needed when HYDRO_ISMAX == true
	double a, y;
	double gamma;
	
		// trace elements of metric
	double g[4];
	
public:
	CCellHelper();
	CCellHelper(CCellHelper& h);
	CCellHelper(CCellHelper* p);
	~CCellHelper();
	
	friend class CCell;
};

class CCell;

struct CCellH5 {
	double x[4];
	double e, P, T, fQGP;
	double v[3];
	double a[5];
	double b, SV, tIS, BV, tISB, S;
};

class CCell { 
	
private:
		// local variables
	
		// cell location
	double x[4];					
		// pointers to cell's neighbors
		// not all neighbors will be active,
		// but mesh always has dead cells at the boundary
	CCell* neighbors[3][2];		
		// set of hydrodynamic quantities for this cell
		// [e; u_x, u_y, u_z; Pi_xx, Pi_yy, Pi_zz, Pi_xy, Pi_xz, Pi_yz]
	double s[11];				
		// time derivative of u
	double dUdT[3];
		// true if cell is active
	bool active;
	
		//working class pointers
#ifdef HYDRO_BOOST_THREADS
	CEos* eos;
	CCellHelper* helper;
#else
	static CEos* eos;
	static CCellHelper* helper;	
#endif
	
	static parameterMap* pMap;
	
		// options and value set by the parameterMap
	static bool mDebug, mSVTrim, mSVTrimInit, mViscNS, mPureBjorken, mBjorken, mSlopeLimit, mLaxFried;
	static bool mLinT, mLogT, mLogSinhT, mISVort, mISMax;
	static int mOctant;
	static double mT0, mSVRatio, mBVRatio, mISAMax, mISBMax, mInitFlow, mInitNS, mInitNSTxxP, mSlopeLimitTheta;
	
protected:
		// functions to fill calc variables using internal information
	
		// fill the augmented coefficient matrix
	void   fillM();
	
		// returns entries to the augmented coefficient matrix
	double getECoeff(int j);
	double getPCoeff(int i, int j);
	double getPiCoeff(int i, int j);
	double getPiCoeff(int j);
	
		// recalculates derivatives
	void   calcDeriv();				
	
		// make all the eos calls for static variables
	void   fillEosVar();
	
		// update helper function
	void fillPixy();
	
		// calcDeriv() helper functions
	void fillDS();
	void fillDU();
	void fillDPixy();
	
public:
		// differential elements of position
	static double dx[4];
	
		// construct dummy cell in order to set the parameterMap
		// probably unnecessary
	CCell(parameterMap* pMap);
	
		// creates cell with 4 coordinates
	CCell(double,double,double,double);
	
		// copy constructor
	CCell(CCell* cell);
	
		// destructor
	~CCell();
	
		// Fill parameters from parameterMap
		// Called by dummy constructor
	void paramFill();
	
		// grabs:
	
	inline bool getActive() {return active;}
	
		// cell position
	double getTau();
	inline double getEta() {return x[3];}
	inline double getZ() {return getTau()*sinh(x[3]);}
	inline double getTime() {return getTau()*cosh(x[3]);}
	inline double getX() {return x[1];}
	inline double getY() {return x[2];}
	inline double getX(int i) {if (i==0) return getTau(); else if(i<4 && i>0) return x[i]; else return 0.;}
	inline double getDx(int i) {if (i<4 && i>=0) return dx[i]; else return 0.;}
	
		//relativistic velocities 
	inline double getU0() {return sqrt(1.+s[1]*s[1]+s[2]*s[2]-helper->g[3]*s[3]*s[3]);}
	inline double getUx() {return s[1];}
	inline double getUy() {return s[2];}
	inline double getUz() {return s[3];}
	inline double getU(int i) {if (i>0) return s[i]; else return getU0();}
	
	inline double getVx() {return s[1]/helper->gamma;}
	inline double getVy() {return s[2]/helper->gamma;}
	inline double getVz() {return s[3]*getTau()/helper->gamma;}
	
		// radial velocity
	inline double getVr() { return sqrt(s[1]*s[1]+s[2]*s[2])/getU0();}
		// lorentz factor due to radial velocity only
	inline double getGammaTrans(){ return sqrt( 1. + s[1]*s[1] + s[2]*s[2]);}
		// longitudinal velocity in the lab frame
	inline double getUzRest() { return sinh(x[3])*(1.+ (s[3]*s[3]/(s[0]+1.))) + cosh(x[3])*s[3];}
		// lorentz factor in lab frame
	inline double getGammaRest() {return s[0]*cosh(x[3])+s[3]*sinh(x[3]);}
		// lorentz factor in lab frame
	inline double getVzRest() { return getUzRest()/getGammaRest();}
		// longitudinal rapidity
	inline double getYL() { return (0.5) * log((1+getVzRest())/(1-getVzRest()));}
	inline double getDUDT(int v) { if (v>0) return dUdT[v-1]; else return -(dUdT[0]*s[1]+dUdT[1]*s[2]+dUdT[2]*s[3])/s[0];}
	double getExpansionRate();
	
		// longitudinal velocity of the mesh at this {eta,tau}
	inline double getUMesh() {return sqrt(1+pow(tanh(x[3]),2))*tanh(x[3]);}
	
		//spatial components of stress energy tensor
	inline double getA(int v) {return s[4+v];}
	inline double getB() {return s[10];}
		// energy density
	inline double getE() {return exp(s[0]);}
	
		// returns spatial part of tensor in fluid frame
	double getTxy(int,int);						
	double getTxy(CCell*,int,int);
	double getPixy(int,int);
	double getPixyBulk(int,int);
	double getTxyLocal(int,int);
	double getPixyLocal(int,int);
	double getTxyCalc(int,int);
	double getPixyCalc(int,int);
	double getTxyLocalCalc(int,int);
	double getPixyLocalCalc(int,int);
	double getPixyNS(int,int);
	double getPixyNSLocal(int,int);
	
	inline double getBPiCalc(){return getISAlphaCalc()*s[10];}
	
		// direct eosCalc calls
	inline double getSCalc() {return eos->getS(getE());}
	inline double getPCalc() {return eos->getP(getE());}
	inline double getCs2Calc() {return eos->getCs2(getE());}
	inline double getSEquilCalc() {return eos->getS(getE());}
	inline double getTISCalc() {if (!mSVTrim) return eos->getTIS(getE()); else return eos->getTIS(getE())/(1.+exp((x[1]*x[1]+x[2]*x[2]-10.)/0.6));}
	inline double getTISBCalc() {return eos->getTISB(getE());}
	inline double getTCalc() {return eos->getT(getE());}
	inline double getSigmaACalc() {return eos->getSigmaA(getE());}  
	inline double getSigmaBCalc() {return eos->getSigmaB(getE());}  
	inline double getISAlphaCalc(){return eos->getISAlpha(getE());}
	inline double getISGammaCalc(){return eos->getISGamma(getE());}
	inline double getISBetaCalc(){return eos->getISBeta(getE());}
	inline double getDISAlphaDECalc(){return eos->getDISAlphaDE(getE());}
	inline double getDISGammaDECalc(){return eos->getDISGammaDE(getE());}
	inline double getSVCalc() {if (!mSVTrim) return eos->getSV(getE()); else return eos->getSV(getE())/(1.+exp((x[1]*x[1]+x[2]*x[2]-10.)/0.6));}
	inline double getBVCalc() {return eos->getBV(getE());}
	inline double getAMaxCalc() { return eos->getAMax(getE());}
	inline double getBMaxCalc() { return eos->getBMax(getE());}
	inline double getDAMaxDECalc() { return eos->getDAMaxDE(getE());}
	inline double getDBMaxDECalc() { return eos->getDBMaxDE(getE());}
	
		// avoids eos call
	inline double getS() {return eos->getS();}
	inline double getP() {return eos->getP();}
	inline double getCs2() {return eos->getCs2();}
	inline double getSEquil() {return eos->getS();}
	inline double getTIS() {return eos->getTIS();}
	inline double getTISB() {return eos->getTISB();}
	inline double getT() {return eos->getT();}
	inline double getSigmaA() {return eos->getSigmaA();}  
	inline double getSigmaB() {return eos->getSigmaB();}  
	inline double getISAlpha(){return eos->getISAlpha();}
	inline double getISGamma(){return eos->getISGamma();}
	inline double getISBeta(){return eos->getISBeta();}
	inline double getDISAlphaDE(){return eos->getDISAlphaDE();}
	inline double getDISGammaDE(){return eos->getDISGammaDE();}
	inline double getSV() {return eos->getSV();}
	inline double getBV() {return eos->getBV();}
	inline double getAMax() { return eos->getAMax();}
	inline double getBMax() { return eos->getBMax();}
	inline double getDAMaxDE() { return eos->getDAMaxDE();}
	inline double getDBMaxDE() { return eos->getDBMaxDE();}
	
		// return elements of s
	inline double getS(int i) {return s[i];}
		// return mesh derivatives of quantities is s
	inline double getDS(int i,int j) { return helper->dS[i][j];}
	
		// Calculates the 4x4 tensor of local derivatives of velocity
		// using the mesh derivatives information in dS[][]
		// *Does not include effects of time derivatives*
	//	double getSpLocDS(int,int);
	//	double getDiPixyLocal(int,int,int);
	
		// returns the vorticity in the local frame 0.5*(d_i v_j - d_j v_i)
		// Uses dULocal so cell must be updated first
	double getVorticity(int i,int j);
	
		// return working class
	inline CEos* getEos() {return eos;} 
	inline CCellHelper* getHelper() {return helper;}
	
	inline CCell* getNeighbor(int m, int n) {return neighbors[m][n];}
	
		// Variable sets:
	
		// cell positions
	inline void setTau(double v) {if (mLinT) x[0]=v; else if (mLogT) x[0] = log(v/mT0); else if (mLogSinhT) x[0] = log(sinh(v));}
	inline void setEta(double v) {x[3] = v;}
	inline void setX(double v)   {x[1] = v;}
	inline void setY(double v)   {x[2] = v;}
	
	inline void setX(double v[4]) {setTau(v[0]); x[1]=v[1]; x[2]=v[2]; x[3]=v[3];}
	
	inline void setDt(double v) {dx[0] = v;}
	inline void setDx(int i, double v) {dx[i]=v;}
	
	inline void setActive(bool v) {	active = v;}
	
		// relativistic velocities  (gamma v)
	void setU(double,double,double);
	inline void setUx(double v) {s[1] = v;}
	inline void setUy(double v) {s[2] = v;}
	inline void setUz(double v) {s[3] = v;}
	
		// spatial components of stress energy tensor
	inline void setA(int v, double v1) {if (v>0 && v<6) s[3+v] = v1;}
	inline void setB(double v) {s[10] = v;}
	inline void setE(double v) {s[0] = log(v);}
	void setS(int i, double v);
	inline void setS(double v[11]){for (int i=0;i<11;i++) s[i]=v[i];}
	
	inline void setDUDT(int i, double v) {dUdT[i]=v;}
	
		//working class
	inline void setEos(CEos* p) {eos = p;}
	inline void setHelper(CCellHelper* p ) {helper = p;}
	
		// set neighbors
	inline void setEtaNeighbors(CCell* p1, CCell* p2) {neighbors[2][0]=p1; neighbors[2][1]=p2;}
	inline void setXNeighbors(CCell* p1, CCell* p2) {neighbors[0][0]=p1; neighbors[0][1]=p2;}
	inline void setYNeighbors(CCell* p1, CCell* p2) {neighbors[1][0]=p1; neighbors[1][1]=p2;}
	
		// updates: up[i], uup
		// calls: calcDeriv(), FillM()
		// called by: forward()
	void update();

		// sets shear elements to the Navier-Stokes Value
		// calls: update.
	void initNS();
	
	void addInitialFlow();
	
		// turns this cell off
		// if adjacent cell is along axis, turns that cell off
	CCell* deactivate();
	
		// integrates the cell forward in time
		// uses derivatives from this cell
		// puts result into passed cell
	void forward(CCell*);
	
		// integrates the cell forward in time
		// uses derivatives from this cell
		// adds result to 1st passed cell and assigns them to 2nd passed cell
		// Essential to Runge-Kutta for 2nd order or higher accuracy in time
	void forward(CCell*, CCell*);
	
		// Outline for the final step in RK4 integration
		// Not fully implemented.
	void forward(CCell*, CCell*, CCell*, CCell*);
	
	void copy(CCell* mC);
	void copyST(CCell* mC);
	void flip(int i);
	
		// prints
	void print();
	void selfPrint();
	
	void selfUpdate();
	void verify();
	
	double getTZR();
	double getTZPhi();
	double getTRR();
	double getTRPhi();
	double getTPhiPhi();
	double getPiZR();
	double getPiZPhi();
	double getPiRR();
	double getPiRPhi();
	double getPiPhiPhi();

	// state of dt coefficient matrix
	void printM();
	void printMm(double Mm[11][12]);
	void getM(double mM[11][12]);
	
	double minmod(double d1, double d2, double d3);
	double max(double d1, double d2, double d3);
	double min(double d1, double d2, double d3);
	
	double max(double d1, double d2);
	double min(double d1, double d2);
	double minmod(double d1, double d2);
	
	friend bool operator!= (CCell &c1, CCell &c2);
};
	//#include "CCell.cpp"
#endif //__CCELL_h_INCLUDE__
