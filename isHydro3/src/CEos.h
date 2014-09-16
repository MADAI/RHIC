#ifndef __CEOS_h_INCLUDE__
#define __CEOS_h_INCLUDE__

/*
 *  CEos.h
 *  Created by Joshua Vredevoogd on 3/4/09.
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sys/stat.h> 
#include <coralutils.h>
#include <cassert>
#include "hydroDef.h"

	// calculates EOS variables
class CEosCalculator{
  protected:
		// equation of state variables
    static double *temp, *ed, *pr, *sd;
		// second derivatives for splines
	static double *t2, *e2, *p2, *s2;
	static int aSize;
	
		// As a test, I've made everything but the searching variable static.
		// This is so that one can make multiple, simultaneous
		// searches of the EoS, as in from different threads.
	int lastAccess;

	static bool mEosUseSpline,mEosJosh,mEosScott;
	static double mSVRatio, mBVRatio, eAtT, sAtT;
	static double svSwitchTemp, mFoTemp, mSvHighT, mSvHighTSlope;
	static double mBVCent, mBVWidth;
	static bool bExpTail, bPowerTail, bGaussTail, bRealisticEtaS, bBagMerge;

		// functions that add entropy production in a temperature region
		// consistent with keeping 0 < cs2 < 1/3
	void addBump(double T0, double width, double height);
	void getMaxAmp(double T0, double width, double &minH, double &maxH );
	void getMaxAmp(float T0, float mWidth, float &minH, float &maxH);
	void getMaxAmp(float T0, float &minH, float &maxH);
	bool bounds (float T0, float mWidth, float height);
	
		// Fills lastAccess with the array index just lower than x
		// Returns -1 for below ed[0], returns aSize-1 for above ed[aSize-1]
	void search(double e);
	
		// Compute the second derivatives via tridiagonalization
		// Called once during construction
		// Required for splint() interpolation 
		// Algorithm from Numerical Recipes for C++ Press et al. (Section 3.3)
		// Inputs are (1,2) coordinates of the data points to be fit
		// (3) second derivative at each data point
		// (4,5) first derivative at the first and last data point
	void spline(double *x, double *y, double *y2, double yp1, double ypA);
	
		// Compute the interpolated value using spline (from NR)
		// Algorithm from Numerical Recipes for C++, Press et al. (Section 3.3)
	double splint(double *y, double *y2, double e);
	
		// Instead of value, returns derivative (dydx) at x=e
		// Useful for interpolating speed of sound squared
		// Algorithm from Numerical Recipes for C++, Press et al. (Section 3.3)
	double splintDeriv(double *y, double *y2, double e);
	
		// Calculates derivative dydx(x=x1) using quadratic ansatz
	double getDeriv(double x1, double x2, double x3, double y1, double y2, double y3);
	
		// verifies a file's existence
		// swiped wholesale from: http://www.techbytes.ca/techbyte103.html
	bool checkFile(string fName);
	
	void InitializeEosJosh(parameterMap* pM);
	void InitializeEosScott(parameterMap* pM);
	
public:
	CEosCalculator(parameterMap* pMap);
	~CEosCalculator();
	CEosCalculator();
	CEosCalculator(CEosCalculator& p);
	CEosCalculator(CEosCalculator* p);
	
	double getCs2(double e);
	double getP(double e);
	double getS(double e);
	double getTIS(double e);
	double getTISB(double e);
	double getSV(double e);
	double getBV(double e);
	double getTA(double e);
	double getT(double e);
	double getSigmaA(double e);
	double getSigmaB(double e);
	double getISAlpha(double e);
	double getISBeta(double e);
	double getISGamma(double e);
	double getDISAlphaDE(double e);
	double getDISGammaDE(double e);
	double getAMax(double e);
	double getBMax(double e);
	double getDAMaxDE(double e);
	double getDBMaxDE(double e);
	
	double getEGivenT(double T);
	double getSGivenT(double T);
	double getEGivenS(double mS);
	
		// dump the values 
	void printEos(const char*);
	void printEosEnergyUnits(const char*);
	void printEos(FILE*);
	void printEosEnergyUnits(FILE*);
	
};

	// container for variables in the equation of state
class CEos {
	
protected:
		// EOS variables - avoids unnecessary EOS calls
	double cS2, P, S, tIS, tISB, SV, BV, T;
	double sigmaA, sigmaB, alphaIS, gammaIS, betaIS, aMax, bMax;
		// derivatives of thermodynamic quantities
	double dAlphaISDE, dGammaISDE, dAMaxDE, dBMaxDE;
		
	static bool mISMax;
	CEosCalculator* mEosC;
	
public:
	CEos();
	CEos(parameterMap* p);
	~CEos();
		// full copy
	CEos(CEos& e);
		// generates a new CEosCalculator
	CEos(CEos* p);
	
	inline double getCs2(double e) {return mEosC->getCs2(e);};
	inline double getP(double e) {return mEosC->getP(e);};
	inline double getS(double e) {return mEosC->getS(e);};
	inline double getTIS(double e) {return mEosC->getTIS(e);};
	inline double getTISB(double e) {return mEosC->getTISB(e);};
	inline double getSV(double e) {return mEosC->getSV(e);};
	inline double getBV(double e) {return mEosC->getBV(e);};
	inline double getTA(double e) {return mEosC->getTA(e);};
	inline double getT(double e) {return mEosC->getT(e);};
	inline double getSigmaA(double e) {return mEosC->getSigmaA(e);};
	inline double getSigmaB(double e) {return mEosC->getSigmaB(e);};
	inline double getISAlpha(double e) {return mEosC->getISAlpha(e);};
	inline double getISBeta(double e) {return mEosC->getISBeta(e);};
	inline double getISGamma(double e) {return mEosC->getISGamma(e);};
	inline double getDISAlphaDE(double e) {return mEosC->getDISAlphaDE(e);};
	inline double getDISGammaDE(double e) {return mEosC->getDISGammaDE(e);};
	inline double getAMax(double e) {return mEosC->getAMax(e);};
	inline double getBMax(double e) {return mEosC->getBMax(e);};
	inline double getDAMaxDE(double e) {return mEosC->getDAMaxDE(e);};
	inline double getDBMaxDE(double e) {return mEosC->getDBMaxDE(e);};
	
	inline double getEGivenT(double e) {return mEosC->getEGivenT(e);};
	inline double getSGivenT(double e) {return mEosC->getSGivenT(e);};
	inline double getEGivenS(double mS) {return mEosC->getEGivenS(mS);};
	
		// dump the values 
	inline void printEos(const char* v) {mEosC->printEos(v);};
	inline void printEosEnergyUnits(const char* v) {mEosC->printEosEnergyUnits(v);};
	inline void printEos(FILE* f) {mEosC->printEos(f);};
	inline void printEosEnergyUnits(FILE* f) {mEosC->printEosEnergyUnits(f);};
												
		// Use CEosCalculator to fill CEos variables
	void fill(double e);
	
		// for cleanliness cells ask CEos to adjust
		// if the shear viscosity is position dependent
		// NOTE :: changes here need to migrate to:
		// CCell::getSVCalc(double) and CCell::getTISCalc(double)
	void trimSV(double r);
	
	inline double getS() {return S;};
	inline double getP() {return P;}
	inline double getCs2() {return cS2;}
	inline double getTIS() {return tIS;}
	inline double getTISB() {return tISB;}
	inline double getT() {return T;}
	inline double getSigmaA() {return sigmaA;}  
	inline double getSigmaB() {return sigmaB;}  
	inline double getISAlpha(){return alphaIS;}
	inline double getISGamma(){return gammaIS;}
	inline double getISBeta(){return betaIS;}
	inline double getISAMax() {return aMax;}
	inline double getISBMax() {return bMax;}
	inline double getDISAlphaDE(){return dAlphaISDE;}
	inline double getDISGammaDE(){return dGammaISDE;}
	inline double getSV(){return SV;}
	inline double getBV(){return BV;}
	inline double getAMax() { return aMax;}
	inline double getBMax() { return bMax;}
	inline double getDAMaxDE() { return dAMaxDE;}
	inline double getDBMaxDE() { return dBMaxDE;}
};

#endif //__CEOS_h_INCLUDE__
