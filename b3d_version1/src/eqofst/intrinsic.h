#ifndef __INTRINSIC_INCLUDE_H__
#define __INTRINSIC_INCLUDE_H__

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "species.h"

// CLASS CIntrinsic ******************************************
class CIntrinsic{
public:
	CSpecies *species;
	double epsilon,T,P,dedT,sdens,eta,B,tauIS_a,tauIS_b,asigma,bsigma,amax,bmax;
	double *mu,*density,*rho,*flux;
	//
	void Print();
	CIntrinsic(CSpecies *species);
	~CIntrinsic();
	void Reset(CSpecies *species);
	void ZeroMu();
	void CopyFrom(CIntrinsic *intr);
	double GetLambdaFact();
};

#endif
