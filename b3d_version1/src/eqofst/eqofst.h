#ifndef __INCLUDE_EQOFST_H__
#define __INCLUDE_EQOFST_H__

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "intrinsic.h"
#include "species.h"
#include "coral.h"
#include "coralutils.h"
#include "constants.h"
#include "coral.h"

using namespace std;

class CEqofst{
public:
	CEqofst();
	CEqofst(CSpecies *species);
	~CEqofst();
string parsfilename;
	CSpecies *species;
	void PECheck(CIntrinsic *intrinsic);
	//
	virtual void Calc_of_TMu(CIntrinsic *intrinsic);
	virtual void Calc_of_EMu(CIntrinsic *intrinsic);
	virtual void Calc_of_SAlpha(CIntrinsic *intrinsic);
	virtual void Calc_of_SRho(CIntrinsic *intrinsic);
	virtual void Calc_of_ERho(CIntrinsic *intrinsic);
	virtual void Calc_of_ERhoOverS(CIntrinsic *intrinsic);
	//
	virtual void CalcVisc(CIntrinsic *intrinsic);
	//
	virtual void InitIntrinsic_EMu0(CIntrinsic *&intrinsic,double epsilon);
	//
	int MaxNcharges;
	double etafact;
	double Bfact;

	void freegascalc_onespecies(double m,double t,double &p,double &e,double &dens,double &sigma2,double &dedt);
	void FreeGasCalc_of_TMu(CIntrinsic *intrinsic);
	double GetLambdaFact(CIntrinsic *intrinsic);
	double GetBfact(double m,double T); // \int d^3p p^4/e^3 exp(-e/T)
	double GetSigmaB(CIntrinsic *intrinsic);
	double GetSigmaB_test(CIntrinsic *intrinsic);

protected:   // these need to use newton's methods and related arrays need to created first 
	void FreeGasCalc_of_EMu(CIntrinsic *intrinsic);
	void FreeGasCalc_of_SAlpha(CIntrinsic *intrinsic);
	void FreeGasCalc_of_SRho(CIntrinsic *intrinsic);
	void FreeGasCalc_of_ERho(CIntrinsic *intrinsic);
	void FreeGasCalc_of_ERhoOverS(CIntrinsic *intrinsic);
	//
	//
	double *y,*x,*rhotry,**M;
	CGSLMatrix_Real **gslmatrix;
	parameterMap eqofstpars;
	void MakeNewtonArrays();
	void DestroyNewtonArrays();
	void FindXwhereYeqMX(int NDIM);
};

class CEqofst_ThreePhase : public CEqofst{
public:
	CEqofst_ThreePhase(CSpecies *hspecies,string parsfilename_set);

	void Calc_of_TMu(CIntrinsic *intrinsic);
	void Calc_of_EMu(CIntrinsic *intrinsic);
	void Calc_of_ERho(CIntrinsic *intrinsic);
	//
	void CalcVisc(CIntrinsic *intrinsic);
	//
	void InitIntrinsic_EMu0(CIntrinsic *&intrinsic,double epsilon);
	//
	CIntrinsic *intrinsicTc;
	CSpecies *hspecies,*qgpspecies;
	//
	double L,T_h,h_h,P_h,epsilon_h,sdens_h,eta_h,sigma2_h;
	double T_qgp,epsilon_qgp,P_qgp,h_qgp,sdens_qgp,eta_qgp,sigma2_qgp,c2qgp,c2mixed;

};


#endif
