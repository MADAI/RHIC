#ifndef __RHICSTAT_H__
#define __RHICSTAT_H__
#include "coralutils.h"
#include "qualifier.h"
using namespace std;

class CRHICStat;

class CRunInfo{  /** this is info that pertains to a specific run */
public:
	CRunInfo(int NX,int NY);
	double *x,*w;
	double *y,*z,*zfiterror;
	double *sigmay;
	double *ylinear,*zlinear,*xlinear,*zfit,*yfit,*xquad;
	double netdiff_exp,netdiff_fit,netdiff_fitexp;
	bool good;
	void Print();
	static CRHICStat *rhicstat;
};

class CZGetter{
public:
	CZGetter();
	CRHICStat *rhicstat;
	virtual void GetZ(double *x,double *z);
	virtual void GetZError(double *x,double *zfiterror);
	bool GETERROR;
};

class CZGetter_QuadFit : public CZGetter {
public:
	CZGetter_QuadFit(CRHICStat *rsptr);
	void GetZ(double *x,double *z);
	void PrintQuadCode();
	void InitQuadFit();
	double ***Aquad,**Bquad,*Cquad;
	void PrintCoefficients();
};

class CZGetter_GP : public CZGetter {
public:
	CZGetter_GP(CRHICStat *rsptr);
	void InitInterpolator();
	void LinearFit();
	void GetZ(double *x,double *z);
	void GetZError(double *x,double *error);
	double GetCov(int iz,CRunInfo *runinfo1,CRunInfo *runinfo2);
	double GetCov(int iz,double *x1,double *x2);
	void CalcHyperPars();
	void CalcCovInv();
	void PrintHyperPars();
	void ReadHyperPars();
	void GetCovDerivs(int iz,double *x1,double *x2,double &C,double *DC,double **DDC);
	double ***Cov,***CovInv,**CovInvDotZ,**alphanorm;
	double **hyperR,*hyperTheta0,*hyperNugget;
	double *Dzsquared;
	double **mlinear,*blinear;
	bool LINEAROFF;
	double hyperPower;
	bool READ_HYPERPARS,NORMALIZE;
	double hyperR_default,hyperNugget_default;
	CGSLMatrix_Real *gslmatrix_NRUNS;
	bool COVINVCALC;
	int NX,NZ,NRUNS,NTESTRUNS; //copied from rhicstat
};

class CZGetter_LocalLinear : public CZGetter {
public:
	CZGetter_LocalLinear(CRHICStat *rpstr);
	void InitInterpolator();
	void InitLocalLinearFit();
	CGSLMatrix_Real *gslmatrix;
	void GetZ(double *x,double *z);
	double *ll_slope,*ll_xz,**ll_xx,**ll_xxinv,*xbar,*weight;
	double ll_delx;
};

class CRHICStat{
public:
	CRHICStat(string statpars_filename);
	parameterMap parmap;
	string FIT_TYPE;
	bool USE_MODEL_ERROR,USE_COSH_LL;
	string PATH_NAME,MODEL_OUTPUT_DIR,PARAMETERS_DIR;
	long long int NMCMC,NBURN;
	CZGetter *zgetter;
	string *yname;
	string *xname;
	double *xmin,*xmax;
	int NX,NY,NZ,NRUNS,NTESTRUNS,NGOODRUNS;
	double SIGMA2_EMULATOR; // squared error of emulator for single z component
	double *xbar,*ybar,*sigmaybar,*xmcmc,**xhatmcmc;
	double **sigmaxx,**sigmayy,**dxdz,**dxdz_inv;
	double **Uytoz,**Uytoz_inv,**Uxtow,**Uxtow_inv;
	double *eigenvalyy,*eigenvalxx;
	double *uncertainty;
	CRunInfo **runinfo,**testinfo;
	CRunInfo *expinfo,*fitinfo,*bestinfo;
	CGSLMatrix_Real *gslmatrix_NY;
	CGSLMatrix_Real *gslmatrix_NX;
	CRandom *randy;
	double MCMC_STEPSIZE;
	bool PCA_CALC_DXDZ;

	void FitExpData();
	void ScaleXY();
	void InitArrays();
	void InitX();
	void InitY();
	void ReadAllX();
	void ReadX(string filename,CRunInfo *runinfo);
	void ReadAllY();
	void ReadY(string filename,CRunInfo *runinfo);
	void PCA();
	void GetZFromY(CRunInfo *runinfo);
	void GetYFromZ(CRunInfo *runinfo);
	void GetXlinearFromZ(CRunInfo *runinfo);
	void GetZlinearFromX(CRunInfo *runinfo);
	void GetZlinearFromXlinear(CRunInfo *runinfo);
	void GetXlinearFromY(CRunInfo *runinfo);
	void GetYlinearFromX(CRunInfo *runinfo);
	void GetXlinearFromZlinear(CRunInfo *runinfo);
	void GetYlinearFromZlinear(CRunInfo *runinfo);
	void GetYfitFromZfit(CRunInfo *runinfo);
	void GetWFromXlinear(CRunInfo *runinfo);
	void GetXlinearFromW(CRunInfo *runinfo);
	void PrintXlinear(CRunInfo *runinfo);
	void PrintX(CRunInfo *runinfo);
	void PrintY(CRunInfo *runinfo);
	void PrintZ(CRunInfo *runinfo);
	void PrintZlinear(CRunInfo *runinfo);
	void PrintYlinear(CRunInfo *runinfo);
	void PrintZfit(CRunInfo *runinfo);
	void PrintYfit(CRunInfo *runinfo);
	void WritePars(string filename);
	void CheckTestRuns();
	void PlotZvsX();
	double GetLL(double *x);
	void Metropolis();
	void MetropolisOpt();
	void Metropolis(long long int nburn,long long int nmcmc,double stepsize);
	void MetropolisOpt(long long int nburn,long long int nmcmc,double stepsize);
	void CalcError(CRunInfo*);
	void CalcNetDiffExp(CRunInfo *runinfo);
	void CalcNetDiffFit(CRunInfo *runinfo);
	void CalcNetDiffFitExp(CRunInfo *runinfo);
	void CalcCovariance();
	void PerformFits();
};



#endif

