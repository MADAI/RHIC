#ifndef __PCA_H__
#define __PCA_H__
#include "coralutils.h"
#include "qualifier.h"
using namespace std;

class CPCA{
public:
	string yname[1000];
	double sigmay[1000];
	int nruns,ny,nnames;
	CQualifiers qualifiers;
	CPCA(int nruns_set);
	double *ybar,*value,**spread;
	void ReadResults();
	double *eigenval;
	double **evec;
	CGSLMatrix_Real *gslmatrix;
	bool namecheck(string varname);
	string pcaname[10000];
	void Calc();
};

class CRHICStat{
	string yname[1000];
	string xname[1000];
	int NX,NY,NRUNS;
	double *xbar,*ybar,**x,**y,**xprime,**yprime;
	CRHICStat(int NRUNS);
	void CalcSensitivity();
};

#endif

