#ifndef __JOSHCONVERT_H__
#define __JOSHCONVERT_H__
#include <cstdio>
#include <vector>
#include <cstdlib>
#include "b3d.h"

using namespace std;

class CRing;

class CjoshConvert{
public:
	parameterMap parmap;
	string run_name,qualifier;
	int NRINGSMAX;
	CResList *reslist;
	double Tf,ETAMAX;
	// will generate particles with -eta_max < eta < eta_max
	double epsilonf,Pf,nhadronsf;
	vector<double> densityf,boseweight;
	void WriteNewFormat();
	double WriteRing(int iring);
	void Init(string run_name);
	void ReadInput(string qualifier);
	void freegascalc_onespecies(double m,double t,double &p,double &e,double &dens,double &sigma2,double &dedt);
	CRandom *randy;
protected:
	int nres,nrings,nprcells,NBOSE;
	FILE *input,*vertex_output,*triangle_output;
	vector<CRing> ring;
	vector<vector<CRing> > ring3d;
	vector<int> nrings3d;
	void GetRingInfo();
	void ReadHeader(FILE *);
};

// This fills out Pi[4][4] given Pi_xx, Pi_yy, Pi_xy, u_x and u_y (all specified in lab frame)
void FillOutPi(double **Pi,double pixx,double piyy,double pixy,double ux,double uy);

// Stores info about a ring that breaks up at specific tau
class CRing{
public:
	CRing();
	void Print();
	void Read(FILE *fptr);
	void Read3D(FILE *fptr);
	void FillOutArrays(double *x,double *y,double *uxi,double *uyi,double *dToverH_xx,double *dToverH_xy,double *dToverH_yy,double *XX);
	const static int nphi=15;
	int nread;
	double tau,etamin,etamax;
	double r[nphi+1],ux[nphi+1],uy[nphi+1],T[nphi+1],Xscale[nphi+1];
	double dToverH[nphi+1][4][4];
	static CjoshConvert *jc;
};

#endif
