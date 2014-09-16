#ifndef __bjmaker_h__
#define __bjmaker_h__

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "b3d.h"

using namespace std;
class CRing; class CResList;

class CBjMaker{
public:
	//CHYDROtoB3D();
	// lambdafact is the ratio s.t., lambda_ij = pi_ij / lambdafact,
	// where pi_ij is the shear tensor in matter frame and lambda_ij tells how the momenta are scaled 
	// during collision time  p_i = ptilde_i + lambda_ij ptilde_j, where ptilde is generated isotropically
	double T,etamax,etagauss,Rx,Ry,tau;
	bool gaussian,balance;
	void Init();
	int GenerateParticles(); // returns nparts
	void GenerateParticles_Gaussian(int nparts); // 
	void GenerateParticles_Gaussian_Balance(int nparts); // makes pairs of pi+pi- or pi0pi0
	// will generate particles with -eta_max < eta < eta_max
	CResList *reslist;
	double epsilon,P,lambdafact;
	int MakeEvent();
	void Init(double T_in,double etamax_in,string inputfilename_in);
	CB3D *b3d;
protected:
	int nwrite,nsample;
	CRandom *randy;
	int nres;
	double *density,*ID;
	double denstot;
	string tmpfilename;
	FILE *tmpfile,*oscarfile,*input;
};


#endif
