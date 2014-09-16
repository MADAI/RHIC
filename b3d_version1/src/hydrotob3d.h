#ifndef __HYDROTOB3D_H__
#define __HYDROTOB3D_H__
//#define VIZWRITE

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <complex>
#include <cstdio>
#include <list>
#include <sys/stat.h>
#include "part.h"
#include "coralutils.h"
#include "H5Cpp.h"
#include "resonances.h"
#include "bjmaker.h"
#include "inelastic.h"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

using namespace std;

class CRing; class CResList; class CPRCell;

class CHYDROtoB3D{
public:
	//CHYDROtoB3D();
	// lambdafact is the ratio s.t., lambda_ij = pi_ij / lambdafact,
	// where pi_ij is the shear tensor in matter frame and lambda_ij tells how the momenta are scaled 
	// during collision time  p_i = ptilde_i + lambda_ij ptilde_j, where ptilde is generated isotropically
	double T,ETAMAX,DETA;
	double meanpt,meanu,etot,vtot;
	int normpt,NETA;
	bool initialization;
	// will generate particles with -eta_max < eta < eta_max
	CResList *reslist;
	double epsilon,P,lambdafact;
	int MakeEvent();
	int MakeEventPR();
	int MakeEvent3D();
	int GenerateParticles(int iring);
	int GenerateParticlesPR(int iprcell);
	int GenerateParticles3D(int ieta,int iring);
	void freegascalc_onespecies(double m,double t,double &p,double &e,double &dens,double &sigma2,double &dedt);
	void Init();
	void InitPR();
	void Init3D();
	void ReadInput();
	void ReadInputPR();
	void ReadInput3D();
	CB3D *b3d;
	void GetLambdaFact();
	void TestLambdaFact();
protected:
	int MC_NWrite;
	double MC_Ntarget,MC_Nbar,nsample,Ncheck;
	CRandom *randy;
	int nres,nrings,nprcells,*nrings3d;
	double *density,*ID;
	string tmpfilename;
	FILE *tmpfile,*input;
	CRing *ring;
	CRing **ring3d;
	CPRCell *prcell;
	void GetRingInfo();
	void GetRingInfo3D();
	void ReadHeader(FILE *);
	void ReadHeaderPR(FILE *);
	void ReadHeader3D(FILE *);
};

// This fills out Pi[4][4] given Pi_xx, Pi_yy, Pi_xy, u_x and u_y (all specified in lab frame)
void FillOutPi(double **Pi,double pixx,double piyy,double pixy,double ux,double uy);

// Stores info about a ring that breaks up at specific tau
class CRing{
public:
	CRing();
	void Read(FILE *fptr);
	void Read3D(FILE *fptr);
	void FillOutArrays(double *x,double *y,double *uxi,double *uyi,double *dToverH_xx,double *dToverH_xy,double *dToverH_yy);
	int nphi;
	int nread;
	double tau,etamin,etamax;
	double *r,*ux,*uy;
	double ***dToverH;
};

class COSUHydrotoB3D;
class CFreezeoutCell{
public:
	double P,epsilon,T,u[4],r[4],Omega[4],**dToverH; //pitilde[alpha][beta]/(epsilon+H);
	double lambda[4][4],wmax;
	bool Read(FILE *file);
	bool Read3D(FILE *file);
	int GenerateParticles();
	int GenerateParticles3D();
	CFreezeoutCell();
	static CB3D *b3d;
	static COSUHydrotoB3D *osuhydrotob3d;
};

class COSUHydrotoB3D{
public:
	bool initialization;
	int HYDRO_OCTANT_SYMMETRY;
	vector<CFreezeoutCell> cell;
	void Init();
	void ReadInput();
	int MakeEvent();
	double MC_Ntarget,MC_Nbar,MC_NWrite;
	int NRES,NCELLS;
	double *density,*ID,netdens;
	CRandom *randy;
	void freegascalc_onespecies(double m,double t,double &p,double &e,double &dens,double &sigma2,double &dedt);
	CResList *reslist;
	double T,epsilon,P,lambdafact;
	int NSAMPLE;
	double ETAMAX;
	CB3D *b3d;
	double meanpt,meanu,etot,vtot;
	long long int normpt;
	void GetLambdaFact();
	void TestLambdaFact();
};

class CPRCell{
public:
	CPRCell();
	static double T,epsilon,P,omega_xyspacing,omega_tauspacing,etamax;
	void Read(FILE *fptr);
	double x,y,ux,uy,eta,tau;
	double Omega_x,Omega_y,Omega_0;
	double dToverH[4][4];
	void GetPiTildeOverH(double pixx,double pixy,double piyy);
};

#endif
