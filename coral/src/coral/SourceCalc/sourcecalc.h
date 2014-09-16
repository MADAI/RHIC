#ifndef __INCLUDE_SOURCECALC_H__
#define __INCLUDE_SOURCECALC_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <fstream>
#include "coral.h"
//#include "arrays.h"
//#include "parametermap.h"
//#include "misc.h"
//#include "random.h"
//#include "constants.h"
//#ifdef __B3D_USE_HDF5__
//#include "H5Cpp.h"
#include "part.h"
//#ifndef H5_NO_NAMESPACE
//using namespace H5;
//#endif
//#endif

using namespace std;

class CSourceCalc{
public:
	parameterMap spars;
	virtual void CalcS(CCHArray *A);
	virtual void CalcS(int lx,int ly,int lz,CCHArray *A);
	virtual void CalcS(CMCList *&lista,CMCList *&listb);
	virtual void CalcS(C3DArray *threed);
	virtual void GaussCFCalc(C3DArray *cf3d);

	void CombineMCLists(CMCList *lista,CMCList *listb,CCHArray *A);
	void CombineMCLists(CMCList *lista,CMCList *listb,CCHArray *A,int NMC);
	void CombineMCLists(CMCList *lista,CMCList *listb,C3DArray *threed);
	
	void CombineMCPRLists(CMCPRList *lista,CMCPRList *listb,CCHArray *A);
	void CombineMCPRLists(CMCPRList *lista,CMCPRList *listb,CCHArray *A,int NMC);
	void CombineMCPRLists(CMCPRList *lista,CMCPRList *listb,C3DArray *threed);
	
	void ReadSPars(char *sparsfilename);
	void NormCheck(CCHArray *A);
	double GetNorm(CCHArray *A);
	void NormCheck(C3DArray *threed);
	double GetNorm(C3DArray *threed);
	void CalcEffGaussPars(CCHArray *A);
	void CalcEffGaussPars(CCHArray *A,double &Rx,double &Ry,double &Rz,
		double &Xoff,double &Yoff,double &Zoff);
	void CalcEffGaussPars(CMCPRList *list,double &Rx,double &Ry,double &Rz,
												double &Xoff,double &Yoff,double &Zoff);
	void CalcEffGaussParsQ2(CMCPRList *list,double &Rx,double &Ry,double &Rz);
	void CalcEffGaussParsPureBose(CMCPRList *list,double &lambda,double &Rx,double &Ry,double &Rz);
	void CalcEffGaussParsPureBose_GetM(double ***C,double *x,int nxyz,double **M);
	void CalcEffGaussParsPureBose_GetYChi2(double ***C,double *x,int nxyz,double *y,double &chi2);
	CSourceCalc();
	CSourceCalc(string sparsfilename);
	CRandom *randy;
	virtual ~CSourceCalc(){};
private:

};

class CSourceCalc_GX1D : public CSourceCalc{
/*
S = { N_G * lambda_G * exp(-r^2/4R^2)
+ N_X * lambda_X * exp(-[r^2/X^2+4R^4/X^4]) } * [r^2/(r^2+a^2)]^(L/2)
where N_x and N_g are constants that make individual contributions
integrate to unity for L=0
*/
public:
	CSourceCalc_GX1D();
	void InitSPars();
	void SetSPars(double lambda,double Xfrac,double R,double X,double a);
	void CalcS(int lx,int ly,int lz,CCHArray *A);
	void CalcS(CCHArray *A);
	void CalcS(CMCList *&lista,CMCList *&listb){
		printf("CSourceCalc_GX1D::CalcS(CMCList *&lista,CMCList *&listb) undefined\n");
		exit(1);
	};
	void CalcS(C3DArray *threed){
		printf("CSourceCalc_GX1D::CalcS(C3DArray *threed) undefined\n");
		exit(1);
	};
	void GaussCFCalc(C3DArray *cf3d){
		printf("CSourceCalc_GX1D::CalcS(C3DArray *cf3d) undefined\n");
		exit(1);
	};
};

class CSourceCalc_Gaussian : public CSourceCalc{
public:
	CSourceCalc_Gaussian();
	void SetSPars(double lambdaset,double Rxset,double Ryset,double Rzset);
	void SetSPars(double lambdaset,double Rxset,double Ryset,double Rzset,
		double Xoffset,double Yoffset,double Zoffset);
	void SetSPars(double lambdaset,double Rxset,double Ryset,double Rzset,
		double Xoffset,double Yoffset,double Zoffset,
		double Euler_phiset,double Euler_thetaset,
		double Euler_psiset);
	void CalcS(CCHArray *A);
	void CalcS(C3DArray *threed);
	void CalcS(int lx,int ly,int lz,CCHArray *A){
		printf("CSourceCalc_Gaussian:CalS(int lx,int ly,int lz,CCHArray *A) undefined\n");
		exit(1);	
	};
	void CalcS(CMCList *&lista,CMCList *&listb){
		printf("CSourceCalc_Gaussian:CalS(CMCList *&lista,CMCList *&listb) undefined\n");
		exit(1);	
	};
	void GaussCFCalc(C3DArray *cf3d);

private:
	void CalcAlpha(double **alpha,CCHArray *A);
	void InitSPars();
	// S ~ exp{-alpha_ij (x_i-off_i)(x_j-off_j)} 
	// Calcs ignore zeroth components
};

class CSourceCalc_EllipticBlast : public CSourceCalc{
public:
	CSourceCalc_EllipticBlast();
	void CalcS(CCHArray *A);
	void CalcS(int lx,int ly,int lz,CCHArray *A){
		printf("CalcS(int lx,int ly,int lz,CCHArray *A) undefined\n");
		exit(1);
	};
	void CalcS(CMCList *&lista,CMCList *&listb){
		printf("CalcS(CMCList *&lista,CMCList *&listb) undefined\n");
		exit(1);
	};
	void CalcS(C3DArray *threed){
		printf("CalcS(C3DArray *threed) undefined\n");
		exit(1);
	};
	void SetSPars(double Rxset,double Ryset,
		double Tauset,
		double BetaXset,double BetaYset,
		double Tset,double Ptset,
		double Phiset,double EtaGset,
		double Maset,double Mbset);
	void SetSPars(double Rset,double Tauset,
		double Betaset,double Tset,double Ptset);
private:
	void Get_r(double *p,int nsample,double **r);
	void InitSPars();
};

class CSourceCalc_OSCAR : public CSourceCalc{
public:
	CSourceCalc_OSCAR();
	CSourceCalc_OSCAR(string sparsfilename);
	//void CalcS(CCHArray *A);
	void CalcS(CMCList *&lista,CMCList *&listb);
	void CalcS(CCHArray *A){
		printf("CSourceCalc_OSCAR:CalcS(CCHArray *A) undefined\n");
		exit(1);
	};
	void CalcS(int lx,int ly,int lz,CCHArray *A){
		printf("CSourceCalc_OSCAR:CalcS(int lx,int ly,int lz,CCHArray *A) undefined\n");
		exit(1);
	};
	void CalcS(C3DArray *threed){
		printf("CSourceCalc_OSCAR:CalcS(C3DArray *threed) undefined\n");
		exit(1);
	};
	void SetSPars(double Pt_set,double delpt_set,double phimin_set,double phimax_set,double ymin_set,double ymax_set);
	void SetIDs(int *idlista,int nida,int *idlistb,int nidb);
private:
	int *idlist_a,*idlist_b,nid_a,nid_b;
	void ReadR(double **ra,int &na,double **rb,int &nb);
	void InitSPars();
	bool IDMatch(int ident,int *idlist,int nid);
	bool Check(double *p,double *r,double m,double **ra,int &n);
};

class CSourceCalc_HDF5 : public CSourceCalc{
public:
	CSourceCalc_HDF5();
	CSourceCalc_HDF5(string sparsfilename);
	//void CalcS(CCHArray *A);
	void CalcS(CMCPRList *&lista,CMCPRList *&listb);
	void CalcS(CCHArray *A){
		printf("CSourceCalc_HDF5:CalcS(CCHArray *A) undefined\n");
		exit(1);
	};
	void CalcS(int lx,int ly,int lz,CCHArray *A){
		printf("CSourceCalc_HDF5:CalcS(int lx,int ly,int lz,CCHArray *A) undefined\n");
		exit(1);
	};
	void CalcS(C3DArray *threed){
		printf("CSourceCalc_HDF5:CalcS(C3DArray *threed) undefined\n");
		exit(1);
	};
	void SetSPars(double Pt_set,double delpt_set,double phimin_set,double phimax_set,double ymin_set,double ymax_set);
	void SetIDs(int *idlista,int nida,int *idlistb,int nidb);
private:
	int *idlist_a,*idlist_b,nid_a,nid_b;
	long long int ReadPR(double **pa,double **ra,int &na,double **pb,double **rb,int &nb);
	void InitSPars();
	bool IDMatch(int ident,int *idlist,int nid);
	bool Check(double *p,double *r,double m,double **pa,double **ra,int &n);
};

class CSourceCalc_HDF5_MultiBin : public CSourceCalc{
public:
	CSourceCalc_HDF5_MultiBin();
	CSourceCalc_HDF5_MultiBin(string sparsfilename);
	//void CalcS(CCHArray *A);
	void CalcS(CMCPRList ***&lista,CMCPRList ***&listb);
	int NPTBINS,NPHIBINS;
	double DELPT,PTMIN,PTMAX,DELPHI;
	void CalcS(CCHArray *A){
		printf("CSourceCalc_HDF5:CalcS(CCHArray *A) undefined\n");
		exit(1);
	};
	void CalcS(int lx,int ly,int lz,CCHArray *A){
		printf("CSourceCalc_HDF5:CalcS(int lx,int ly,int lz,CCHArray *A) undefined\n");
		exit(1);
	};
	void CalcS(C3DArray *threed){
		printf("CSourceCalc_HDF5:CalcS(C3DArray *threed) undefined\n");
		exit(1);
	};
	void SetSPars(double Pt_set,double delpt_set,double phimin_set,double phimax_set,double ymin_set,double ymax_set);
	void SetIDs(int *idlista,int nida,int *idlistb,int nidb);
	//CompType *ptype;
private:
	int *idlist_a,*idlist_b,nid_a,nid_b;
	void ReadPR(double ****pa,double ****ra,int **&na,double ****pb,double ****rb,int **&nb);
	void InitSPars();
	bool IDMatch(int ident,int *idlist,int nid);
	bool Check(double *p,double *r,double m,double **pa,double **ra,int &n);
};

class CSourceCalc_OSCAR_MultiBin : public CSourceCalc{
public:
	CSourceCalc_OSCAR_MultiBin();
	CSourceCalc_OSCAR_MultiBin(string sparsfilename);
	string OSCARfilename;
	//void CalcS(CCHArray *A);
	void CalcS(CMCPRList ***&lista,CMCPRList ***&listb);
	int NPTBINS,NPHIBINS;
	double DELPT,PTMIN,PTMAX,DELPHI;
	void CalcS(CCHArray *A){
		printf("CSourceCalc_OSCAR_MultiBin:CalcS(CCHArray *A) undefined\n");
		exit(1);
	};
	void CalcS(int lx,int ly,int lz,CCHArray *A){
		printf("CSourceCalc_OSCAR_MultiBin:CalcS(int lx,int ly,int lz,CCHArray *A) undefined\n");
		exit(1);
	};
	void CalcS(C3DArray *threed){
		printf("CSourceCalc_OSCAR_MultiBin:CalcS(C3DArray *threed) undefined\n");
		exit(1);
	};
	void SetSPars(double Pt_set,double delpt_set,double phimin_set,double phimax_set,double ymin_set,double ymax_set);
	void SetIDs(int *idlista,int nida,int *idlistb,int nidb);
	bool B3D_BINARY_FORMAT; // use binary format from B3D instead of oscar
	//CompType *ptype;
private:
	int *idlist_a,*idlist_b,nid_a,nid_b;
	long long int ReadPR(double ****pa,double ****ra,int **&na,double ****pb,double ****rb,int **&nb);
	void InitSPars();
	bool IDMatch(int ident,int *idlist,int nid);
	bool Check(double *p,double *r,double m,double **pa,double **ra,int &n);
};


class CSourceCalc_Blast : public CSourceCalc{
public:
	CSourceCalc_Blast();
	void CalcS(CCHArray *A);
	void CalcS(CMCList *lista,CMCList *listb);
	void CalcS(int lx,int ly,int lz,CCHArray *A){
		printf("CSourceCalc_Blast::CalcS(int lx,int ly,int lz,CCHArray *A) undefined\n");
		exit(1);
	};
	void CalcS(C3DArray *threed){
		printf("CSourceCalc_Blast::CalcS(C3DArray *threed) undefined\n");
		exit(1);
	};
	void SetSPars(double lambdaset,double Rset,double Tauset,double DelTauset,double Betaset,double Tset,double Ptset,double EtaGset,double Maset,double Mbset);
	void SetSPars(double lambdaset,double Rset,double Tauset,double DelTauset);
private:
	void GetMCList(double *p,CMCList *mclist);
	void InitSPars();
	double GetTau(double tau0,double deltau);
};


#endif
