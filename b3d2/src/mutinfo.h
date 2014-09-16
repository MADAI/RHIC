#ifndef __MUCALC_H__
#define __MUCALC_H__
#include "b3d.h"
using namespace std;

class CMuTInfo{
public:
	CMuTInfo();
	// vectors hold information for different times
	vector<double> Pxpi,Pypi,Epi,Mpi,PxK,PyK,EK,MK;
	vector<int> Npi,NK;
	// this is information for latest time
	double Tpi,mupi,uxpi,uypi,TK,muK,uxK,uyK;
	
	void Zero(); // sets arrays to zero
	void UpdateNMPE(CB3DCell *cell);
	void MuTCalc();
	void Print();
	void FindMuTUxUy(double tau,int N,double E,double M,double Px,double Py,double degen,double &T,double &mu,double &ux,double &uy);
	void FindMuTInfo_pi(int itau);
	void FindMuTInfo_K(int itau);
	//void MuTCalc_PionsWithBose();
	static CB3D *b3d;
	static double DELTAU;
	static int NTAU;
	static int NETEVENTS;
};

#endif