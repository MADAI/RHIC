#ifndef __REGEN_H__
#define __REGEN_H__
#include "b3d.h"
#include "resonances.h"
using namespace std;

typedef map<double,CResInfo *> CDensMap;
typedef pair<double, CResInfo*> CDensPair;

class CMuTInfo;

// Information for specific temperature

class CRegenerate{
public:
	double MC_Nbar,MC_Ntarget;
	CRandom *randy;
	CResInfoMap baryonmap;
	CRegenerate();
	CResList *reslist;
	int NT;
	double Tmin,delT,TAU0,DELTAU;
	vector< array <double,4> > BDens;  // vector cover different temperatures
	double sigmavmax;
	vector< array <CDensMap,4> > BDensMap;
	
	// These are workspace variables, and refer specific cells, but are not saved
	CPart *part1,*part2;
	array<CPart*,5> product;
	double dNMax[4][4];
	double dNMaxNet,sigmavrel;
	int NK;
	
	void CalculateBDens(int iT);
	void GetdNMax(CMuTInfo *muTinfo);
	void GetP1P2SigmaV(CMuTInfo *muTinfo);
	
	bool CheckForRegeneration(CB3DCell *cell,CMuTInfo *muTinfo);
	bool GetBBbarResInfo(CMuTInfo *muTinfo,CResInfo *&resinfo1,CResInfo *&resinfo2);
	void FillOutPartPRInfo(CPart *part,CB3DCell *cell,CMuTInfo *muTinfo);

 	static CB3D *b3d;
};

#endif