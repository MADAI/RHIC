#ifndef __MUCALC_H__
#define __MUCALC_H__
#include "b3d.h"
using namespace std;

class CLocalInfo;
class CLocalSpeciesInfo;
typedef multimap<int,CLocalSpeciesInfo *> CLSIMap;

class CLocalSpeciesInfo{
public:
	string name;
	vector<double> T,ur,mu,Pr,E,M;
	vector<int> N;
	CLocalSpeciesInfo(string name_set);
	void Zero();
	void MuCalc();
	void MuCalc_PionsWithBose();
	void MuPrint();
	static CB3D *b3d;
};

class CLocalInfo{
public:
	static int NRBINS,NETEVENTS,IMUCALC;
	static double DELR2;
	static CB3D *b3d;
	static bool printing;
	int itau;
	CLSIMap lsimap;
	CLocalSpeciesInfo *pion;
	CLocalSpeciesInfo *kaon;
	CLocalSpeciesInfo *nucleon;
	CLocalSpeciesInfo *delta;
	CLocalSpeciesInfo *lambda;
	CLocalSpeciesInfo *sigma;
	CLocalSpeciesInfo *sigmastar;
	CLocalSpeciesInfo *xsi;
	CLocalSpeciesInfo *xsistar;
	CLocalSpeciesInfo *omega;
	double GetRegenFactor(CPart *part1,CPart *part2);
	CLocalInfo(int itau);
	void MuCalc();
	void Zero();
	void CheckMap();
};

#endif