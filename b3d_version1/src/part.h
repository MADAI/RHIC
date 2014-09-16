#ifndef __PART_H__
#define __PART_H__

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <list>
#include <map>
#include <sys/stat.h>
#include "H5Cpp.h"
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

using namespace std;

class CAction; class CB3D; class CB3DCell; class CResInfo; class CPart; class CPartH5;

typedef multimap<int,CPart *> CPartMap;
typedef multimap<double,CAction *> CActionMap;
typedef pair<int,CPart*> CPartPair;
typedef pair<double,CAction*> CActionPair;

//!A particle in the CB3D model.
/*!
 \version 1.0
 \author Scott Pratt
 \date March 2011
 
 This class generates the particles used in the CB3D model. Note that resonance information (stored in CResInfo objects) is different than the actual particles used; instead, a CPart object contains a pointer to a CResInfo object corresponding to the resonance it represents. In addition, the CPart object contains information about is coordinate and momentum 4 vectors, as well as its rapidity. Finally, it contains methods to add and remove actions (CAction objects) for the particle.
 
 In the CB3D model, the total number of particles is fixed (it is set by the CB3D::NPARTSMAX) parameter. As an attempt to improve performance, all of the particles are allocated in one memory block in the CB3D constructor, and stored in the "dead" particle map (CB3D::DeadPartMap). During model function, particles are intialized by setting their relevant parameters and moving them to the "live" particle map (CB3D::PartMap). Once the particle moves outside the outer boundary of the model space, it is transferred to the output particle map (CB3D::FinalPartMap).
 */

class CPart{
public:
	CPart();
	~CPart();
	CB3DCell *cell,*nextcell;
	double tau0,tau_lastint,tauexit,taudecay;
	double y,eta;
	double *p,*r,mass;
	int listid;
	int actionmother; //refers to action from which particle was created
	CResInfo *resinfo;

	void Propagate(double tau);
	void FindDecay();
	void FindCellExit();
	void Init(int ID,double x,double y,double tau,double eta,double px,double py,double mass,double rapidity);
	void Init_NoColls(int ID,double x,double y,double tau,double eta,double px,double py,double mass,double rapidity);
	void InitH5(CPartH5*);
	void Setp0();
	void Copy(CPart *part);
	void Copy(CPartH5 *parth5);
	double GetMass();
	void CyclicReset();
	int key;
	void SetInitialKey();

	void Kill();
	void FindActions();
	void SubtractAction(CAction *actionptr);
	void AddAction(CAction *actionptr);
	void KillActions();
	void Print();
	void CheckMap(CPartMap *expectedpartmap);
	void ChangeMap(CPartMap *newmap);
	void BjorkenTranslate();
	void BjorkenUnTranslate();
	void Boost(double *u);
	void BoostP(double *u);
	void BoostR(double *u);
	//~CPart();

	// These are the actions involving these particles
	CActionMap actionmap;

	CPartMap *currentmap; // PartList for a Cell, or b3d->DeadPartList
	CB3DCell *FindCell();

	static CB3D *b3d;
	double GetEta(double tau);
	double GetPseudoRapidity();
	double GetMT();
	void SetY();
	void SetEta(double neweta);
	void FindCollisions();

	CPartMap::iterator GetPos(CPartMap *pmap);
	CPartMap::iterator DeleteFromMap(CPartMap *partmap);
	void ChangeCell(CB3DCell *newcell);
	void RemoveFromCell();

	void ChangePartMap(CPartMap *newmap);
	CPartMap::iterator DeleteFromCurrentMap();
	void AddToMap(CPartMap *newmap);
	void AddToMap(CPartMap::iterator guess,CPartMap *newmap);
	bool active;
};


class CPartH5{
public:
	int listid,ID;
	double px,py,rapidity,mass,x,y,eta,tau;
};

#endif
