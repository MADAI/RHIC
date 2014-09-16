#ifndef __INCLUDE_SPECIES_H__
#define __INCLUDE_SPECIES_H__

#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

// CLASS CSpecies *********************************************
class CSpecies{
public:
	double *m,*degen,**Q;
	int Nspecies,Ncharges,*NID,**ID;
	void Print();
protected:
	void InitStandardHadrons();
	//CSpecies();
};

class CSpecies_StandardHadrons_Equil : public CSpecies {
public:
	CSpecies_StandardHadrons_Equil();
};

class CSpecies_StandardHadrons : public CSpecies {
public:
	CSpecies_StandardHadrons();
};

class CSpecies_StandardHadrons5Q : public CSpecies {
public:
	CSpecies_StandardHadrons5Q();
};

class CSpecies_PionsOnly : public CSpecies {
public:
	CSpecies_PionsOnly();
};

class CSpecies_RelGas : public CSpecies {
public:
	CSpecies_RelGas();
};

class CSpecies_QGP : public CSpecies {
public:
	CSpecies_QGP();
};

#endif
