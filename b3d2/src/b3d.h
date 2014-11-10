#ifndef __B3D_H__
#define __B3D_H__
//#define __B3D_USE_HDF5__

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <complex>
#include <cstdio>
#include <sys/stat.h>
#include <ctime>
#include <vector>
#include <array>
#include <unordered_map>
#include "part.h"
#include "coralutils.h"
#include "resonances.h"
#include "inelastic.h"
#include "sampler.h"
#include "joshconvert.h"
#include "hydrotob3d.h"
#include "mutinfo.h"
#include "seinfo.h"
#include "regen.h"
#include "cell.h"

using namespace std;

class CAction; class CB3D; class CB3DCell; class Csampler; class CHYDROtoB3D; class CLocalInfo;
class CRegenerate; class CSEInfo;

typedef unordered_map<long int,CPart *> CPartMap;
typedef multimap<double,CAction *> CActionMap;
typedef pair<long int,CPart*> CPartPair;
typedef pair<double,CAction*> CActionPair;

class CB3D{
public:
	parameterMap parmap;
	CPartMap DeadPartMap;
	CPartMap PartMap;		//!< A C++ map for active CPart objects in the model.
	CPartMap FinalPartMap;	//!< A C++ map that stores information about particles that have left the collision mesh).
	//!A C++ map for CAction objects
	/*!
	This map is used to schedule and organize the various actions that the model must perform in time order. It contains all actions (as CAction objects) that have yet to occur, and the map's key is the boost-invariant time \f$\tau\f$ at which the action is scheduled to occur.
	*/
	CActionMap ActionMap;
	CActionMap DeadActionMap; // action objects not in line to be processed
	CResList *reslist;	//!< The CResList instance for the mucalodel (dynamically allocated).
	CInelasticList *inelasticlist;	//!< The CInelasicList instance for the model (dynamically allocated).
	
	int NXY;	//!< Determines size of mesh. The mesh size is \f$(2NXY,2NXY, 2NETA)\f$.
	int NETA;
	int NRINGSMAX;
	int NPRCELLSMAX;
	int COOPERFRYE_WMAX;
	bool COOPERFRYE_CREATENEGPARTS;
	bool USE_OLD_SAMPLER;
	int HYDRO_OCTANT_SYMMETRY;
	double XYMAX,ETAMAX,DXY,DETA;
	Csampler *sampler;
	CHYDROtoB3D *hydrotob3d;
	bool BJORKEN,COLLISIONS,INELASTIC,HYDRO_PURE_BJORKEN,DENSWRITE,BARYON_ANNIHILATION,MUTCALC,SECALC;
	double ANNIHILATION_SREDUCTION;  // reduces annihilation cross section based on amount of strangeness
	int DENSWRITE_NTAU;
	int NBOSE;
	double DENSWRITE_DELTAU,MUTCALC_DELTAU;
	vector<vector<vector<CB3DCell *> > > cell;
	vector<vector<CMuTInfo *> > muTinfo;
	vector<double> annihilation_array;
	CSEInfo *SEinfo;
	void ReadHydroInput();
	CPart *GetDeadPart();
	void GetDeadParts(CPart *&part1,CPart *&part2);
	void GetDeadParts(array<CPart*,5> &product);
	CAction *GetDeadAction();
	vector<CLocalInfo *> localinfo;
	
	//!Constructor.
	/*!
	The constructor, which initializes all other elements of the model run.
	
	\param[in] run_name_set This is the "run name" read in from the command line.
	*/
	CB3D(); // this is a constructor which does nothing but create an object
	CB3D(string run_name_set); // this gets all arrays ready
	void InitCascade();
	void InitAnalysis();
	~CB3D();	//!< Destructor.
	double tau,TAUCOLLMAX;
	int itau;
	//
	// READ IN FROM PARAMETER FILE
	int NACTIONSMAX;
	int NPARTSMAX,nbaryons,npartstot,nactionstot;
	int NSAMPLE,NSCATT_MAX;
	int DELNPARTSTOT,DELNACTIONSTOT;
	bool BINARY_RW;
	double SIGMAMAX,SIGMADEFAULT, SIGMAINELASTIC, Q0; // cross sections in sq. fm
	string input_dataroot;
	string output_dataroot;
	string run_name,qualifier;

	string oscarfilename;
	FILE *oscarfile;
	//
	void SetQualifier(string qualifier_set);
	void MovePartsToFinalMap();
	double WriteOSCAR(int ievent);  // returns dnch/deta
	void WriteDens();
	void WriteAnnihilationData();
	
	void FindAllCollisions();
	void FindAllCellExits();
	void PerformAllActions();
	void Reset();
	void KillAllActions();
	void KillAllParts();
	void SplitPart(CPart *part1,CPart *part2);

	void AddAction_Activate(CPart *part);
	void AddAction_Decay(CPart *part,double taudecay);
	void AddAction_Collision(CPart *part1,CPart *part2,double tau,double pibsquared);
	void AddAction_ResetCollisions(double taureset);
	//void AddAction_SwallowParticles(double tau_breakup);
	void AddAction_ExitCell(CPart *part);
	void AddAction_DensCalc(double tauwrite);
	void AddAction_MuTCalc(double taucalc);
	void AddAction_SECalc(double taucalc);

	void ListFutureCollisions();
	void PrintPartList();
	void PrintMuTInfo();
	void WriteWeights();
	void IncrementWeightArrays();

	bool FindCollision(CPart *part1,CPart *part2,double &taucoll);
	void Decay(CPart *mother,int &nbodies,array<CPart *,5> &daughter);
	double CalcSigma(CPart *part1,CPart *part2);

	CRandom *randy;

	void PrintActionMap(CActionMap *actionmap);

	double GetPiBsquared(CPart *part1,CPart *part2);
	int Collide(CPart *part1,CPart *part2,int &nproducts,array<CPart*,5> &product,double pibsquared); // will collide if sigma>scompare
	void Scatter(CPart *part1,CPart *part2,CPart *part3,CPart *part4);
	bool Merge(CPart *part1,CPart *part2,CPart *part3,CResInfo *resinfo);
	void InelasticScatter(CPart *part1, CPart *part2,CPart *part3,CPart *part4,CInelasticInfo inelinfo);
	double GetAnnihilationSigma(CPart *part1,CPart *part2,double &vrel);
	int Annihilate(CPart *part1,CPart *part2,int &nproducts,array<CPart *,5> &product);

	void CheckActions();
	bool ERROR_PRINT;

	long long int nscatter,n ,nmerge,nswallow,npass,nexit,nactivate,nannihilate,nregenerate,nactionkills;
	long long int nactions,ninelastic, ncollisions,ndecay,ncheck;

	void freegascalc_onespecies(double m,double t,double &p,double &e,double &dens,double &sigma2,double &dedt);
	
	bool BALANCE_DECAY,BALANCE_CALC;
	int ibalmax;
	void ReadBalanceParts();
	
	// These are used for Analysis
	CPart **partarray;
	bool CALCGARRAYS;
	bool STAR_ACCEPTANCE;
	CRegenerate *regen;
	int ReadOSCAR(int ievent);
	void ReadOSCARHeader();
	int DecayParts(int nparts);
	double CalcSpectra_PHENIX();
	void CalcSpectra_PHENIXppbar();
	double CalcSpectra_ALICE();
	void Calc3DSpectra();
    void CalcBalance();
	void CalcV2_STAR();
	void CalcV2_ALICE();
	void CalcHBT_STAR();
	void CalcHBT_ALICE();
	void CalcRealityDiff();
	void CalcGamma();
	double legendre(int ell,double x);
	void Consolidate(string run_name);
};
//!An action in the CB3D model.
/*!
\version 1.0
\author Scott Pratt
\date March 2011

This class handles any actions that the model takes during execution. Examples of "actions" that the model takes are a resonance decaying, a particle crossing a cell boundary, a collision, new particles being generated, etc. In this way, a complex system of interacting particles is reduced to a scheduled list of actions. Scheduling is handled using a C++ map container of CAction objects, keyed by the boost-invariant time tau (\f$\tau\f$) at which they are scheduled to occur. Note that this map is revised consistently, as future actions often change dramatically as a result of the current action.

Actions are allocated and are moved from the map of future actions (CB3D::ActionMap) to the list of completed actions (CB3D::DeadActionMap) once they have been performed.
*/
class CAction{
public:
	double tau;
	double pibsquared;
	int listid;
	double key;
	int type; // =0 for activation, 1 for decay, 2 for collision, ....  6 for ExitCell
	// These are the particles in the action
	CPartMap partmap;

	void Kill();
	void AddPart(CPart *partptr);
	void Print();

	void Perform();
	void PerformDensCalc();
	void PerformMuTCalc();
	void PerformSECalc();
	void PerformActivate();
	void PerformExitCell();
	void PerformDecay();
	void PerformCollide();
	void PerformCollide_BALANCE();
	void PerformResetCollisions();
	array <CResInfo *,5> daughterresinfo;
	//void PerformSwallowParticles();
	CAction();
	CAction(int keyset);
	~CAction();

	static CB3D *b3d;

	CActionMap::iterator GetPos(CActionMap *actionmap);
	void MoveToActionMap();
	void RemoveFromActionMap();
	void AddToMap(CActionMap *newmap);
	void AddToMap(CActionMap::iterator guess,CActionMap *newmap);
	void CheckPartList();
	CActionMap *currentmap;
	array<CPart *,5> product;
};



#endif
