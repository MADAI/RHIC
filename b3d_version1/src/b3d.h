#ifndef __B3D_H__
#define __B3D_H__
#define __B3D_USE_HDF5__

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <complex>
#include <cstdio>
#include <list>
#include <sys/stat.h>
#include <ctime>
#include "part.h"
#include "coralutils.h"
#include "resonances.h"
#include "bjmaker.h"
#include "inelastic.h"
#include "hydrotob3d.h"

#ifdef __B3D_USE_HDF5__
#include "H5Cpp.h"
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif
#endif

using namespace std;

class CAction; class CB3D; class CB3DCell; class CHYDROtoB3D; class COSUHydrotoB3D;

typedef multimap<int,CPart *> CPartMap;
typedef multimap<double,CAction *> CActionMap;
typedef pair<int,CPart*> CPartPair;
typedef pair<double,CAction*> CActionPair;

//!The main model routine.
/*!
\version 1.0
\author Kevin Novak, Scott Pratt
\date March 2011

This is the main routine for the running of the model.

\todo More detail in the general description
\todo FinalPartMap
*/

class CB3D{
public:
	//!A C++ map containing all information about the parameter set.
	/*!
	The parameterMap type is a custom version of the generalized C++ map container. It contains various methods for storing and returning almost all data types. The parameter map is designed to be implemented using the fixed.param and stats.param convention discussed in the User's Manual.
	*/
	parameterMap parmap;
	//!A C++ map for dead CPart objects.
	/*!
	This map is populated by "dead" particles; that is, particles which have been removed from the functional particle map.
	\sa PartMap
	*/
	CPartMap DeadPartMap;
	CPartMap PartMap;		//!< A C++ map for active CPart objects in the model.
	CPartMap FinalPartMap;	//!< A C++ map that stores information about particles that have left the model (hit the outer edge).
	//!A C++ map for CAction objects
	/*!
	This map is used to schedule and organize the various actions that the model must perform in time order. It contains all actions (as CAction objects) that have yet to occur, and the map's key is the boost-invariant time \f$\tau\f$ at which the action is scheduled to occur.
	*/
	CActionMap ActionMap;
	//!A C++ map for CAction objects that have already occured.
	/*!
	\sa ActionMap
	*/
	CActionMap DeadActionMap;
	CResList *reslist;	//!< The CResList instance for the model (dynamically allocated).
	CInelasticList *inelasticlist;	//!< The CInelasicList instance for the model (dynamically allocated).
	//!A dynamically allocated array of CPart objects.
	/*!
	This array is allocated to be the size of the parameter NPARTSMAX (the highest number of particles in the model), and is used only to allocate the memory that the CPart objects use. It was thought that if the total number of particles (dead, alive, and final) remains constant and allocated at the initialization of the model would be a more efficient method.
	*/
	CPart **partarray;
	//!A dynamically allocated array of CAction objects.
	/*!
	Much like the part array, this array is dynamically allocated (using NACTIONSMAX) and used to allocate the memory for all CAction objects during the model's running.
	*/
	CAction **actionarray;
	
	int NXY;	//!< Determines size of mesh. The mesh size is \f$(2NXY,2NXY, 2NETA)\f$.
	int NETA;
	int NRINGSMAX;
	int NPRCELLSMAX;
	bool OSUHYDRO,PRHYDRO;
	int HYDRO_OCTANT_SYMMETRY;
	double XYMAX,ETAMAX,DXY,DETA;
	CHYDROtoB3D *hydrotob3d;
	COSUHydrotoB3D *osuhydrotob3d;
	CBjMaker bjmaker;
	bool BJORKEN,COLLISIONS,INELASTIC,HYDRO_PURE_BJORKEN,DENSWRITE,ANNIHILATION_CHECK;
	double ANNIHILATION_SREDUCTION;  // reduces annihilation cross section based on amount of strangeness
	int DENSWRITE_NTAU;
	double DENSWRITE_DELTAU;
	CB3DCell ****cell;
	double *annihilation_array;
	void ReadHydroInput();
	int HydrotoB3D();
	
	//!Constructor.
	/*!
	The constructor, which initializes all other elements of the model run.
	
	\param[in] run_name_set This is the "run name" read in from the command line.
	*/
	CB3D(); // this is a constructor which does nothing but create an object
	CB3D(string run_name_set); // this gets all arrays ready
	void InitArrays();
	~CB3D();	//!< Destructor.
	double tau,TAUCOLLMAX;
	int itau;
	int ievent_write,ievent_read;
	//
	// READ IN FROM PARAMETER FILE
	int NACTIONSMAX;
	int NPARTSMAX,nbaryons;
	double SIGMAMAX,SIGMADEFAULT, SIGMAINELASTIC, Q0; // cross sections in sq. fm
	string input_dataroot;
	string output_dataroot;
	string run_name,qualifier;
	CompType *ptype;

	#ifdef __B3D_USE_HDF5__
	string outfilename;
	H5File *h5outfile, *h5infile;
	int ReadDataH5(int ievent);
	double WriteDataH5(); // returns dnch/deta
	#endif
	string oscarfilename;
	FILE *oscarfile;
	int NACTIONS;
	int NSAMPLE;
	//
	void SetQualifier(string qualifier_set);
	void MovePartsToFinalMap();
	double WriteOSCAR();  // returns dnch/deta
	void WriteDens();
	void WriteAnnihilationData();
	
	void FindAllCollisions();
	void FindAllCellExits();
	void PerformAllActions();
	void Reset();
	void KillAllActions();
	void KillAllParts();

	void AddAction_Activate(CPart *part);
	void AddAction_Decay(CPart *part,double taudecay);
	void AddAction_Collision(CPart *part1,CPart *part2,double tau);
	void AddAction_ResetCollisions(double taureset);
	//void AddAction_SwallowParticles(double tau_breakup);
	void AddAction_ExitCell(CPart *part);
	void AddAction_DensCalc(double tauwrite);

	void ListFutureCollisions();
	void PrintPartList();

	bool FindCollision(CPart *part1,CPart *part2,double &taucoll);
	void Decay(CPart *&mother,int &nbodies, CPart **&daughter);
	double CalcSigma(CPart *part1,CPart *part2);

	CRandom *randy;

	void PrintActionMap(CActionMap *actionmap);

	double GetPiBsquared(CPart *part1,CPart *part2);
	int Collide(CPart *part1,CPart *part2); // will collide if sigma>scompare
	void Scatter(CPart *part1,CPart *part2);
	bool Merge(CPart *part1,CPart *part2,CResInfo *resinfo);
	void InelasticScatter(CPart *part1, CPart *part2, CInelasticInfo inelinfo);
	int Annihilate(CPart *part1,CPart *part2);

	void CheckActions();
	bool ERROR_PRINT;

	long long int nscatter,ndecay,nmerge,nswallow,npass,nexit,nactivate,ncheck,nannihilate,nactionkills;
	long long int nactions,ninelastic, ncollisions;

	void freegascalc_onespecies(double m,double t,double &p,double &e,double &dens,double &sigma2,double &dedt);
};
//!An action in the CB3D model.
/*!
\version 1.0
\author Scott Pratt
\date March 2011

This class handles any actions that the model takes during execution. Examples of "actions" that the model takes are a resonance decaying, a particle crossing a cell boundary, a collision, new particles being generated, etc. In this way, a complex system of interacting particles is reduced to a scheduled list of actions. Scheduling is handled using a C++ map container of CAction objects, keyed by the boost-invariant time tau (\f$\tau\f$) at which they are scheduled to occur. Note that this map is revised consistently, as future actions often change dramatically as a result of the current action.

Much like particles and CPart objects, the total number of actions is also a constant (set by CB3D::NACTIONSMAX). Actions are allocated in one memory block in the CB3D constructor, and are moved from the map of future actions (CB3D::ActionMap) to the list of completed actions (CB3D::DeadActionMap) once they have been performed.
*/
class CAction{
public:
	double tau;
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
	void PerformActivate();
	void PerformExitCell();
	void PerformDecay();
	void PerformCollide();
	void PerformResetCollisions();
	//void PerformSwallowParticles();
	CAction();
	~CAction();

	static CB3D *b3d;

	CActionMap::iterator GetPos(CActionMap *actionmap);
	void MoveToActionMap();
	void RemoveFromActionMap();
	void AddToMap(CActionMap *newmap);
	void AddToMap(CActionMap::iterator guess,CActionMap *newmap);
	void CheckPartList();
	CActionMap *currentmap;
};
//!A cell in the expanding cell mesh
/*!
\version 1.0
\author Scott Pratt
\date March 2011

In the CB3D model, the model space is expressed as a mesh grid of cells that expand as time propogates. The mesh is populated by cells, which are CB3DCell objects. This class keeps track of the particles populating it (in a particle map), as well as its spatial dimensions and neighbors. The neighbors are especially relevant, as actions such as collisions are scheduled by checking against particles inside its the current cell, as well as all neighboring cells.
*/

class CB3DCell{
public:
	class CB3DCell *neighbor[3][3][3];
	int ix,iy,ieta;
	class CB3DCell *creflection;
	int ireflection;
	double xmin,xmax,ymin,ymax,etamin,etamax;
	CPartMap partmap;
	void PrintPartMap(CPartMap *partmap);
	void KillAllParts();
	void ReKeyAllParts();
	void Print();
	double *dens;
	
	CB3DCell(double xmin,double xmax,double ymin,double ymax,double etamin,double etamax);
	static CB3D *b3d;
};


#endif
