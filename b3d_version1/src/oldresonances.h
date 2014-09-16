#ifndef __RESONANCES_H__
#define __RESONANCES_H__

#include "b3d.h"

using namespace std;

class CResInfo;

class CB3D;
//!Information about the possible decay channels of a resonance		
/*!
\author Scott Pratt
\version 1.0
\date March 2011

This class is used to store information about one of the various decay branches of a resonance when read in from the file. This is mostly a intermediate class, in between the reading in of decays and the creation of CMerge objects. As of version 1, it is restricted to (up to) five particle decays.

\sa CMerge
*/
class CBranchInfo{
public:
	CBranchInfo *nextbptr;	//!< Pointer to next CBranchInfo object (the next decay channel.)
	int nbodies;	//!< Number of bodies in the decay. (2-5)
	CResInfo *resinfoptr[5];	//!< Array of pointers to the CResInfo objects in the decay.
	double branching;	//!< Branching fraction of decay channel.
	CBranchInfo();	//!< Constructor
};
//!Information about an excitation of a resonance (two particles merging into one).
/*!
\author Scott Pratt, Kevin Novak
\version 1.0
\date March 2011

This class stores the information about the excitation of a resonance, one of the possible actions that occur during the collision routine. The are developed by considering the decay channels of the various resonances (as CBranchInfo objects), and exploiting the time reversibility of decays (i.e. if a resonance can decay into two particles, the two particles can excite the same resonance). These classes are used to during the collision and decay routines.
*/
class CMerge{
public:
	CMerge(CResInfo *resinfo,double branching, int L);	//!< Constructor
	CResInfo *resinfo;	//!< Pointer to the CResInfo object corresponding to the excited resonance.
	int L;	//!< Angular momentum value of the exit channel. (Used in collide routine.)
	double branching;	//!< Branching fraction of excitation channel.
	CMerge *next;	//!< pointer to next CMerge object in the MergeArray routine.
};
//!Information about one of the resonances in CB3D
/*!
\author Scott Pratt, Kevin Novak
\version 1.0
\date March 2011

This class is used to provide a more detailed structure for the resonances in CB3D. Note that there are 319 resonances included in the model; the resonances are read in and stored in a singly linked list, and then referenced by pointer throughout the rest of the model function. Note that both stable and decaying particles are treated identically.

\todo DecayGetResInfoptr, DecayGetResInfoptr_minmass
*/
class CResInfo{
public:
	//!Resonance number
	/*!
	The resonances are, by convention, numbered from 0 to 318 in the order in which I read in. These serve to allow the merging exit channel array (MergeArray) and inelastic exit channel array (if used) (InelasticArray) to be referenced using the ires values on the incoming particles.
	 
	\sa CResList::MergeArray
	\sa CInelasticList::InelasticArray
	*/
	int ires;
	double mass;	//!< Pole mass of resonance.
	double spin;	//!< Spin of resonance.
	double width;	//!< Pole width of resonance.
	double minmass;	//!< Minimum invariant mass required to excite resonance.
	string name;	//!< Name of resonance. (as read in from file)
	int code;		//!< Particle code of resonance. (from Particle Data Group datafile)
	int charge;		//!< Charge of particle resonance.
	int strange;	//!< Strangeness of resonance.
	int baryon;		//!< Bayon number of resonance.
	int count;		//!< Depreciated.
	int G_Parity;	//!< G Parity eigenvalue, or NULL otherwise.
	bool decay;		//!< Whether the resonance decays or not.
	CResInfo *nextResInfoptr;	//!< Pointer to next CResInfo object in linked list.
	CBranchInfo *firstbptr;	//!< Pointer to first decay channel of particle.
	CBranchInfo	*bptr_minmass;	//!< Pointer to decay channel with minimum mass.
	void Print();	//!< Print out information.
	void DecayGetResInfoptr(int &nbodies,CResInfo **&daughterresinfoptr);
	void DecayGetResInfoptr_minmass(int &nbodies,CResInfo **&daughterresinfoptr);
	bool CheckForDaughters(int code);
	CResInfo();	//!< Constructor.
	static CRandom *ranptr;	//!< A dynamically allocated random number generator.
	
};
//!This class is used for storing and organizing CResInfo objects.
/*!
\version 1.0
\author Kevin Novak
\date March 2011

Much like the CInelasticList class is used to store and organize CInelasticInfo objects for inelastic exit channels, the CResList class is used to  organize the CResInfo linked list, and store and organize the merge exit channels (CMerge objects).

\todo freegascalc_onespecies, CalcEoS, PrintYields
*/
class CResList{
public:
	int NResonances;	//!< The number of resonances read in.
	CResList();	//!< Constructor.
	CResList(parameterMap* parmap_in);
	//!Pointer to the first CResInfo pointer (the head of the linked list).
	CResInfo *GfirstResInfoptr;
	//!Given a particle code, gets a pointer to the appropriate CResInfo object.
	/*!
	As the CResInfo objects are stored in a singly-linked list, this method is mostly used to reference specific resonances in the list, especially when the exact order of the list is not known. Given a particle code and a pointer, the list is iterated over until the pointer is now pointed at the correct resonance.
	\sa CResInfo::code
	\param[in] code Particle code to be found.
	\param[in,out] resinfoptr CResInfo pointer to be pointed.
	*/
	void GetResInfoptr(int code,CResInfo *&resinfoptr);
	//!Read in the resonance and decay file information
	/*!
	This method serves to read in the resonance and decay information from the external file. Additionally, it creates the linked list of CResInfo objects and creates and organizes the array of CMerge objects (MergeArray).
	
	\sa MergeArray
	*/
	void ReadResInfo();
	void PrintYields();
	void CalcEoS(double T_0,double T_F,double delT);
	void CalcEoSandChi(double T);
	void freegascalc_onespecies(double m,double t,double &p,double &e,double &dens,double &sigma2,double &dedt);
	parameterMap *parmap;	//!< Pointer to the CB3D parameter map.
	//!An array of CMerge objects, storing merge exit channels.
	/*!
	This is a dynamically allocated 2D array of linked lists of CMerge objects, corresponding to the possible merging exit channels. The array is referenced by the \a ires values of the two incoming particles.
	
	\sa CMerge::ires
	*/
	CMerge ***MergeArray;
	static CB3D *b3d;	//!< Pointer to the CB3D object the CResList object belongs to.
	
};


#endif
