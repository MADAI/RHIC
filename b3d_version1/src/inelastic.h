#ifndef __INELASTIC_H__
#define __INELASTIC_H__

#include "b3d.h"

using namespace std;
//!Stores information about inelastic scattering exit channels
/*!
\author Kevin Novak
\version 1.0
\date March 2011

This class stores information about inelastic scattering channels for collisons in the CB3D routine. It contains the resonance information for two exit particles, and is designed to be referenced by the net charge, strangeness, and baryon number of the combined particles.
*/
class CInelasticInfo{
public:
	//!Constructor
	CInelasticInfo(CResInfo *resinfo_1, CResInfo *resinfo_2, int type);
	//!A CResInfo pointer 
	/*!
	A CResInfo pointer that points to resonance information for one of the exit particles. By convention,
	the resonance with the smaller pole mass is stored in resinfo_1
	*/
	CResInfo *resinfo_1;
	CResInfo *resinfo_2;	/*!< A CResInfo pointer. \sa resinfo_1 */
	//!Collision type (Depreciated)
	/*!
	An integer that can be used to determine the collision type. Right now, the model has no method of differentiating
	between exit channels for which more data/certainty exists, and the type integer can be used to solve this. However,
	this has yet to be implemented.
	*/
	int type;
	//!The minimum invariant mass required to implent exit channel
	/*!
	This is the minimum invariant mass required for the process to occur. The collision routine compares this to
	\f$\sqrt{s}\f$ when considering whether to include the exit channel.
	*/
	double min_mass;
	int net_q;	/*!< the net charge of the combined particles */
	int net_s;	/*!< the net strangeness of the combined particles */
	int net_b;	/*!< the net baryon number of the combined particles */
	//!Print information to screen
	/*!
	A default print routine which allows for detailed information to be printed to the screen. Mostly used for debugging.
	*/
	void Print();
};
//!A class for handling the aggregated CInelasticInfo objects
/*!
\author Kevin Novak
\version 1.0
\date March 2011

This is a helper class that handles the CInelasticInfo objects. The class provides a list structure similar to the CResList class, with similar functionality. It handles reading in, organizing, and storing the various CInelasticInfo objects.
*/
class CInelasticList{
public:
	//!Filename where inelastic data is stored.
	string filename;
	//!Number of resonances in CB3D class. Repetitive, but simpler.
	int NResonances;
	//!Constructor
	CInelasticList();
	//!CResInfo pointer that points to the first CResInfo object in the list.
	CResInfo *GfirstResInfoptr;
	//!Reads in Inelastic information.
	/*!
	This method generates the inelastic information list. It will either read in the information to generate the inelastic arrays from /parameters/inelastic.dat (if FromFile is true), or create a new array (if FromFile is false) and store it to /parameters/inelastic.dat. This is to improve model function time, as the creation of the array (especially if UseInelasticArray = true) can be quite computationally intensive.
	\param[in] FromFile whether to read in inelastic information from a file, or create new information and store to a file.
	*/
	void ReadInelasticInfo(bool FromFile);
	//!Adds a CInelasticInfo object to the inelastic list.
	/*!
	By convention (and simplicity), both the CResInfo and CInelasticInfo lists are ordered in terms of increasing energy. For that reason, it is necessary to include a sorting method (SortedAdd). This method inserts a CInelasticInfo object into a CInelasticInfo list (either a element of ThermalArray or InelasticArray) at the appropriate minimum required energy.
	\sa CInelasticInfo::min_mass
	\param[in] list_in A C++ list container of CInelasticInfo objects.
	\param[in] inelastic_in the CInelasticInfo object to be added to the list.
	*/
	void SortedAdd(list<CInelasticInfo> &list_in, CInelasticInfo inelastic_in);
	//!Checks to see if a inelastic exit channel is valid.
	/*!
	This method is used for checking whether a inelastic exit channel of \f$1,2 \leftrightarrow 3,4\f$ is to be included in the InelasticArray.Besides checking that for conserved quantities (charge, strangeness, baryon number, and G parity, where applicable), this allows for more detailed selection of exit channels to be included.
	\param[in] res1 A CResInfo object corresponding to one incoming particle.
	\param[in] res2 A CResInfo object corresponding to the second incoming particle.
	\param[in] res3 A CResInfo object corresponding to the first outgoing particle.
	\param[in] res4 A CResInfo object corresponding to the second outgoing particle.
	*/
	bool AddToArrayCheck(CResInfo res1, CResInfo res2, CResInfo res3, CResInfo res4);
	parameterMap *parmap;	/*!< The parameter map for the given CB3D run. */
	//!A more detailed array of CInelasticInfo objects, designed to be referenced by CResInfo resonance number.
	/*!
	This is an optional array, dynamically allocated to be a NResonances by NResonances array of C++ list containers of CInelasticInfo objects. Similar to the MergeArray array in the CResList class, this is designed to be referenced by the resonance numbers of the incoming particles, and should (in theory) allow for completely individualized defining of exit channels. However, unlike the MergeArray, the sheer number of CInelasticInfo objects in the array mean that it is extremely computationally intensive to create and store the array. As such, this method is basically depreciated, and will be addressed in future revisions.
	
	\sa CResList::MergeArray
	*/
	list<CInelasticInfo>** InelasticArray;
	//!Themal Inelastic Scattering array.
	/*!
	This array stores all the inelastic scattering exit channels. The indices are in referenced in alphabetical order; \f$|b|\f$ (baryon number), \f$|q|\f$ (charge), \f$|s|\f$ (strangeness), \f$Sign(b)\f$ baryon sign, \f$Sign(q)\f$ charge sign, and \f$Sign(s)\f$ strangeness sign.
	*/
	list<CInelasticInfo> ThermalArray[3][5][7][2][2][2];
	static CB3D *b3d;	/*!< A pointer to the CB3D object the list belongs to. */
	
	static bool UseFile;	/*!< Determines whether or not to read in info from file. Set in CB3D constructor from parmap. */
	static bool UseInelasticArray;	/*!< Determines whether or not to define InelasticArray. Set in CB3D constructor from parmap. */
};

#endif
