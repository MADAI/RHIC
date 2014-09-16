#ifndef __OCTMESH_h_INCLUDE__
#define __OCTMESH_h_INCLUDE__
/*  Cmesh.h
 *  isHydro3
 *  Created by Joshua Vredevoogd on 2/11/09.
 */

#include "hydroDef.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <list>
#include <string.h>
#ifdef HYDRO_BOOST_THREADS
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/tss.hpp>
#endif
#include "coralutils.h"
#include "CCell.h"
#include "CEos.h"
#include "CMesh.h"
#include "cornelius.h"

class octMesh : public CMesh {
	
public:
	octMesh(parameterMap *pM);
	octMesh(parameterMap* pM, const char *fName);
	octMesh(CMesh* mesh);
	virtual ~octMesh();
	
	void klnInput();
	void fillActiveCells();
	void average(CMesh*,CMesh*); //puts average of the two meshs here
	void copyMesh(CMesh* m);
	void fixEntropyMesh();
	inline void initBGK() {std::cout << "BGK requires fullMesh." << std::endl; exit(1);}
	
	inline CCell* getCell(int eta, int x, int y){return mCells[eta+1][x+1][y+1];}	
	void update(int eta, int x, int y);
	void selfUpdate(int eta, int x, int y);
	void deaden();
	
	void flipEdge();
	void flipEdge(int,CCell*,CCell*,CCell*);
	
	void getFOS(double mFOS[XSIZE+YSIZE][2], double mFOSigma[XSIZE+YSIZE][2], int &foSize);
	bool getFOS(double mFOS[XSIZE+YSIZE][2], double mFOSigma[XSIZE+YSIZE][2], 
				double mFOVelo[XSIZE+YSIZE][2], double mFODiss[XSIZE+YSIZE][6], int &foSize); 
	bool getFOS(int eta, double mFOS[XSIZE+YSIZE][2], double mFOSigma[XSIZE+YSIZE][2], 
				double mFOVelo[XSIZE+YSIZE][2], double mFODiss[XSIZE+YSIZE][6], int &foSize);
	
	bool goodFOS(double mFOS[XSIZE+YSIZE][2], int foSize);
	int genFOS(CMesh* m);
	
	void checkAzimuthalSymmetry();
	void writeIEPvEta();
	
	double getdNd3p(double px, double py, double pz, double m);
	void calcdNdY(double m, double d, double *output, double dRap, int nRap);
	double getF( int x, int y, int n, double pt, double phi, double rap, double m);
	
};
#endif 