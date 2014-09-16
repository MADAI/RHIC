#ifndef __FULLMESH_h_INCLUDE__
#define __FULLMESH_h_INCLUDE__

#include "hydroDef.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <string.h>
#include <coralutils.h>
#include "CCell.h"
#include "CEos.h"
#include "CMesh.h"

class fullMesh : public CMesh {
public:
	fullMesh(parameterMap*);
	fullMesh(parameterMap* pM, const char *fName);
	fullMesh(CMesh* mesh);
	virtual ~fullMesh();	
	
	void klnInput();
	void bgkInput();
	void fillActiveCells();
	void average(CMesh*,CMesh*);  //puts average of the two meshs here
	void copyMesh(CMesh* m);
	void fixEntropyMesh();
	void initBGK();
	
	void update(int eta, int x, int y);
	void selfUpdate(int eta, int x, int y);
	inline CCell* getCell(int eta, int x, int y){return mCells[eta+mNSizeOrig+1][x+mXSizeOrig+1][y+mYSizeOrig+1];}
	void deaden();
	
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