#ifndef __CELL_H__
#define __CELL_H__
#include "b3d.h"

using namespace std;

class CB3D; class CB3DCell; class CMuTInfo;

typedef unordered_map<long int,CPart *> CPartMap;
typedef multimap<double,CAction *> CActionMap;
typedef pair<long int,CPart*> CPartPair;
typedef pair<double,CAction*> CActionPair;

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
	vector<double> dens;
	
	CB3DCell(double xmin,double xmax,double ymin,double ymax,double etamin,double etamax);
	static CB3D *b3d;
};

//!A cell in the expanding cell mesh
/*!
\version 1.0
\author Scott Pratt
\date March 2011

In the CB3D model, the model space is expressed as a mesh grid of cells that expand as time propogates. The mesh is populated by cells, which are CB3DCell objects. This class keeps track of the particles populating it (in a particle map), as well as its spatial dimensions and neighbors. The neighbors are especially relevant, as actions such as collisions are scheduled by checking against particles inside its the current cell, as well as all neighboring cells.
*/

#endif