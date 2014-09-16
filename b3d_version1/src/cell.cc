#ifndef __PART_CC__
#define __PART_CC__

#include "b3d.h"
using namespace std;

CB3D *CB3DCell::b3d=NULL;

CB3DCell::CB3DCell(double xminset,double xmaxset,double yminset,double ymaxset,double etaminset,double etamaxset){
	xmin=xminset; xmax=xmaxset; ymin=yminset; ymax=ymaxset; etamin=etaminset; etamax=etamaxset;	
	ireflection=0;
	creflection=NULL;
}

void CB3DCell::Print(){
	printf("___ CELL INFO _____\n");
	printf("ix=%d, iy=%d, ieta=%d, xmin=%g, xmax=%g, ymin=%g, ymax=%g, etamin=%g, etamax=%g\n", ix,iy,ieta,xmin,xmax,ymin,ymax,etamin,etamax);
	printf("%d parts in cell\n",int(partmap.size()));
	printf("---------------------\n");
}

#endif
