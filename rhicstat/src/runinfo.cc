#ifndef __RUNINFO_CC__
#define __RUNINFO_CC__
#include "rhicstat.h"
using namespace std;
CRunInfo::CRunInfo(int NX,int NY){
	x=new double[NX];
	w=new double[NX];
	y=new double[NY];
	z=new double[NY];
	zfiterror=new double[NY];
	sigmay=new double[NY];
	xlinear=new double[NX];
	ylinear=new double[NY];
	zlinear=new double[NY];
	zfit=new double[NY];
	yfit=new double[NY];
	good=true;
}

void CRunInfo::Print(){
	int ix,iz,iy;
	int NX=rhicstat->NX,NY=rhicstat->NY,NZ=rhicstat->NZ;
	double yreal,yfitreal;
	for(ix=1;ix<NX;ix++)
	printf("%10.3e =%s\n",rhicstat->xbar[ix]+x[ix]*(rhicstat->xmax[ix]-rhicstat->xmin[ix])/sqrt(12.0),rhicstat->xname[ix].c_str());
	printf(" i     z      zfit     zfiterror   zlinear\n");
	for(iz=0;iz<NZ;iz++){
		printf("%2d %8.5f %8.5f+/-%8.5f %8.5f\n",iz,z[iz],zfit[iz],zfiterror[ix],zlinear[iz]);
	}
	printf(" i     y         yfit     sigma_y  name\n");
	for(iy=0;iy<NY;iy++){
		yreal=rhicstat->ybar[iy]+y[iy]*rhicstat->sigmaybar[iy];
		yfitreal=rhicstat->ybar[iy]+yfit[iy]*rhicstat->sigmaybar[iy];
		printf("%2d %10.3e %10.3e %10.3e  %s\n",iy,yreal,yfitreal,rhicstat->sigmaybar[iy],rhicstat->yname[iy].c_str());
	}
	printf("----------------------------------------\n");
}

#endif

