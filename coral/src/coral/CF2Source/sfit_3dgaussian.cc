#ifndef __INCLUDE_SFIT_3DGAUSSIAN_CC__
#define __INCLUDE_SFIT_3DGAUSSIAN_CC__
#include "sfit.h"

class C3DArray;
class CKernel;
class CSourceCalc;

CCF2SFit_3DGaussian::CCF2SFit_3DGaussian(CSourceCalc *scset,
					 C3DArray *cexpset,
					 C3DArray *cerrorset,
					 C3DArray *ctheory3Dset,
					 CCHArray *ctheoryset,
					 CCHArray *sourceset,
					 CKernel *kernelset){
  // npars is different for different subclasses
  //
  sourcecalc=scset;
  cexp3D=cexpset;
  cerror3D=cerrorset;
  ctheoryCH=ctheoryset;
  ctheory3D=ctheory3Dset;
  sourceCH=sourceset;
  kernel=kernelset;
  Init();

  // initialization of pars is also unique to given subclass

  AddPar("lambda",parameter::getD(sourcecalc->spars,"lambda",0.5),
	 0.02,0.0,1.5);
  AddPar("Rx",parameter::getD(sourcecalc->spars,"Rx",5),
	 0.2,1.0,12.0);
  AddPar("Ry",parameter::getD(sourcecalc->spars,"Ry",5),
	 0.2,1.0,12.0);
  AddPar("Rz",parameter::getD(sourcecalc->spars,"Rz",5),
	 0.2,1.0,12.0);
  AddPar("Xoff",parameter::getD(sourcecalc->spars,"Xoff",0.0),
	     0.2,-10.0,10.0);

  InitErrorMatrix();

}
#endif

