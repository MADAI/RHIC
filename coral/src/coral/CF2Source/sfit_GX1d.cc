#ifndef __INCLUDE_SFIT_GX1D_CC__
#define __INCLUDE_SFIT_GX1D_CC__
#include "sfit.h"

CCF2SFit_GX1D::CCF2SFit_GX1D(CSourceCalc *scset,
			     CCHArray *cexpset,
			     CCHArray *cerrorset,
			     CCHArray *ctheoryset,
			     CCHArray *sourceset,
			     CKernel *kernelset){
  int i,j;

  //
  sourcecalc=scset;
  cexpCH=cexpset;
  cerrorCH=cerrorset;
  ctheoryCH=ctheoryset;
  sourceCH=sourceset;
  kernel=kernelset;
  Init();

  // initialization of pars is also unique to given subclass

  AddPar("lambdaG",parameter::getD(sourcecalc->spars,"lambdaG",0.3),
	     0.02,0.0,1.5);
  AddPar("R",parameter::getD(sourcecalc->spars,"R",5),
	     0.2,1.0,12.0);
  AddPar("lambdaX",parameter::getD(sourcecalc->spars,"lambdaX",0.3),
	     0.02,0.0,1.5);
  AddPar("X",parameter::getD(sourcecalc->spars,"X",10.0),
	     0.4,1.0,25.0);
  AddPar("a",parameter::getD(sourcecalc->spars,"a",5.0),
	     0.2,1.0,20.0);

  for(i=0;i<nfreepars;i++){
    for(j=0;j<nfreepars;j++){
      StepMatrix[i][j]=0.0;
      ErrorMatrix[i][j]=0.0;
      if(i==j){
	StepMatrix[i][j]=par[i]->error;
	ErrorMatrix[i][j]=par[i]->error*par[i]->error;
      }
    }
  }

}
#endif

