#ifndef __INCLUDE_SFIT_BLAST_CC__
#define __INCLUDE_SFIT_BLAST_CC__
#include "sfit.h"

CCF2SFit_Blast::CCF2SFit_Blast(CSourceCalc *scset,C3DArray *cexp3Dset,C3DArray *cerror3Dset, C3DArray *ctheory3Dset,CCHArray *ctheoryset,CCHArray *sourceset,CKernel *kernelset){
	int i,j;
	
	// nfreepars is different for different subclasses
	
	sourcecalc=scset;
	cexp3D=cexp3Dset;
	cerror3D=cerror3Dset;
	ctheoryCH=ctheoryset;
	ctheory3D=ctheory3Dset;
	sourceCH=sourceset;
	kernel=kernelset;
	Init();
	
	// initialization of pars is also unique to given subclass
	
	AddPar("lambda",parameter::getD(sourcecalc->spars,"lambda",1),
		0.01,0.0,1.5);
	AddPar("R",parameter::getD(sourcecalc->spars,"R",13),
		0.05,5.0,20.0);
	AddPar("Tau",parameter::getD(sourcecalc->spars,"Tau",12),
		0.05,5.0,20.0);
	AddPar("DelTau",parameter::getD(sourcecalc->spars,"DelTau",5),
		0.05,0.0,20.0);
	
	
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

