#ifndef __INCLUDE_SFIT_MINUIT_BLAST_CC__
#define __INCLUDE_SFIT_MINUIT_BLAST_CC__

CCF2S_Minuit_Blast::CCF2S_Minuit_Blast(CSourceCalc *scset,C3DArray *cexpset,
				       C3DArray *cerrorset,C3DArray *ctheory3Dset,
				       CCHArray *ctheoryset,CCHArray *sourceset,
				       CKernel *kernelset){
  ndim=3;

  // npars is different for different subclasses
  npars=4;
  //
  sourcecalc=scset;
  cexp3D=cexpset;
  cerror3D=cerrorset;
  ctheory=ctheoryset;
  ctheory3D=ctheory3Dset;
  source=sourceset;
  kernel=kernelset;
  if(pars!=NULL) delete [] pars;
  pars=new CMNPars[npars];
  if(xval!=NULL) delete [] xval;
  xval=new double[npars];

  // initialization of pars is also unique to given subclass
  pars[0].Set("lambda",parameter::getD(sourcecalc->spars,"lambda",1),1.0,0,2.0);
  pars[1].Set("R",parameter::getD(sourcecalc->spars,"R",13.0),1.0,0.0,20.0);
  pars[2].Set("Tau",parameter::getD(sourcecalc->spars,"Tau",12.0),1.0,0.0,20.0);
  pars[3].Set("DelTau",parameter::getD(sourcecalc->spars,"DelTau",5),1.0,0.0,20.0);
  InitMinuit();
  //FixPar(6);
  //FixPar(5);
  //FixPar(4);
}
#endif

