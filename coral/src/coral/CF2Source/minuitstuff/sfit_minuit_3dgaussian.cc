#ifndef __INCLUDE_SFIT_MINUIT_3DGAUSSIAN_CC__
#define __INCLUDE_SFIT_MINUIT_3DGAUSSIAN_CC__

CCF2S_Minuit_3DGaussian::CCF2S_Minuit_3DGaussian(CSourceCalc *scset,C3DArray *cexpset,
				       C3DArray *cerrorset,C3DArray *ctheory3Dset,
				       CCHArray *ctheoryset,CCHArray *sourceset,
				       CKernel *kernelset){
  ndim=3;

  // npars is different for different subclasses
  npars=5;
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
  pars[0].Set("lambda",parameter::getD(sourcecalc->spars,"lambda",1),1.0,0,0);
  pars[1].Set("Rx",parameter::getD(sourcecalc->spars,"Rx",5),1.0,0.0,20.0);
  pars[2].Set("Ry",parameter::getD(sourcecalc->spars,"Ry",5),1.0,0.0,20.0);
  pars[3].Set("Rz",parameter::getD(sourcecalc->spars,"Rz",5),1.0,0.0,20.0);
  pars[4].Set("Xoff",parameter::getD(sourcecalc->spars,"Xoff",0),1.0,0,0);
  //pars[5].Set("Yoff",0.0,1.0,0,0);
  //pars[6].Set("Zoff",0.0,1.0,0,0);

  InitMinuit();
  //FixPar(6);
  //FixPar(5);
  //FixPar(4);
}
#endif

