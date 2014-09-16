#ifndef __INCLUDE_S2C_H
#define __INCLUDE_S2C_H
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include "coral.h"

using namespace std;

class CSourceCalc;

namespace S2CF{
  void s2c(C3DArray *s,CWaveFunction *wf,C3DArray *cf);
  void s2c(C3DArray *s,CKernelWF *kernel,C3DArray *cf);
  void s2c(C3DArray *s,CWaveFunction *wf,C3DArray *cf);
  void s2c(CCHArray *s,CKernel *kernel,CCHArray *cf);
  void s2c(int Lx,int Ly,int Lz,CCHArray *s,CKernel *kernel,CCHArray *cf);
  void s2c(CMCList *lista,CMCList *listb,CWaveFunction *wf,C3DArray *cf);
  void s2c(CMCList *lista,CMCList *listb,CKernel *kernel,C3DArray *cf);
  void s2c(CMCList *lista,CMCList *listb,CKernelWF *kernel,C3DArray *cf,int NMC);
  void s2c(CMCList *lista,CMCList *listb,CKernelWF *kernel,C3DArray *cf3d);
	void s2c_gauss(CSourceCalc *sourcecalc,CKernelWF *kernel,C3DArray *cf3d);
	void s2c_bowlersinyukov(CSourceCalc *sourcecalc,CKernel *kernel,C3DArray *cf3d);
	void s2c_bosons(CMCList *list,C3DArray *cf3d);
};

#endif
