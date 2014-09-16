#ifndef __INCLUDE_INTERPOLATOR_H__
#define __INCLUDE_INTERPOLATOR_H__
#include "eqofst.h"

/*
 *  interpolator.h
 *  Created by Cyrus Faroughy on summer09.
 */

class CEos{
    
 private:
  double step ;
  double e_max ;
  int N ; 
  
  double* e ;
  double* p ;
  double* t ;
  double* sdens;              
  double* tauIS_a;
  double* B;
  double* eta;
  double* dpde;
  double* dedT;
  double* Cs2;

 public:
  CEos();
  ~CEos(); 
  double getP(double e);    // get interpolated value of P (Mev) 
  double getT(double e);    // get interpolated value of T (Mev)           
  double getS(double e);    // get interpolated value of sdens

  double getCs2(double e);  // returns speed of sound (1/c) using Cs2 = sdens/dedT
  double getTIS(double e);  // returns israel stewart relaxation time (fm/c)
  double getSV(double e);   // returns shear viscosity (1/fm)
  double getBV(double e);   // returns bulk viscosity (1/fm)

}; 
#endif
