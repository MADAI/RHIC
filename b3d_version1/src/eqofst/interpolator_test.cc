#ifndef __INCLUDE_INTERPOLATOR_Test_CC__
#define __INCLUDE_INTERPOLATOR_Test_CC__

/*
 *  interpolator_test.cc
 *  Created by Cyrus Faroughy on summer09.
 */

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#define __FIRST_STATIC_DEF__
#include <gsl/gsl_sf.h>
#include <eqofst.h>
#include <coral.h>
#include <interpolator.h>
#include<iostream>
using namespace std;

int main(){
//Declaration section:

double  e = 1873;         //value of energy density to interpolate
double  y_pres = 0;       //interpolated value of pressure
double y_temp = 0;        //interpolated value for temperature
double y_sdens = 0;       //interpolated value for sdens

double Cs2;              //speed of sound (1/c)
double tauIS_a;          //israel stewart relaxation time (fm/c)
double eta;              //shear viscosity (1/fm)
double B;                //bulk viscosity (1/fm)


//interpolator section:
CEos* interpolator;
interpolator=new CEos();
y_pres = interpolator->getP(e);    //get interpolated value for pressure
y_temp = interpolator->getT(e); //get interpolated value for temperature
y_sdens = interpolator->getS(e);    //get interpolated value for sdens

Cs2 = interpolator->getCs2(e); // returns speed of sound (1/c)
tauIS_a = interpolator->getTIS(e);  // returns israel stewart relaxation time (fm/c)
eta = interpolator->getSV(e);   // returns shear viscosity (1/fm)
B = interpolator->getBV(e);   // returns bulk viscosity (1/fm)

cout<<"================================================"<<endl;
cout<<"For an Energy density of "<<e<<" we get: "<<endl;
cout<<" "<<endl;
cout<<"Pressure = "<<y_pres<<endl;
cout<<"Temperature = "<<y_temp<<endl;
cout<<"Entropy density  = "<<y_sdens<<endl;



cout<<"Speed of Sound (1/c) = "<<Cs2<<endl;
cout<<"Relaxation time (fm/c) = "<<tauIS_a<<endl;
cout<<"Shear viscosity (1/fm) = "<<eta<<endl;
cout<<"Bulk viscosity (1/fm) = "<<B<<endl;
cout<<"================================================="<<endl;

return(0);
}
#endif
