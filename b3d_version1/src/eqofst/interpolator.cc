#ifndef __INCLUDE_INTERPOLATOR_CC__
#define __INCLUDE_INTERPOLATOR_CC__

#include <cmath>
#include "interpolator.h"
#include <iostream>
#include <fstream>
#include "def.h"
using namespace std;

/*
 *  interpolator.cc
 *  Created by Cyrus Faroughy on summer09.
 */

//=====================================================================================


// Constructor that creates mesh to use in the interpolator function
CEos::CEos(){

double rhoovers;

//Chose stepsize and maximum energy density in mesh from def.h:
step = INTERP_STEP;
e_max = INTERP_EMAX;

int N = static_cast<int>(e_max/step);
string parsfilename="eqofstpars.dat";



// (1) CREATE OBJECTS:

//Select specie from preprocessing SPECIE in def.h                                   
SPECIE* cspec;
cspec = new SPECIE();

CIntrinsic* intr;
intr = new CIntrinsic(cspec);

CIntrinsic* intr_Tc;
intr_Tc = new CIntrinsic(cspec);
 
CEqofst_ThreePhase*  eqofst;
eqofst = new CEqofst_ThreePhase(cspec,parsfilename);
 
// (2) BEGIN BODY OF CONSTRUCTOR:

//..............Uncomment if you want to create an output file with all the mesh points to be................
 /*
ofstream outStream;
 outStream.open("output.dat");
 if(outStream.fail())                 //Tell the user if the file cannot be opened.
  {
    cout << "Output file opening failed!\n";
    exit(1);
  }
// */ //......................................................................................................


//Allocate arrays
e = new double [N+1];       //epsilon array
p = new double [N+1];       //pressure array
t = new double [N+1];       //temperature array
sdens = new double [N+1];   //entropy density array 
tauIS_a = new double [N+1]; //relaxation time array
B = new double [N+1];       //bulk viscosity array 
eta = new double [N+1];     //shear viscosity array
dpde = new double [N+1];    //seep of sound squared
dedT = new double [N+1];    //dedT
Cs2 = new double [N+1];     //using Cs2=sdens/dedT


//Initialize
e[0]=p[0]=t[0]=sdens[0]=Cs2[0]=0;
intr->T=100;      
intr->ZeroMu();
eqofst->FreeGasCalc_of_TMu(intr); 

//rhoovers stuff
intr_Tc->T=eqofst->T_h;
intr_Tc->ZeroMu();
eqofst->FreeGasCalc_of_TMu(intr_Tc);
intr->CopyFrom(intr_Tc);

//Construct arrays
for(int i=1, epsilon=10;epsilon<=e_max;epsilon+=step, i++){
intr->epsilon=epsilon;

//calc rhoovers 
if(intr->epsilon<=intr_Tc->epsilon){
	for(int iQ=0;iQ<intr->species->Ncharges;iQ++){
		intr->rho[iQ]=(intr_Tc->rho[iQ]/intr_Tc->sdens)*intr->sdens;
	}
	eqofst->Calc_of_ERhoOverS(intr);
}

//intr->ZeroMu();
eqofst->Calc_of_EMu(intr);
eqofst->Calc_of_TMu(intr);
e[i]=intr->epsilon;
p[i]=intr->P;
t[i]=intr->T;
sdens[i]=intr->sdens;
tauIS_a[i]=intr->tauIS_a;
B[i]=intr->B;
eta[i]=intr->eta;
dedT[i]= intr->dedT;
Cs2[i]= sdens[i]/dedT[i];



//printf("e=%g  P=%g  s=%g  T=%g  tau_a=%g  B=%g  eta=%g Cs2=%g dedT=%g\n",e[i],p[i],sdens[i],t[i],tauIS_a[i],B[i],eta[i],Cs2[i],dedT[i]);   
	}

//get dP/de
//for(int i=1, epsilon=10;epsilon<e_max;epsilon+=step, i++){
//dpde[i]=(p[i+1]-p[i])/(e[i+1]-e[i]);
//printf("e=%g  P=%g  s=%g  T=%g  tau_a=%g  B=%g  eta=%g dpde=%g Cs2=%g dedT%g\n",e[i],p[i],sdens[i],t[i],tauIS_a[i],B[i],eta[i],dpde[i],Cs2[i],dedT[i]); 
//}


//......................Also uncomment this if you want to create an output file with all the mesh points to plot..................................... 
  /* 
 for(int i=0; i < N; i++)
    {
    outStream << " " <<e[i]<<" "<<p[i]<<" "<<sdens[i]<<" "<<t[i]<<" "<<tauIS_a[i]<<" "<<B[i]<<" "<<eta[i]<<" "<<dpde[i]<<" "<<Cs2[i]<<" "<<dedT[i]<<endl;  //Output  
      }
    outStream << endl;
 outStream.close(); // */                                       //Close the outStream
//..........................................................................................................................................
}   

//============================================================================================================================================
//DESTRUCTOR TO DEALLOCATE ARRAYS:

CEos::~CEos(){		
delete[ ] e;
delete[ ] p;
delete[ ] t;
delete[ ] sdens;
delete[ ] tauIS_a;
delete[ ] B;
delete[ ] eta;
delete[ ] dpde;
delete[ ] dedT;
delete[ ] Cs2;

}

//=====================================================================================
//Interpolator function for pressure:
// X is the input energy density. The function returns an interpolated value for PRESSURE.
// Value of energy density in the interval [0, e_max].

double CEos::getP(double x) {

double y_pres; //interpolated value for pressure to return
int i_int; //value of i right below x

if(x<0) cout<<"epsilon is negative..."<<endl;
    
else if (x <= e_max){
 
// find value of i right below x:
i_int = static_cast<int>( floor(x/step) );   

double factor = (x - e[i_int]) / (e[i_int+1] - e[i_int]);
y_pres = p[i_int] + factor*(p[i_int+1] - p[i_int]);

}

else if (x>e_max) {
cout<<"In getP: epsilon>e_max not allowed"<<endl;
exit(1);
}

return(y_pres);               
}

//=====================================================================================
//Interpolator function for temperature:
// X is the input energy density. The function returns an interpolated value for TEMPERATURE.


double CEos::getT(double x) {

double y_temp;//interpolated value for temperature to return
int i_int; //value of i right below x

if(x<0) cout<<"epsilon is negative..."<<endl;
    
else if (x <= e_max){
 
// find value of i right below x:
i_int = static_cast<int>( floor(x/step) );   


double factor = (x - e[i_int]) / (e[i_int+1] - e[i_int]);
y_temp = t[i_int] + factor*(t[i_int+1] - t[i_int]);
               
}

else if (x>e_max) {
cout<<"In getT: epsilon>e_max not allowed"<<endl;
exit(1);
}

return(y_temp);               
}


//=====================================================================================
//Interpolator function for sdens:
// X is the input energy density. The function returns an interpolated value for sdens.


double CEos::getS(double x) {

double y_sdens;//interpolated value for temperature to return
int i_int; //value of i right below x

if(x<0) cout<<"epsilon is negative..."<<endl;
    
else if (x <= e_max){
 
// find value of i right below x:
i_int = static_cast<int>( floor(x/step) );   


double factor= (x - e[i_int]) / (e[i_int+1] - e[i_int]);
y_sdens = sdens[i_int] + factor*(sdens[i_int+1] - sdens[i_int]);
               
}

else if (x>e_max) {
cout<<"In getS: epsilon>e_max not allowed"<<endl;
exit(1);
}

return(y_sdens);               
}

//=====================================================================================
//Interpolator function for tauIS_a:
// X is the input energy density. The function returns an interpolated value for relaxation time tauIS_a.

double CEos::getTIS(double x) {

double y_tauIS;//interpolated value for tauIS_a to return
int i_int; //value of i right below x

if(x<0) cout<<"epsilon is negative..."<<endl;
    
else if (x <= e_max){
 
// find value of i right below x:
i_int = static_cast<int>( floor(x/step) );  


double factor= (x - e[i_int]) / (e[i_int+1] - e[i_int]);
y_tauIS = tauIS_a[i_int] + factor*(tauIS_a[i_int+1] - tauIS_a[i_int]);
               
}

else if (x>e_max) {
cout<<"In getTIS: epsilon>e_max not allowed"<<endl;
exit(1);
}

return(y_tauIS);
}

//=====================================================================================
//Interpolator function for B:
// X is the input energy density. The function returns an interpolated value for bulk viscosity B.

double CEos::getBV(double x) {

double y_B;//interpolated value for B to return
int i_int; //value of i right below x

if(x<0) cout<<"epsilon is negative..."<<endl;
    
else if (x <= e_max){
 
// find value of i right below x:
i_int = static_cast<int>( floor(x/step) );  


double factor= (x - e[i_int]) / (e[i_int+1] - e[i_int]);
y_B = B[i_int] + factor*(B[i_int+1] - B[i_int]);
               
}

else if (x>e_max) {
cout<<"In getBV: epsilon>e_max not allowed"<<endl;
exit(1);
}

return(y_B);
}

//=====================================================================================
//Interpolator function for eta:
// X is the input energy density. The function returns an interpolated value for shear viscosity eta.

double CEos::getSV(double x) {

double y_eta;//interpolated value for eta to be returned
int i_int; //value of i right below x

if(x<0) cout<<"epsilon is negative..."<<endl;
    
else if (x <= e_max){
 
// find value of i right below x:
i_int = static_cast<int>( floor(x/step) );   


double factor= (x - e[i_int]) / (e[i_int+1] - e[i_int]);
y_eta = eta[i_int] + factor*(eta[i_int+1] - eta[i_int]);
               
}

else if (x>e_max) {
cout<<"In getSV: epsilon>e_max not allowed"<<endl;
exit(1);
}

return(y_eta);
}


//=====================================================================================
//Interpolator function for Cs2:
// X is the input energy density. The function returns an interpolated value for Cs2.

double CEos::getCs2(double x) {

double y_Cs2;//interpolated value for Cs2 to be returned
int i_int; //value of i right below x

if(x<0) cout<<"epsilon is negative..."<<endl;
    
else if (x <= e_max){
 
// find value of i right below x:
i_int = static_cast<int>( floor(x/step) );   


double factor= (x - e[i_int]) / (e[i_int+1] - e[i_int]);
y_Cs2 = Cs2[i_int] + factor*(Cs2[i_int+1] - Cs2[i_int]);
               
}

else if (x>e_max) {
cout<<"In getCs2: epsilon>e_max not allowed"<<endl;
exit(1);
}

return(y_Cs2);
}  

#endif
