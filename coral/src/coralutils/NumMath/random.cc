#ifndef __INCLUDE_CRANDOM_CC
#define __INCLUDE_CRANDOM_CC

#include "random.h"

using namespace std;

CRandom::CRandom(int seed){
  //Choose the algorithm, see gsldoc.ps (page 178)
  //randy=gsl_rng_alloc(gsl_rng_taus); // crappy but fast
  //randy=gsl_rng_alloc(gsl_rng_ranlxd2); // for anal-retentive types
  //randy=gsl_rng_alloc(gsl_rng_ran2); // Not really double precision numbers
  //randy=gsl_rng_alloc(gsl_rng_knuthran2); // sort of slow
  randy=gsl_rng_alloc(gsl_rng_ranlxd1); // Just right

  gsl_rng_set(randy,seed);
}

void CRandom::reset(int seed){
  gsl_rng_set(randy,seed);
}

double CRandom::ran(void){
  return gsl_rng_uniform(randy);
}

long unsigned int CRandom::iran(unsigned long int imax){
  return gsl_rng_uniform_int(randy,imax);
}

double CRandom::gauss(void){
  return gsl_ran_ugaussian(randy);
}

double CRandom::ran_exp(void){
  return -log( ran() );
}

void CRandom::gauss2(double *randy1,double *randy2){
  double x,y,r2,r,c,s;
TRY_AGAIN:
  x=1.0-2.0*gsl_rng_uniform(randy);
  y=1.0-2.0*gsl_rng_uniform(randy);
  r2=x*x+y*y;
  if(r2>1.0) goto TRY_AGAIN;
  r=sqrt(r2);
  c=x/r;
  s=y/r;
  *randy1=c*sqrt(-2.0*log(r2));
  *randy2=(s/c)**randy1;
}

/*

void CRandom::generate_boltzmann(double mass,double T,double *p){
  const double PI=4.0*atan(1.0);
  double r1,r2,r3,a,b,c;
  double pmag,ctheta,stheta,phi,pgauss;
  if(T/mass>0.6){
  GB_TRYAGAIN:
    r1=ran();
    r2=ran();
    r3=ran();
		a=-log(r1); b=-log(r2); c=-log(r3);
		pmag=T*(a+b+c);
    p[0]=sqrt(pmag*pmag+mass*mass);
    if(ran()>exp((pmag-p[0])/T)) goto GB_TRYAGAIN;
    ctheta=(a-b)/(a+b);
    stheta=sqrt(1.0-ctheta*ctheta);
    phi=T*T*pow(a+b,2)/(pmag*pmag);
    phi=2.0*PI*phi;
    p[3]=pmag*ctheta;
    p[1]=pmag*stheta*cos(phi);
    p[2]=pmag*stheta*sin(phi);
  }
  else generate_boltzmann_alt(mass,T,p);
}

void CRandom::generate_boltzmann_alt(double mass,double T,double *p){
  const double PI=4.0*atan(1.0);
  double r1,r2,r3,r0,I1,I2,I3,Itot;
  double pmag,E,ctheta,stheta,phi,pgauss,K;
GB_TRYAGAIN:
	r0=ran();
	I1=mass*mass;
	I2=2.0*mass*T;
	I3=2.0*T*T;
	Itot=I1+I2+I3;
	if(r0<I1/Itot){
		r1=ran();
		K=-T*log(r1);
	}
	else if(r0<(I1+I2)/Itot){
		r1=ran();
		r2=ran();
		K=-T*log(r1*r2);
	}
	else{
		r1=ran();
		r2=ran();
		r3=ran();
		K=-T*log(r1*r2*r3);
	}
	E=K+mass;
	pmag=sqrt(E*E-mass*mass);
	r0=ran();
	if(r0>pmag/E) goto GB_TRYAGAIN;
	phi=2.0*PI*ran();
	ctheta=1.0-2.0*ran();
	stheta=sqrt(1.0-ctheta*ctheta);
	p[3]=pmag*ctheta;
	p[1]=pmag*stheta*cos(phi);
	p[2]=pmag*stheta*sin(phi);
	p[0]=E;
}

*/

void CRandom::generate_boltzmann(double mass,double T,FourVector &p){
  const double PI=4.0*atan(1.0);
  double r1,r2,r3,a,b,c;
  double pmag,ctheta,stheta,phi,pgauss;
  if(T/mass>0.6){
  GB_TRYAGAIN:
    r1=ran();
    r2=ran();
    r3=ran();
		a=-log(r1); b=-log(r2); c=-log(r3);
		pmag=T*(a+b+c);
    p[0]=sqrt(pmag*pmag+mass*mass);
    if(ran()>exp((pmag-p[0])/T)) goto GB_TRYAGAIN;
    ctheta=(a-b)/(a+b);
    stheta=sqrt(1.0-ctheta*ctheta);
    phi=T*T*pow(a+b,2)/(pmag*pmag);
    phi=2.0*PI*phi;
    p[3]=pmag*ctheta;
    p[1]=pmag*stheta*cos(phi);
    p[2]=pmag*stheta*sin(phi);
  }
  else generate_boltzmann_alt(mass,T,p);
}

void CRandom::generate_boltzmann_alt(double mass,double T,FourVector &p){
  const double PI=4.0*atan(1.0);
  double r1,r2,r3,r0,I1,I2,I3,Itot;
  double pmag,E,ctheta,stheta,phi,pgauss,K;
GB_TRYAGAIN:
	r0=ran();
	I1=mass*mass;
	I2=2.0*mass*T;
	I3=2.0*T*T;
	Itot=I1+I2+I3;
	if(r0<I1/Itot){
		r1=ran();
		K=-T*log(r1);
	}
	else if(r0<(I1+I2)/Itot){
		r1=ran();
		r2=ran();
		K=-T*log(r1*r2);
	}
	else{
		r1=ran();
		r2=ran();
		r3=ran();
		K=-T*log(r1*r2*r3);
	}
	E=K+mass;
	pmag=sqrt(E*E-mass*mass);
	r0=ran();
	if(r0>pmag/E) goto GB_TRYAGAIN;
	phi=2.0*PI*ran();
	ctheta=1.0-2.0*ran();
	stheta=sqrt(1.0-ctheta*ctheta);
	p[3]=pmag*ctheta;
	p[1]=pmag*stheta*cos(phi);
	p[2]=pmag*stheta*sin(phi);
	p[0]=E;
}

int CRandom::GetNPoissonian(double eta){
  return gsl_ran_poisson(randy,eta);
}

#endif
