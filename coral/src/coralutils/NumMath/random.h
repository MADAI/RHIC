#ifndef __INCLUDE_CRANDOM_H__
#define __INCLUDE_CRANDOM_H__

#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
typedef double FourVector[4];

using namespace std;

class CRandom{
 public:
  //! returns from zero to 1
  double ran(void);
  //! returns integer from zero to imax-1
  unsigned long int iran(unsigned long int imax); 
  //! weights with exp(-x^2/2)
  double gauss(void);
  void gauss2(double *g1,double *g2);
  //! weights with exp(-x)
  double  ran_exp(void);
	CRandom(int seed);
	void reset(int seed);
	void generate_boltzmann(double mass,double T,FourVector &p);
	//void generate_boltzmann(double mass,double T,double *p);
	void generate_boltzmann_alt(double mass,double T,FourVector &p);
	//void generate_boltzmann_alt(double mass,double T,double *p);
	int GetNPoissonian(double eta);
private:
  gsl_rng *randy;
};

#endif
