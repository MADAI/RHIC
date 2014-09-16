#ifndef __INCLUDE_KERNEL_H
#define __INCLUDE_KERNEL_H

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <string>
#include "sf.h"
//#include "utils.h"
#include "constants.h"
#include "wavefunction.h"
#include "parametermap.h"
using namespace std;

string getKernelFilename( string datadir, int ell, double q );

/**
  *
  *
  */
class CKernel{

 public:
 
  CKernel( string kparsfilename="" );
  virtual ~CKernel();
  virtual double GetValue( int ell, double q, double r ) const;
  double GetValue( int ell, int iq, int ir ) const;
  bool Read( const parameterMap& parameters );
  bool Write( parameterMap& parameters );
  void ReadData( string datadirname );
  void WriteData( string datadirname );
	void WriteData_Paul( string datadirname );
  void Print();
  int GetLMAX();
  double GetDELR();
  double GetDELQ();
  int GetNQMAX();
  int GetNRMAX();
  double GetQ(int iq);
  virtual void Calc( CWaveFunction *wf );
	void Calc_Paul( CWaveFunction *wf );
  void Calc_ClassCoul( double ma, double mb, int zazb );
  void Calc_PureHBT();
  bool GetIDENTICAL();
  double GetPsiSquared( int iq, int ir, double ctheta );
  double GetPsiSquared( int iq, double r, double ctheta );
  double GetPsiSquared( double q, double r, double ctheta);

 private:

  bool IDENTICAL;
  int ellmax;
  int nrmax,nqmax;
  double delr,delq;
  double ***kernel;
  double *P;
  void ParsInit( string kparsfilename );
  double CalcPsiSquared( int iq, int ir, double ctheta );
  double CalcPsiSquared( int iq, double r, double ctheta );
  void CalcP( double ctheta );
};

class CKernelExactHBT: public CKernel {
  public:
  
    CKernelExactHBT( string kparsfilename="" ): CKernel( kparsfilename ){}
    double GetValue( int ell, double q, double r ) const 
        {return pow(-1.0,ell)*Bessel::jn(ell,2.0*q*r/HBARC);}
};


class CKernelWF{

 public:

  double GetPsiSquared(int iq,int ir,int ictheta);
  double GetPsiSquared(int iq,int ir,double ctheta);
  double GetPsiSquared(int iq,double r,double ctheta);
  double GetPsiSquared(double q,double r,double ctheta);
  void Calc(CWaveFunction *wf);
  void ReadData( string datadirname );
  void WriteData( string datadirname );
  double GetDELR();
  double GetDELQ();
  double GetDELCTHETA();
  int GetNQMAX();
  int GetNRMAX();
  int GetNCTHETA();
  bool GetIDENTICAL();
  CKernelWF( string kparsfilename );
  ~CKernelWF();

 private:

  bool IDENTICAL;
  int nctheta,nrmax,nqmax;
  double delr,delq,delctheta;
  double ***wfarray;
  void ParsInit( string kparsfilename );
  
};

#endif
