#ifndef __CWAVEFUNCTION_WF_PIPLUSPIMINUS_NOSTRONG_CC__
#define __CWAVEFUNCTION_WF_PIPLUSPIMINUS_NOSTRONG_CC__
#include "wavefunction.h"

CWaveFunction_pipluspiminus_nostrong::CWaveFunction_pipluspiminus_nostrong(string  parsfilename){
	int iq,ichannel,*I;
	double q;

	ParsInit(parsfilename);

	m1=MPI;
	m2=MPI;
	IDENTICAL=0;
	q1q2=-1;
	if(COULOMB==false) q1q2=0;
	nchannels=0;

	InitArrays();
	InitWaves();

	printf("pipluspiminus_nostrong wf initialized\n");
}

double CWaveFunction_pipluspiminus_nostrong::CalcPsiSquared(int iq,double r,double ctheta){
	double psisquared;
	double delta_s,delta_p,delta_d;
	complex<double> psi0;

	if(iq>=nqmax){
		psisquared=1.0;
	}
	else{
		psi0=planewave[iq]->planewave(r,ctheta);
		psisquared=real(psi0*conj(psi0));
		//psisquared*=RelativisticCorrection(r,iq);
	}
	return psisquared;
}

#endif
