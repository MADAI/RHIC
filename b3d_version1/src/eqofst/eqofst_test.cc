#include "eqofst.h"

using namespace std;

int main(){

	int iQ;
	double rhotot,sigma2,zeta;
	FILE *fptr=fopen("giorgio.dat","w");

	CEqofst_ThreePhase *eqofst;
	CSpecies_StandardHadrons5Q *hspecies=new CSpecies_StandardHadrons5Q();
	//CSpecies_StandardHadrons_Equil *hspecies=new CSpecies_StandardHadrons_Equil();
	//CSpecies_PionsOnly *hspecies=new CSpecies_PionsOnly();
	CIntrinsic *intr=new CIntrinsic(hspecies);
	CIntrinsic *intr_Tc=new CIntrinsic(hspecies);
	eqofst=new CEqofst_ThreePhase(hspecies,"eqofstpars.dat");
	intr_Tc->T=eqofst->T_h;
	intr_Tc->ZeroMu();
	eqofst->Calc_of_TMu(intr_Tc);
	intr->T=eqofst->T_h;
	intr->ZeroMu();
	eqofst->Calc_of_TMu(intr);

	/*
	for(intr->T=3;intr->T<20;intr->T+=3){
	//intr->T=0.2*139.57
	intr->ZeroMu();
	eqofst->Calc_of_TMu(intr);
	eqofst->GetSigmaB(intr);
	}
	printf("_____________________________________________________\n");
	*/

	printf("^^^^^^^^^ Ncharges=%d ^^^^^^^^^^^^^^^^^\n",hspecies->Ncharges);
	fprintf(fptr,"    T      epsilon     sdens      sigma2   sigma2/P*T   zeta   4*pi*zeta/hbar*s\n");

	for(intr->sdens=intr_Tc->sdens;intr->sdens>0.01*intr_Tc->sdens;intr->sdens=pow(10,-0.1)*intr->sdens){
		rhotot=0.0;
		for(iQ=0;iQ<hspecies->Ncharges;iQ++){
			intr->rho[iQ]=intr_Tc->rho[iQ]*intr->sdens/intr_Tc->sdens;
			rhotot+=intr->rho[iQ];
		}
		eqofst->Calc_of_SRho(intr);
		sigma2=eqofst->GetSigmaB(intr);
		zeta=(sigma2/intr->T)/(rhotot*2.0);
		fprintf(fptr,"%8.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
			intr->T,intr->epsilon,intr->sdens,sigma2,sigma2/(intr->P*intr->T),zeta,4.0*PI*zeta/(197.326*intr->sdens));
	}
	fclose(fptr);


	return 0;
}
