#ifndef __INTRINSIC_INCLUDE_CC__
#define __INTRINSIC_INCLUDE_CC__
#include "intrinsic.h"

CIntrinsic::CIntrinsic(CSpecies *species_set){
  int iQ,ispecies;
	species=species_set;
	density=new double[16];
	flux=new double[16];
	rho=new double[5];
	mu=new double[5];
	for(iQ=0;iQ<5;iQ++){
		rho[iQ]=0.0;
		mu[iQ]=0.0;
	}
	for(ispecies=0;ispecies<16;ispecies++){
		density[ispecies]=flux[ispecies]=0.0;
	}
}

void CIntrinsic::Reset(CSpecies *species_set){
  int iQ,ispecies,oldNcharges,oldNspecies;
	oldNcharges=species->Ncharges;
	oldNspecies=species->Nspecies;
  species=species_set;
	if(oldNcharges!=species->Ncharges){
		//delete [] rho;
		//delete [] mu;
		//rho=new double[species->Ncharges];
		//mu=new double[species->Ncharges];
		for(iQ=0;iQ<species->Ncharges;iQ++) rho[iQ]=mu[iQ]=0.0;
	}
	if(oldNspecies!=species->Nspecies){
		//delete [] density;
		//density=new double[species->N
		for(ispecies=0;ispecies<species->Nspecies;ispecies++) density[ispecies]=flux[ispecies]=0.0;
	}
}

void CIntrinsic::CopyFrom(CIntrinsic *intr){
	if(species!=intr->species) Reset(intr->species);
	int iQ,ispecies;
	for(iQ=0;iQ<species->Ncharges;iQ++){
		rho[iQ]=intr->rho[iQ];
		mu[iQ]=intr->mu[iQ];
	}
	for(ispecies=0;ispecies<species->Nspecies;ispecies++){
		density[ispecies]=intr->density[ispecies];
		flux[ispecies]=intr->flux[ispecies];
	}
	epsilon=intr->epsilon;
	T=intr->T;
	P=intr->P;
	dedT=intr->dedT;
	sdens=intr->sdens;
	eta=intr->eta;
	B=intr->B;
	
	tauIS_a=intr->tauIS_a;
	tauIS_b=intr->tauIS_b;
}

void CIntrinsic::Print(){
  int iQ,i;
	double denstot=0.0;
  printf("________ INTRINSIC PROPERTIES _________\n");
  printf("T=%g, epsilon=%g, P=%g, sdens=%g, dedT=%g\n",
	 T,epsilon,P,sdens,dedT);
  printf("Dens : ");
  for(i=0;i<species->Nspecies;i++){
		printf("%9.3e ",density[i]);
		denstot+=density[i];
	}
  printf("\ndenstot=%g\n",denstot);
	printf("Flux : ");
  for(i=0;i<species->Nspecies;i++) printf("%9.3e ",flux[i]);
  printf("\n");
  printf("Rho  : ");
  for(iQ=0;iQ<species->Ncharges;iQ++) printf("%10.3e ",rho[iQ]);
  printf("\n");
  printf("Mu   : ");
  for(iQ=0;iQ<species->Ncharges;iQ++) printf("%10.3e ",mu[iQ]);
  printf("\n_____________________________________\n");
}

void CIntrinsic::ZeroMu(){
  int iQ,Ncharges=species->Ncharges;
  for(iQ=0;iQ<Ncharges;iQ++) mu[iQ]=0.0;
}

CIntrinsic::~CIntrinsic(){
  if(rho!=NULL) delete [] rho;
	rho=NULL;
  if(density!=NULL) delete [] density;
	density=NULL;
	if(flux!=NULL) delete [] flux;
	flux=NULL;
  if(mu!=NULL) delete [] mu;
	mu=NULL;
}





#endif
