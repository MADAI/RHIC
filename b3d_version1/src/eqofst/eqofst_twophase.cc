#ifndef __INCLUDE_EQOFST_TWO_PHASE_CC__
#define __INCLUDE_EQOFST_TWO_PHASE_CC__
#include "eqofst.h"
#include "species.h"

CEqofst_TwoPhase::CEqofst_TwoPhase(CSpecies *hspecies_set,double T_c_set,double c2qgp_set){
  qgpspecies=new CSpecies_QGP();
  hspecies=hspecies_set;
  species=hspecies;
  MaxNcharges=hspecies->Ncharges;
  MakeNewtonArrays();
  c2qgp=c2qgp_set;
  T_c=T_c_set;
	etafact=1.0/(4.0*PI);
	
  intrinsicTc=new CIntrinsic(species);
  intrinsicTc->T=T_c;
  for(int iQ=0;iQ<hspecies->Ncharges;iQ++) intrinsicTc->mu[iQ]=0.0;
  FreeGasCalc_of_TMu(intrinsicTc);
  epsilon_c=intrinsicTc->epsilon;
  P_c=intrinsicTc->P;
  h_c=P_c+epsilon_c;
  sdens_c=h_c/T_c;
  
}

void CEqofst_TwoPhase::Calc_of_TMu(CIntrinsic *intr){
  double h;
  if(intr->T>T_c){
    if(intr->species!=qgpspecies) intr->Reset(qgpspecies);
    h=(P_c+epsilon_c)*pow(intr->T/T_c,1.0+(1.0/c2qgp));
    intr->P=(P_c+c2qgp*(h-epsilon_c))/(1.0+c2qgp);
    intr->epsilon=h-intr->P;
    intr->sdens=h/intr->T;
		intr->dedT=h/(intr->T*c2qgp);
  }
  else{    
    if(intr->species!=hspecies){
      intr->Reset(hspecies);
      intr->ZeroMu();
    }
    FreeGasCalc_of_TMu(intr);
  }
}

void CEqofst_TwoPhase::Calc_of_EMu(CIntrinsic *intr){
  double dfact,kappa,h;
  if(intr->epsilon>epsilon_c){
    dfact=1.0/(1.0+c2qgp);
    kappa=P_c-c2qgp*epsilon_c;
    intr->P=c2qgp*intr->epsilon+kappa;
    h=intr->epsilon+intr->P;
    intr->sdens=sdens_c*pow(h/h_c,dfact);
    intr->T=h/intr->sdens;
    if(intr->species!=qgpspecies) intr->Reset(qgpspecies);
  }
  else{    
    if(intr->species!=hspecies){
      intr->Reset(hspecies);
      intr->ZeroMu();
    }
    FreeGasCalc_of_EMu(intr);
  }
}

void CEqofst_TwoPhase::Calc_of_ERho(CIntrinsic *intr){
  double dfact,kappa,h;
	int iQ;
  if(intr->epsilon>epsilon_c){
    if(intr->species!=qgpspecies) intr->Reset(qgpspecies);
    dfact=1.0/(1.0+c2qgp);
    kappa=P_c-c2qgp*epsilon_c;
    intr->P=c2qgp*intr->epsilon+kappa;
    h=intr->epsilon+intr->P;
    intr->sdens=sdens_c*pow(h/h_c,dfact);
    intr->T=h/intr->sdens;
  }
  else{
    if(intr->species!=hspecies){
      intr->Reset(hspecies);
			for(iQ=0;iQ<hspecies->Ncharges;iQ++)
      intr->rho[iQ]=intrinsicTc->rho[iQ]*intr->epsilon/epsilon_c;
    }
		FreeGasCalc_of_ERho(intr);
  }
}

void CEqofst_TwoPhase::Calc_of_SAlpha(CIntrinsic *intr){
  double dfact,kappa,h;
	int iQ;
  if(intr->sdens>sdens_c){
    kappa=P_c-c2qgp*epsilon_c;
    h=h_c*pow(intr->sdens/sdens_c,1.0+c2qgp);
    intr->T=h/intr->sdens;
    intr->epsilon=(h-kappa)/(1.0+c2qgp);
    intr->P=h-intr->epsilon;
    if(intr->species!=qgpspecies) intr->Reset(qgpspecies);
  }
  else{    
    if(intr->species!=hspecies){
      intr->Reset(hspecies);
			intr->ZeroMu();
    }
		FreeGasCalc_of_SAlpha(intr);
  }
}


void CEqofst_TwoPhase::Calc_of_SRho(CIntrinsic *intr){
  double dfact,kappa,h;
  if(intr->sdens>sdens_c){
    kappa=P_c-c2qgp*epsilon_c;
    h=h_c*pow(intr->sdens/sdens_c,1.0+c2qgp);
    intr->T=h/intr->sdens;
    intr->epsilon=(h-kappa)/(1.0+c2qgp);
    intr->P=h-intr->epsilon;
    if(intr->species!=qgpspecies) intr->Reset(qgpspecies);
  }
  else{    
    if(intr->species!=hspecies){
      intr->Reset(hspecies);
      intr->ZeroMu();
    }
    FreeGasCalc_of_SRho(intr);
  }
}

void CEqofst_TwoPhase::CalcVisc(CIntrinsic *intr){
  //intr->B=Bfact*intr->sdens;
	intr->B=0.0;
  intr->eta=etafact*intr->sdens;
	intr->tauIS=0.5;
}

void CEqofst_TwoPhase::InitIntrinsic_EMu0(CIntrinsic *&intrinsic,double epsilon){
  int iQ;
  if(intrinsic==NULL){
		intrinsic=new CIntrinsic(qgpspecies);
		intrinsic->T=250.0;
	}
	intrinsic->epsilon=epsilon;
  if(epsilon>epsilon_c){
		if(intrinsic->species!=qgpspecies) intrinsic->Reset(qgpspecies);
    for(iQ=0;iQ<qgpspecies->Ncharges;iQ++) intrinsic->mu[iQ]=0.0;
    Calc_of_EMu(intrinsic);
  }
  else{
    if(intrinsic->species!=hspecies) intrinsic->Reset(hspecies);
    for(iQ=0;iQ<hspecies->Ncharges;iQ++)
      intrinsic->rho[iQ]=intrinsicTc->rho[iQ]*epsilon/epsilon_c;
    Calc_of_ERho(intrinsic);
  }
	
}


#endif
