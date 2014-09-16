#ifndef __INCLUDE_EQOFST_ThreePhase_CC__
#define __INCLUDE_EQOFST_ThreePhase_CC__
#include "eqofst.h"

CEqofst_ThreePhase::CEqofst_ThreePhase(CSpecies *hspecies_set,string parsfilename_set){
	
	parsfilename=parsfilename_set;
	parameterMap eqofstpars;
	double kappa,dfact;
	parameter::ReadParsFromFile(eqofstpars,parsfilename.c_str());
  qgpspecies=new CSpecies_QGP();
  hspecies=hspecies_set;
  species=hspecies;
  MaxNcharges=hspecies->Ncharges;
  MakeNewtonArrays();
	
	c2qgp=parameter::getD(eqofstpars,"C2QGP",0.3);
	c2mixed=parameter::getD(eqofstpars,"C2MIXED",0.05);
	L=parameter::getD(eqofstpars,"L",500.0);
	T_h=parameter::getD(eqofstpars,"TC",165.0);
	etafact=parameter::getD(eqofstpars,"ETAFACT",0.2);
	Bfact=parameter::getD(eqofstpars,"BFACT",0.0);
	
  intrinsicTc=new CIntrinsic(hspecies);
  intrinsicTc->T=T_h;
  for(int iQ=0;iQ<hspecies->Ncharges;iQ++) intrinsicTc->mu[iQ]=0.0;
  FreeGasCalc_of_TMu(intrinsicTc);
  epsilon_h=intrinsicTc->epsilon;
  P_h=intrinsicTc->P;
  h_h=P_h+epsilon_h;
  sdens_h=h_h/T_h;
	sigma2_h=intrinsicTc->asigma*intrinsicTc->asigma;

	double denstot=0.0;
	for(int ispecies=0;ispecies<intrinsicTc->species->Nspecies;ispecies++) denstot+=intrinsicTc->density[ispecies];
	double tauIS_a=1.0/(denstot*2.5);  // assumes 25 mb cross section
	eta_h=pow(intrinsicTc->asigma,2)*tauIS_a/T_h;
	
	epsilon_qgp=epsilon_h+L;
	P_qgp=P_h+c2mixed*L;
	h_qgp=epsilon_qgp+P_qgp;
	dfact=1.0/(1.0+c2mixed);
	sdens_qgp=sdens_h*pow(h_qgp/h_h,dfact);
  T_qgp=h_qgp/sdens_qgp;
	sigma2_qgp=T_qgp*T_qgp*sdens_qgp/5.0;

	eta_qgp=etafact*sdens_qgp*HBARC/(4.0*PI);
	
}

void CEqofst_ThreePhase::Calc_of_TMu(CIntrinsic *intr){
  double h;
  if(intr->T>T_qgp){
    if(intr->species!=qgpspecies) intr->Reset(qgpspecies);
    h=(P_qgp+epsilon_qgp)*pow(intr->T/T_qgp,1.0+(1.0/c2qgp));
    intr->P=(P_qgp+c2qgp*(h-epsilon_qgp))/(1.0+c2qgp);
    intr->epsilon=h-intr->P;
    intr->sdens=h/intr->T;
		intr->dedT=h/(intr->T*c2qgp);
  }
	else if(intr->T>T_h){
		if(intr->species!=qgpspecies) intr->Reset(qgpspecies);
    h=(P_h+epsilon_h)*pow(intr->T/T_h,1.0+(1.0/c2mixed));
    intr->P=(P_h+c2mixed*(h-epsilon_h))/(1.0+c2mixed);
    intr->epsilon=h-intr->P;
    intr->sdens=h/intr->T;
		intr->dedT=h/(intr->T*c2mixed);
	}
  else{
    if(intr->species!=hspecies){
      intr->Reset(hspecies);
      intr->ZeroMu();
			printf("why am i here?\n");
    }
    FreeGasCalc_of_TMu(intr);
		printf("check check: T=%g, sigma2/PT=%g\n",intr->T,intr->asigma*intr->asigma/(intr->P*intr->T));
  }
	CalcVisc(intr);
}

void CEqofst_ThreePhase::Calc_of_EMu(CIntrinsic *intr){
  double dfact,kappa,h;
  if(intr->epsilon>epsilon_qgp){
		if(intr->species!=qgpspecies) intr->Reset(qgpspecies);
		dfact=1.0/(1.0+c2qgp);
		intr->P=P_qgp+c2qgp*(intr->epsilon-epsilon_qgp);
		h=intr->epsilon+intr->P;
		intr->sdens=sdens_qgp*pow(h/h_qgp,dfact);
		intr->T=h/intr->sdens;
	}
	else if(intr->epsilon>epsilon_h){
		if(intr->species!=qgpspecies) intr->Reset(qgpspecies);
		dfact=1.0/(1.0+c2mixed);
		intr->P=P_h+c2mixed*(intr->epsilon-epsilon_h);
		h=intr->epsilon+intr->P;
		intr->sdens=sdens_h*pow(h/h_h,dfact);
		intr->T=h/intr->sdens;
	}
	else{
    if(intr->species!=hspecies){
      intr->Reset(hspecies);
      intr->ZeroMu();
    }
    FreeGasCalc_of_EMu(intr);
  }
	CalcVisc(intr);
}

void CEqofst_ThreePhase::Calc_of_ERho(CIntrinsic *intr){
  double dfact,kappa,h;
	int iQ;
  if(intr->epsilon>epsilon_qgp){
		if(intr->species!=qgpspecies) intr->Reset(qgpspecies);
		dfact=1.0/(1.0+c2qgp);
		intr->P=P_qgp+c2qgp*(intr->epsilon-epsilon_qgp);
		h=intr->epsilon+intr->P;
		intr->sdens=sdens_qgp*pow(h/h_qgp,dfact);
		intr->T=h/intr->sdens;
	}
	else if(intr->epsilon>epsilon_h){
		if(intr->species!=qgpspecies) intr->Reset(qgpspecies);
		dfact=1.0/(1.0+c2mixed);
		intr->P=P_h+c2mixed*(intr->epsilon-epsilon_h);
		h=intr->epsilon+intr->P;
		intr->sdens=sdens_h*pow(h/h_h,dfact);
		intr->T=h/intr->sdens;
	}
  else{
    if(intr->species!=hspecies){
			intr->Reset(hspecies);
			for(iQ=0;iQ<hspecies->Ncharges;iQ++) intr->rho[iQ]=intrinsicTc->rho[iQ]*intr->epsilon/epsilon_h;
		}
    FreeGasCalc_of_ERho(intr);
  }
	CalcVisc(intr);
}

void CEqofst_ThreePhase::InitIntrinsic_EMu0(CIntrinsic *&intr,double epsilon){
  int iQ,i;
	double emin,xfactor;
  if(intr==NULL){
		intr=new CIntrinsic(qgpspecies);
		intr->T=200.0;
	}
	intr->epsilon=epsilon;
  if(epsilon>epsilon_h){
		if(intr->species!=qgpspecies) intr->Reset(qgpspecies);
    Calc_of_EMu(intr);
  }
  else{
		if(intr->species!=hspecies) intr->Reset(hspecies);

		for(iQ=0;iQ<hspecies->Ncharges;iQ++){
			intr->rho[iQ]=intrinsicTc->rho[iQ]*epsilon/epsilon_h;
		}
		Calc_of_ERho(intr);
  }
	CalcVisc(intr);
	
}

void CEqofst_ThreePhase::CalcVisc(CIntrinsic *intr){
	double w,Bmin=1.0E-4,Bmax,denstot;
	double deltaPmax_a=intr->P;
	double deltaPmax_b=intr->P;

	if(etafact>0.999){
		if(intr->epsilon<epsilon_qgp && intr->epsilon>epsilon_h){
			//intr->eta=etafact*intr->sdens*HBARC/(4.0*PI);
			w=(intr->epsilon-epsilon_h)/(epsilon_qgp-epsilon_h);
			intr->eta=w*eta_qgp+(1.0-w)*eta_h;
			intr->asigma=w*sigma2_qgp+(1.0-w)*sigma2_h;
			intr->asigma=sqrt(intr->asigma);		

			intr->tauIS_a=intr->eta*intr->T/(intr->asigma*intr->asigma);
			intr->amax=deltaPmax_a*0.5*sqrt(3.0);

			Bmax=Bfact*intr->sdens*HBARC/(4.0*PI);
			w=fabs(intr->epsilon-0.5*(epsilon_qgp+epsilon_h))/(epsilon_qgp-epsilon_h);
			intr->B=(1.0-w)*Bmax;
			if(intr->B<Bmin) intr->B=Bmin;
			intr->tauIS_b=HBARC*(1.25/PI)/intr->T;
			intr->bsigma=sqrt(intr->B*intr->T/intr->tauIS_b);
			intr->bmax=deltaPmax_b;
		}
		else if(intr->epsilon>=epsilon_qgp){
			intr->eta=etafact*intr->sdens*HBARC/(4.0*PI);
			intr->asigma=0.2*intr->T*intr->T*intr->sdens;
			intr->asigma=sqrt(intr->asigma);
			intr->tauIS_a=intr->eta*intr->T/(intr->asigma*intr->asigma);
			intr->amax=deltaPmax_a*0.5*sqrt(3.0);

			intr->B=Bmin;
			intr->tauIS_b=HBARC*(1.25/PI)/intr->T;
			intr->bsigma=sqrt(intr->B*intr->T/intr->tauIS_b);
			intr->bmax=deltaPmax_b;		
		}
		else{
			denstot=0.0;
			for(int ispecies=0;ispecies<intr->species->Nspecies;ispecies++) denstot+=intr->density[ispecies];
			intr->tauIS_a=1.0/(denstot*2.5);  // assumes 25 mb cross section
			intr->eta=(intr->asigma*intr->asigma)*intr->tauIS_a/intr->T;
			//if(intr->eta/intr->P>10.0){
				//intr->eta=10.0*intr->P;
				//intr->tauIS_a=10.0;
			//}
			//intr->asigma=sqrt(intr->eta*intr->T/intr->tauIS_a);
			intr->amax=deltaPmax_a*0.5*sqrt(3.0);

			intr->B=Bmin;
			intr->tauIS_b=0.5;
			intr->bsigma=sqrt(intr->B*intr->T/intr->tauIS_b);
			intr->bmax=deltaPmax_b;
		}
	}
	else{
		intr->eta=etafact*intr->sdens*HBARC/(4.0*PI);
		intr->tauIS_a=(1.25/PI)*HBARC/intr->T;
		intr->asigma=sqrt(intr->eta*intr->T/intr->tauIS_a);
		intr->amax=deltaPmax_a*0.5*sqrt(3.0);

		intr->B=Bmin;
		intr->tauIS_b=(1.25/PI)*HBARC/intr->T;
		intr->bsigma=sqrt(intr->B*intr->T/intr->tauIS_b);
		intr->bmax=deltaPmax_b;
	}
}
	
/*void CEqofst_ThreePhase::CalcVisc(CIntrinsic *intr){
	double w,Bmin=1.0E-4,Bmax;
	double deltaPmax_a=intr->P;
	double deltaPmax_b=intr->P;

	intr->B=0.0;
	if(intr->epsilon<epsilon_qgp && intr->epsilon>epsilon_h){
		Bmax=Bfact*intr->sdens*HBARC/(4.0*PI);
		w=fabs(intr->epsilon-0.5*(epsilon_qgp+epsilon_h))/(epsilon_qgp-epsilon_h);
		intr->B=(1.0-w)*Bmax;
	}
	if(intr->B<Bmin) intr->B=Bmin;

	intr->eta=etafact*intr->sdens*HBARC*(intr->T/T_h)/(4.0*PI);
	intr->tauIS_a=HBARC/intr->T;
	intr->tauIS_b=HBARC/intr->T;
		//if(intr->epsilon<epsilon_qgp&& intr->epsilon>epsilon_h)	intr->tauIS_b=4.0;

	intr->asigma=sqrt(intr->eta*intr->T/intr->tauIS_a);
	intr->bsigma=sqrt(intr->B*intr->T/intr->tauIS_b);

	intr->bmax=deltaPmax_b;
	intr->amax=deltaPmax_a*0.5*sqrt(3.0);

}*/


#endif
