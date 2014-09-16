#ifndef __INCLUDE_EQOFST_CC__
#define __INCLUDE_EQOFST_CC__
#include "eqofst.h"

CEqofst::CEqofst(){
	rhotry=NULL;
	y=NULL;
	x=NULL;
	M=NULL;
}

CEqofst::CEqofst(CSpecies *species_set){
	species=species_set;
	MaxNcharges=species->Ncharges;
	MakeNewtonArrays();
}

CEqofst::~CEqofst(){
	DestroyNewtonArrays();
}

void CEqofst::MakeNewtonArrays(){
	int iQ;
	gslmatrix=new CGSLMatrix_Real *[MaxNcharges+2];
	x=new double[MaxNcharges+1];
	y=new double[MaxNcharges+1];
	rhotry=new double[MaxNcharges+1];
	M=new double*[MaxNcharges+1];
	for(iQ=0;iQ<=MaxNcharges;iQ++){
		x[iQ]=y[iQ]=0.0;
		M[iQ]=new double[MaxNcharges+1];
		gslmatrix[iQ+1]=NULL;
	}
}

void CEqofst::DestroyNewtonArrays(){
	int iQ;
	if(rhotry!=NULL) delete [] rhotry;
	if(y!=NULL) delete [] y; if(x!=NULL) delete [] x;
	if(M!=NULL){
		for(iQ=0;iQ<=MaxNcharges;iQ++) delete [] M[iQ];
		delete [] M;
	}
	rhotry=NULL; x=NULL; y=NULL; M=NULL;
	delete(gslmatrix);
}

void CEqofst::FindXwhereYeqMX(int NDIM){
	int iQ,jQ;


	if(gslmatrix[NDIM]==NULL){
		gslmatrix[NDIM]=new CGSLMatrix_Real(NDIM);
	}
	gslmatrix[NDIM]->SolveLinearEqs(y,M,x);
}

void CEqofst::Calc_of_TMu(CIntrinsic *intr){
	if(intr==NULL || intr->species!=species){
		if(intr!=NULL) delete(intr);
		intr=new CIntrinsic(species);
	}
	FreeGasCalc_of_TMu(intr);
}

void CEqofst::FreeGasCalc_of_TMu(CIntrinsic *intr){
	double pi,ei,dedTi,fi,mui,sigma2i,sigma2=0.0;
	species=intr->species;
	int i,iQ,Ncharges=species->Ncharges,Nspecies=species->Nspecies;
	intr->epsilon=intr->P=intr->dedT=intr->sdens=0.0;
	for(iQ=0;iQ<Ncharges;iQ++){
		intr->rho[iQ]=0.0;
	}
	for(i=0;i<Nspecies;i++){
		mui=0.0;
		for(iQ=0;iQ<Ncharges;iQ++) mui+=species->Q[iQ][i]*intr->mu[iQ];
		fi=exp(mui/intr->T)*species->degen[i];
		freegascalc_onespecies(species->m[i],intr->T,pi,ei,intr->density[i],sigma2i,dedTi);
		intr->epsilon+=ei*fi;
		intr->P+=pi*fi;
		intr->density[i]*=fi;
		intr->flux[i]=fi*2.0*PI*pow(intr->T,3)*exp(-species->m[i]/intr->T);
		intr->dedT+=fi*dedTi;
		intr->sdens+=(fi*pi+fi*ei-mui*intr->density[i])/intr->T;
		sigma2+=fi*sigma2i;
		for(iQ=0;iQ<Ncharges;iQ++){
			intr->rho[iQ]+=intr->density[i]*species->Q[iQ][i];
		}
	}
	intr->asigma=sqrt(sigma2);
}

void CEqofst::Calc_of_EMu(CIntrinsic *intr){
	if(intr==NULL || intr->species!=species){
		if(intr!=NULL) delete(intr);
		intr=new CIntrinsic(species);
	}
	FreeGasCalc_of_EMu(intr);
}

void CEqofst::FreeGasCalc_of_EMu(CIntrinsic *intr){
	double pi,ei,dedTi,fi,mui,e,dT,sigma2i;
	int itry=0,itrymax=10;
	species=intr->species;
	int i,iQ,Ncharges=species->Ncharges,Nspecies=species->Nspecies;
	do{
		e=intr->dedT=0.0;

		for(i=0;i<Nspecies;i++){
			if(itry>itrymax){
				printf("did not converge in %d iterations\n",itrymax);
				intr->Print();
				exit(1);
			}
			mui=0.0;
			for(iQ=0;iQ<Ncharges;iQ++) mui+=species->Q[iQ][i]*intr->mu[iQ];
			fi=exp(mui/intr->T)*species->degen[i];
			freegascalc_onespecies(species->m[i],intr->T,pi,ei,intr->density[i],sigma2i,dedTi);
			e+=ei*fi;
			intr->dedT+=dedTi*fi;
		}
		dT=(intr->epsilon-e)/intr->dedT;
		while(fabs(dT)>0.5*intr->T) dT=0.5*dT;
		intr->T+=dT;
		itry+=1;
	} while(fabs(intr->epsilon-e)/intr->epsilon>1.0E-6);
	FreeGasCalc_of_TMu(intr);
}

void CEqofst::Calc_of_SAlpha(CIntrinsic *intr){
	if(intr==NULL || intr->species!=species){
		if(intr!=NULL) delete(intr);
		intr=new CIntrinsic(species);
	}
	FreeGasCalc_of_SAlpha(intr);
}

void CEqofst::FreeGasCalc_of_SAlpha(CIntrinsic *intr){
	double pi,ei,dedTi,alphai,e,dT,s,dsdT,T0,mui,fi,sigma2i;
	double *f=NULL,*alpha=NULL;
	species=intr->species;
	int i,iQ,Ncharges=species->Ncharges,Nspecies=species->Nspecies;
	int itry=0,itrymax=10;
	T0=intr->T;
	f=new double[Nspecies];
	alpha=new double[Nspecies];
	for(i=0;i<Nspecies;i++){
		alpha[i]=0.0;
		for(iQ=0;iQ<Ncharges;iQ++)
			alpha[i]-=species->Q[iQ][i]*intr->mu[iQ]/intr->T;
		f[i]=species->degen[i]*exp(-alpha[i]);
	}
	do{
		if(itry>itrymax){
			printf("did not converge in %d iterations\n",itrymax);
			intr->Print();
			exit(1);
		}
		s=dsdT=0.0;
		for(i=0;i<Nspecies;i++){
			freegascalc_onespecies(species->m[i],intr->T,pi,ei,intr->density[i],sigma2i,dedTi);
			intr->density[i]*=f[i];
			dedTi*=f[i];
			pi*=f[i];
			ei*=f[i];
			s+=(ei+pi)/intr->T;
			s+=alpha[i]*intr->density[i];
			dsdT+=(dedTi+(alpha[i]/intr->T)*ei)/intr->T;
		}
		dT=(intr->sdens-s)/dsdT;
		while(fabs(dT)>0.5*intr->T) dT=0.5*dT;
		intr->T+=dT;
	} while(fabs(intr->sdens-s)/intr->sdens>1.0E-6);

	for(iQ=0;iQ<Ncharges;iQ++) intr->mu[iQ]=alpha[iQ]*intr->T;

	FreeGasCalc_of_TMu(intr);
	if(alpha!=NULL) delete [] alpha;
	if(f!=NULL) delete [] f;
}

void CEqofst::Calc_of_SRho(CIntrinsic *intr){
	if(intr==NULL || intr->species!=species){
		if(intr!=NULL) delete(intr);
		intr=new CIntrinsic(species);
	}
	FreeGasCalc_of_SRho(intr);
}

void CEqofst::FreeGasCalc_of_SRho(CIntrinsic *intr){
	double stry,error,delerror,delm,si,pi,ei,xarg,X,dedti,sigma2i;
	int itry,ntries,itrymax=20;
	species=intr->species;
	int i,iQ,jQ,Ncharges=species->Ncharges,Nspecies=species->Nspecies;
	ntries=0;
	do{
		stry=intr->epsilon=intr->P=intr->dedT=0.0;
		for(iQ=0;iQ<=Ncharges;iQ++)
			for(jQ=0;jQ<=Ncharges;jQ++) M[iQ][jQ]=0.0;
		for(iQ=0;iQ<Ncharges;iQ++) rhotry[iQ]=0.0;

		for(i=0;i<Nspecies;i++){
			X=species->degen[i];
			xarg=0.0;
			for(iQ=0;iQ<Ncharges;iQ++) xarg+=species->Q[iQ][i]*intr->mu[iQ]/intr->T;
			X*=exp(xarg);
			freegascalc_onespecies(species->m[i],intr->T,pi,ei,intr->density[i],sigma2i,dedti);
			pi*=X;
			ei*=X;
			dedti*=X;
			intr->density[i]*=X;
			si=(pi+ei)/intr->T;
			for(iQ=0;iQ<Ncharges;iQ++){
				si-=intr->density[i]*species->Q[iQ][i]*intr->mu[iQ]/intr->T;
				rhotry[iQ]+=intr->density[i]*species->Q[iQ][i];
			}
			stry+=si;
			intr->dedT+=dedti;
			intr->P+=pi;
			intr->epsilon+=ei;
			for(iQ=0;iQ<Ncharges;iQ++){
				M[Ncharges][iQ]+=species->Q[iQ][i]*(si-pi/intr->T)/intr->T;
				M[Ncharges][Ncharges]+=intr->mu[iQ]*species->Q[iQ][i]*intr->density[i]/(intr->T*intr->T);
				M[Ncharges][Ncharges]-=(intr->mu[iQ]/(intr->T*intr->T))*species->Q[iQ][i]*(si+ei/intr->T);
				for(jQ=0;jQ<=iQ;jQ++)
					M[iQ][jQ]+=intr->density[i]*species->Q[iQ][i]*species->Q[jQ][i]/intr->T;
			}
			M[Ncharges][Ncharges]+=dedti/intr->T;

		}

		for(iQ=0;iQ<Ncharges;iQ++)
			for(jQ=iQ+1;jQ<=Ncharges;jQ++) M[iQ][jQ]=M[jQ][iQ];

		y[Ncharges]=stry-intr->sdens;
		for(iQ=0;iQ<Ncharges;iQ++) y[iQ]=rhotry[iQ]-intr->rho[iQ];

		FindXwhereYeqMX(Ncharges+1);

		intr->T-=x[Ncharges];
		for(iQ=0;iQ<Ncharges;iQ++) intr->mu[iQ]-=x[iQ];

		error=0.0;
		for(iQ=0;iQ<=Ncharges;iQ++){
			delerror=y[iQ]*intr->T/intr->epsilon;
			error+=delerror*delerror;
		}

		ntries+=1;
	} while(fabs(error)>1.0E-8 && ntries<itrymax);
	if(ntries==itrymax){
		printf("did not converge in %d iterations, error=%g\n",
			ntries,error);
		intr->Print();
		exit(1);
	}

	FreeGasCalc_of_TMu(intr);
}

void CEqofst::Calc_of_ERho(CIntrinsic *intr){
	if(intr==NULL || intr->species!=species){
		if(intr!=NULL) delete(intr);
		intr=new CIntrinsic(species);
	}
	FreeGasCalc_of_ERho(intr);
}

void CEqofst::FreeGasCalc_of_ERho(CIntrinsic *intr){
	double etry,error,pi,ei,dedti,xarg,X,delerror,oldT,sigma2i;
	CSpecies *species=intr->species;
	int i,iQ,jQ,Ncharges=species->Ncharges,Nspecies=species->Nspecies;
	int itry,ntries,itrymax=40;
	ntries=0;
	double emin=0.0;
	for(i=0;i<Nspecies;i++) emin+=intr->density[i]*species->m[i];
	//if(emin>intr->epsilon){
		//printf("CEqofst::FreeGasCalc_of_ERho,emin=%g, epsilon=%g\n",emin,intr->epsilon);
		//exit(1);
	//}
	do{
		etry=intr->dedT=0.0;
		for(iQ=0;iQ<Ncharges;iQ++) rhotry[iQ]=0.0;
		for(iQ=0;iQ<=Ncharges;iQ++){
			for(jQ=0;jQ<=Ncharges;jQ++) M[iQ][jQ]=0.0;
		}

		for(i=0;i<Nspecies;i++){
			X=species->degen[i];
			xarg=0.0;
			for(iQ=0;iQ<Ncharges;iQ++) xarg+=species->Q[iQ][i]*intr->mu[iQ]/intr->T;
			X*=exp(xarg);
			freegascalc_onespecies(species->m[i],intr->T,pi,ei,intr->density[i],sigma2i,dedti);
			intr->density[i]*=X;
			for(iQ=0;iQ<Ncharges;iQ++){
				rhotry[iQ]+=intr->density[i]*species->Q[iQ][i];
				M[Ncharges][iQ]-=ei*X*species->Q[iQ][i];
				for(jQ=0;jQ<=iQ;jQ++)
					M[iQ][jQ]-=species->Q[iQ][i]*species->Q[jQ][i]*intr->density[i];
			}
			intr->dedT+=dedti*X;
			etry+=ei*X;
		}

		for(iQ=0;iQ<Ncharges;iQ++)
			for(jQ=iQ+1;jQ<=Ncharges;jQ++) M[iQ][jQ]=M[jQ][iQ];
		M[Ncharges][Ncharges]=-intr->T*intr->T*intr->dedT;

		for(iQ=0;iQ<Ncharges;iQ++) y[iQ]=rhotry[iQ]-intr->rho[iQ];
		y[Ncharges]=etry-intr->epsilon;

		FindXwhereYeqMX(Ncharges+1);		

		for(iQ=0;iQ<Ncharges;iQ++) intr->mu[iQ]=(-intr->mu[iQ]/intr->T)-x[iQ];
		oldT=intr->T;
		intr->T=1.0/((1.0/intr->T)-x[Ncharges]);
		if(intr->T<1.0){
			//printf("looking for very small T=%g, epsilon=%g\n",intr->T,intr->epsilon);
			intr->T=0.5*oldT;
				//exit(1);
		}
		for(iQ=0;iQ<Ncharges;iQ++) intr->mu[iQ]=-intr->mu[iQ]*intr->T;

		error=0.0;
		for(iQ=0;iQ<=Ncharges;iQ++){
			if(iQ<Ncharges) delerror=y[iQ]*intr->T/intr->epsilon;
			else delerror=y[iQ]/intr->epsilon;
			error+=delerror*delerror;
		}
		//printf("error=%g\n",error);

		ntries+=1;
	} while(fabs(error)>1.0E-8 && ntries<itrymax);
	if(ntries==itrymax && fabs(error)>1.0E-6){
		printf("In FreeGasCalc_of_ERho, did not converge in %d iterations, error=%g\n",ntries,error);
		intr->Print();
		exit(1);
	}
	FreeGasCalc_of_TMu(intr);
}

void CEqofst::Calc_of_ERhoOverS(CIntrinsic *intr){
	if(intr==NULL || intr->species!=species){
		if(intr!=NULL) delete(intr);
		intr=new CIntrinsic(species);
	}
	FreeGasCalc_of_ERhoOverS(intr);
}

void CEqofst::FreeGasCalc_of_ERhoOverS(CIntrinsic *intr){
	double etry,error,ps,es,ds,dedts,xarg,Xs,delerror,oldT,sigma2s,sdens,T,dedt,ri,stry;
	double e_target;
	double **B,*A,*C;
	CSpecies *species=intr->species;
	int s,iQ,jQ,Ncharges=species->Ncharges,Nspecies=species->Nspecies;
	int itry,ntries,itrymax=20;
	double *rhotry=new double[Ncharges],*alpha=new double[Ncharges];
	A=new double[Ncharges];
	C=new double[Ncharges];
	B=new double *[Ncharges];
	double *rhoovers_target=new double [Ncharges];
	for(iQ=0;iQ<Ncharges;iQ++) B[iQ]=new double[Ncharges];
	for(iQ=0;iQ<Ncharges;iQ++){
		rhoovers_target[iQ]=intr->rho[iQ]/intr->sdens;
	}
	e_target=intr->epsilon;
	for(iQ=0;iQ<Ncharges;iQ++) alpha[iQ]=-intr->mu[iQ]/intr->T;
	T=intr->T;
	ntries=0;
	do{
		/* printf("T=%g",T);
		for(iQ=0;iQ<Ncharges;iQ++) printf(", alpha[%d]=%g",iQ,alpha[iQ]);
		printf("\n"); */
		for(iQ=0;iQ<Ncharges;iQ++){
			A[iQ]=C[iQ]=rhotry[iQ]=0.0;
			for(jQ=0;jQ<Ncharges;jQ++) B[iQ][jQ]=0.0;
		}
		etry=stry=dedt=0.0;
		for(s=0;s<Nspecies;s++){
			Xs=species->degen[s];
			xarg=0.0;
			for(iQ=0;iQ<Ncharges;iQ++) xarg-=species->Q[iQ][s]*alpha[iQ];
			Xs*=exp(xarg);
			freegascalc_onespecies(species->m[s],T,ps,es,ds,sigma2s,dedts);
			ps*=Xs;
			es*=Xs;
			ds*=Xs;
			dedts*=Xs;
			dedt+=dedts;
			etry+=es;
			stry+=(ps+es)/T;
			for(iQ=0;iQ<Ncharges;iQ++){
				rhotry[iQ]+=ds*species->Q[iQ][s];
				C[iQ]-=es*species->Q[iQ][s];
				for(jQ=0;jQ<Ncharges;jQ++){
					A[iQ]-=ds*species->Q[iQ][s]*species->Q[jQ][s]*alpha[jQ];
					B[iQ][jQ]-=ds*species->Q[iQ][s]*species->Q[jQ][s];
				}
			}
		}
		for(iQ=0;iQ<Ncharges;iQ++) stry+=rhotry[iQ]*alpha[iQ];

		for(iQ=0;iQ<Ncharges;iQ++){
			M[Ncharges][iQ]=C[iQ];
			ri=rhotry[iQ]/(stry*stry);
			M[iQ][Ncharges]=(C[iQ]/stry)+ri*T*dedt;
			for(jQ=0;jQ<Ncharges;jQ++){
				M[iQ][Ncharges]-=ri*alpha[jQ]*C[jQ];
				M[iQ][jQ]=(B[iQ][jQ]/stry)-ri*((C[jQ]/T)+A[jQ]);
			}
		}
		M[Ncharges][Ncharges]=-T*T*dedt;

		for(iQ=0;iQ<Ncharges;iQ++) y[iQ]=(rhotry[iQ]/stry)-rhoovers_target[iQ];
		y[Ncharges]=etry-e_target;

		FindXwhereYeqMX(Ncharges+1);		

		oldT=T;
		T=1.0/((1.0/T)-x[Ncharges]);
		for(iQ=0;iQ<Ncharges;iQ++){
			alpha[iQ]-=x[iQ];
		}
		/*
		printf("epsilon=%g",etry);
		for(iQ=0;iQ<Ncharges;iQ++) printf(", rhoovers[%d]=%g",iQ,rhotry[iQ]/stry);
		printf("\n");
		printf("TARGETS: epsilon=%g",e_target);
		for(iQ=0;iQ<Ncharges;iQ++) printf(", rhoovers[%d]=%g",iQ,rhoovers_target[iQ]);
		printf("\n");
		printf("--------------------------------------------------------------------------------\n");*/
		if(T<1.0){
			//printf("looking for very small T=%g, epsilon=%g\n",intr->T,intr->epsilon);
			T=0.5*oldT;
				//exit(1);
		}

		error=0.0;
		for(iQ=0;iQ<=Ncharges;iQ++){
			if(iQ<Ncharges) delerror=y[iQ];
			else delerror=y[iQ]/e_target;
			error+=delerror*delerror;
		}
		//printf("error=%g\n",error);

		ntries+=1;
	} while(fabs(error)>1.0E-12 && ntries<itrymax);
	if(ntries==itrymax && fabs(error)>1.0E-6){
		printf("In FreeGasCalc_of_ERhoOverS, did not converge in %d iterations, error=%g\n",ntries,error);
		intr->Print();
		exit(1);
	}

	intr->T=T;
	for(iQ=0;iQ<Ncharges;iQ++) intr->mu[iQ]=-intr->T*alpha[iQ];
	FreeGasCalc_of_TMu(intr);
	//intr->Print();

	delete [] rhotry; delete [] alpha; delete rhoovers_target; delete [] A; delete [] C;
	for(iQ=0;iQ<Ncharges;iQ++) delete [] B[iQ];
	delete [] B;

}

void CEqofst::CalcVisc(CIntrinsic *intr){
	if(intr==NULL || intr->species!=species){
		if(intr!=NULL) delete(intr);
		intr=new CIntrinsic(species);
	}
	intr->B=intr->eta=0.0;  
}

// This calculates intrinsic properties of a 1-species ideal gas of
// classical particles of mass m at temperature T
void CEqofst::freegascalc_onespecies(double m,double t,double &p,double &e,double &dens,double &sigma2,double &dedt){
	const double prefactor=1.0/(2.0*PI*PI*pow(HBARC,3));
	double k0,k1,z,k0prime,k1prime,m2,m3,m4,t2,t3,I1,I2,Iomega;
	m2=m*m;
	m3=m2*m;
	m4=m2*m2;
	t2=t*t;
	t3=t2*t;
	z=m/t;
	if(z>1000.0){
		p=e=dens=dedt=0.0;
		printf("z is huge=%g, m=%g, t=%g\n",z,m,t);
	}
	else{
		if(z<0.0){
			printf("___z=%g,m=%g,T=%g ___\n",z,m,t);
			exit(1);
		}
		k0=Bessel::K0(z);
		k1=Bessel::K1(z);
		p=prefactor*(m2*t2*k0+2.0*m*t3*k1);
		e=prefactor*(3.0*m2*t2*k0+(m3*t+6.0*m*t3)*k1);
		dens=p/t;
		k0prime=-k1;
		k1prime=-k0-k1/z;
		dedt=prefactor*(6.0*m2*t*k0+(m3+18.0*m*t2)*k1-3.0*m3*k0prime-((m4/t)+6.0*m2*t)*k1prime);
		Iomega=exp(-m/t)/(30.0*PI*PI*HBARC*HBARC*HBARC);
		I1=pow(m,1.5)*pow(t,3.5)*7.5*sqrt(2.0*PI);
		I2=24.0*pow(t,5);
		sigma2=Iomega*(I1+I2+0.5*sqrt(I1*I2));  // this is an approximation (+/-2%) to messy integral
	}
}

void CEqofst::PECheck(CIntrinsic *intr){
	double Pcheck=0.0,echeck=0.0;
	double pmag,E,dp=1.0,mui,mass;
	double prefactor,fact2;
	species=intr->species;
	int i,iQ,Ncharges=species->Ncharges,Nspecies=species->Nspecies;
	for(i=0;i<Nspecies;i++){
		mui=0.0;
		mass=species->m[i];
		for(iQ=0;iQ<Ncharges;iQ++) mui+=species->Q[iQ][i]*intr->mu[iQ];
		prefactor=(species->degen[i]/(2.0*PI*PI*HBARC*HBARC*HBARC));
		prefactor*=exp(mui/intr->T);
		for(pmag=0.5*dp;pmag<10000;pmag+=dp){
			E=sqrt(pmag*pmag+mass*mass);
			fact2=pmag*pmag*dp*exp(-E/intr->T);
			Pcheck+=prefactor*fact2*pmag*pmag/(3.0*E);
			echeck+=prefactor*fact2*E;
		}
	}
	printf("Check: epsilon=%g =? %g,   P=%g =? %g\n",intr->epsilon,echeck,intr->P,Pcheck);
}

void CEqofst::InitIntrinsic_EMu0(CIntrinsic *&intrinsic,double epsilon){
	int iQ;
	if(intrinsic!=NULL) delete(intrinsic);
	intrinsic=new CIntrinsic(intrinsic->species);
	intrinsic->epsilon=epsilon;
	for(iQ=0;iQ<intrinsic->species->Ncharges;iQ++) intrinsic->mu[iQ]=0.0;
	Calc_of_EMu(intrinsic);
}

//
// \lambda_{ij}=\delta T_{ij}/(\epsilon+P) \frac{1}{\lambda_{\rm fact}} 

double CEqofst::GetLambdaFact(CIntrinsic *intr){
	CSpecies *species=intr->species;
	int ispecies,iQ,n,i;
	const int nmax=70;
	double G[nmax+5];
	double lambdafact,z,Ipp=0.0,dIpp,J,nfact,sign,alpha;
	for(ispecies=0;ispecies<species->Nspecies;ispecies++){
		z=species->m[ispecies]/intr->T;
		alpha=0.0;
		for(iQ=0;iQ<species->Ncharges;iQ++){
			alpha+=intr->mu[iQ]*species->Q[iQ][ispecies]/intr->T;
		}

		G[0]=gsl_sf_gamma_inc(5,z)*pow(z,-5);
		for(int i=1;i<nmax+5;i++){
			n=5.0-2.0*i;
			if(n!=-1)	G[i]=(-exp(-z)/n)+(G[i-1]*z*z-z*exp(-z))/((n+1.0)*n);
			else G[i]=gsl_sf_gamma_inc(-1,z)*z;
		}

		J=0.0;
		nfact=1.0;
		sign=1.0;
		for(n=0;n<nmax;n+=1.0){
			if(n>0) sign=-1.0;
			J+=sign*nfact*(G[n]-2.0*G[n+1]+G[n+2]);
			nfact=nfact*0.5/(n+1.0);
			if(n>0) nfact*=(2.0*n-1.0);
		}
		Ipp=Ipp+species->degen[ispecies]*exp(alpha)
			*pow(species->m[ispecies],4)*(-z*J+15.0*gsl_sf_bessel_Kn(2,z)/(z*z));	
	}
	lambdafact=Ipp/(60.0*PI*PI*HBARC*HBARC*(intr->epsilon+intr->P));
	return lambdafact;	
}

// \int d^3p p^4/e^3 exp(-e/T)
double CEqofst::GetBfact(double m,double T){
	const int nmax=200;
	double nfact=1.0,Itest,a,I1test,I2test,I3test,atest;
	double x,dx=0.004,z=m/T;
	double XX=exp(-z),XXX,alpha;
	double Gminus[2*nmax+1];
	double prefactor=1.0/(2.0*PI*PI*HBARC*HBARC*HBARC);
	double k0=Bessel::K0(z);
	double k1=Bessel::K1(z);
	double k2=Bessel::Kn(2,z);
	double pressure=prefactor*(m*m*T*T*k0+2.0*m*T*T*T*k1);
	int n;
	//Gminus[m]=Gamma(-m,z)
	Gminus[0]=gsl_sf_expint_E1(z);
	Itest=0.0;
	for(x=1.0+0.5*dx;x<100/z;x+=dx){
		Itest+=dx*exp(-z*x)/x;
	}
	printf("Gamma(-%d,%g)=%g =? %g\n",0,z,Gminus[0],Itest);
	for(n=1;n<2*nmax+1;n++){
		Gminus[n]=(-1.0/double(n))*(Gminus[n-1]-XX/pow(z,n));
		Itest=0.0;
		for(x=1.0+0.5*dx;x<100;x+=dx){
			Itest+=dx*exp(-z*x)/pow(x,n+1);
		}
		printf("Gamma(-%d,%g)=%g =? %g\n",n,z,Gminus[n],Itest);
	}

	double I1,I2,I3;
	I1=3.0*m*m*Bessel::Kn(1,z)/z;
	I3=-m*m*m/T*Bessel::Kn(2,z)/z;
	I2=m*m*XX;
	alpha=-0.5*m*m;
	for(n=1;n<nmax;n++){
		I2+=alpha*pow(z,2*n)*Gminus[2*n-1];
		printf("n=%d, alpha=%g, Gminus=%g\n",n,alpha,Gminus[2*n-1]);
		alpha=alpha*double(n-0.5)/double(n+1.0);
	}
	a=I1+I2+I3;

	double dp=0.004*T,p,e;
	I1test=I2test=I3test=Itest=atest=0.0;
	for(p=0.5*dp;p<40.0*T+20.0*sqrt(m*T);p+=dp){
		e=sqrt(m*m+p*p);
		XXX=exp(-e/T);
		I1test+=dp*p*p*XXX/e;
		I2test+=dp*p*p*XXX/(e*e);
		I3test+=dp*p*p*XXX;
		Itest+=dp*XXX*pow(p,4)/pow(e,3);
		atest+=dp*p*p*((2.0*p*p/e)-(p*p*p*p/(e*e*e)))*XXX;
	}
	I1test=I1test*3.0;
	I2test=m*m*I2test/T;
	I3test=-I3test/T;
	printf("I1=%g =? %g, I2=%g =? %g, I3=%g =? %g, a=%g =? %g =? %g\n",I1,I1test,I2,I2test,I3,I3test,a,Itest,I1test+I2test+I3test);

	a=(m*m*prefactor*a/3.0)+pressure;
	atest=atest*prefactor/3.0;
	printf("---------  a=%g, atest=%g, pressure=%g\n",a,atest,pressure);
	return a;
}

double CEqofst::GetSigmaB(CIntrinsic *intrinsic){
	double pi,*e,*rho,*e2,p2i,mi,mui,K0,K1,K2,K3,m3,m4,m5,alpha,dfact,prefactor=1.0/(2.0*PI*PI*HBARC*HBARC*HBARC);
	double etest=0.0,e2test=0.0,ptest=0.0,p2test=0.0,dp2i,rhotest=0.0,dp,p,pmax,z,fact1,alt_zetaovertau=0.0,stot=0.0;
	double p2=0.0,e2tot=0.0,denstot=0.0,P=0.0,etot=0.0,energy;
	const int nmax=100;
	double T=intrinsic->T;
	CSpecies *species=intrinsic->species;
	int ispecies,jspecies,n,Nspecies=species->Nspecies,iQ,Ncharges=species->Ncharges;
	double Gminus[2*nmax+1],XX,XXX,EE,*E,*dPdrho,dPde;
	E=new double[Nspecies];
	e=new double[Nspecies];
	rho=new double[Nspecies];
	e2=new double[Nspecies];
	EE=0.0;
	dPdrho=new double[Nspecies];
	for(jspecies=0;jspecies<Nspecies;jspecies++) E[jspecies]=dPdrho[jspecies]=0.0;
	for(ispecies=0;ispecies<Nspecies;ispecies++){
		mi=species->m[ispecies];
		dfact=prefactor*species->degen[ispecies];
		mui=0.0;
		for(iQ=0;iQ<Ncharges;iQ++) mui+=intrinsic->mu[iQ]*species->Q[iQ][ispecies];
		dfact*=exp(mui/T);
		//printf("mi=%g, degen=%g, dfact=%g\n",mi,species->degen[ispecies],dfact);
		m3=pow(mi,3);
		m4=m3*mi;;
		m5=m4*mi;
		z=mi/T;
		K0=Bessel::K0(z);
		K1=Bessel::K1(z);
		K2=K0+2.0*K1/z;
		K3=K1+4.0*K2/z;
		e[ispecies]=m4*dfact*((K1/z)+3.0*(K2/(z*z)));
		pi=m4*dfact*((K0/(z*z))+2.0*K1/(z*z*z));
		rho[ispecies]=pi/T;
		e2[ispecies]=m5*dfact*(z*K2+3.0*K3)/(z*z);
		p2i=e2[ispecies]-2.0*mi*mi*rho[ispecies];
		XX=exp(-z);
		fact1=m5*dfact/z;
		p2i+=fact1*XX;
		Gminus[0]=gsl_sf_expint_E1(z);
		for(n=1;n<2*nmax+1;n++){
			Gminus[n]=(-1.0/double(n))*(Gminus[n-1]-XX/pow(z,n));
			//printf("Gminus[%d]=%g\n",n,Gminus[n]);
		}
		alpha=-0.5*fact1;
		for(n=1;n<nmax;n++){
			dp2i=alpha*pow(z,2*n)*Gminus[2*n-1];
			p2i+=dp2i;
			//printf("alpha[%d]*Gminus[]*z^2n=%g\n",n,alpha*Gminus[2*n-1]*pow(z,2*n));
			alpha*=double(n-0.5)/double(n+1);
		}
		p2i+=dp2i*nmax*2.0/3.0;
		p2i=p2i/9.0;
		etest=0.0;e2test=0.0;ptest=0.0;p2test=0.0;rhotest=0.0;
		dp=0.001*(sqrt(mi*T)+T);
		pmax=50.0*(sqrt(mi*T)+T);
		for(p=0.5*dp;p<pmax;p+=dp){
			energy=sqrt(p*p+mi*mi);
			XXX=exp(-energy/T);
			rhotest+=dfact*p*p*dp*XXX;
			ptest+=dfact*p*p*p*p*dp*XXX/(3.0*energy);
			etest+=dfact*p*p*energy*dp*XXX;
			e2test+=dfact*p*p*energy*energy*dp*XXX;
			p2test+=dfact*pow(p,6)*dp*XXX/(9.0*energy*energy);
		}
		p2i=p2test;
		P+=pi; etot+=e[ispecies]; denstot+=rho[ispecies]; e2tot+=e2[ispecies]; p2+=p2i;
		stot+=(pi+e[ispecies]-mui*rho[ispecies])/T;
		//printf("m=%g | rho=%g=?%g | P=%g=?%g | e=%g=?%g | e2=%g=?%g, ratio=%g | p2=%g=?%g, ratio=%g\n", mi,rho[ispecies],rhotest,pi,ptest,e[ispecies],etest,e2[ispecies],e2test,e2[ispecies]/e2test,p2i,p2test,p2i/p2test);
	}
	alt_zetaovertau=p2;
	EE=0.0;
	for(ispecies=0;ispecies<Nspecies;ispecies++){
		for(jspecies=0;jspecies<Nspecies;jspecies++){
			if(jspecies!=ispecies) E[ispecies]+=e[jspecies]*e[jspecies]/rho[jspecies];
		}
		EE+=e[ispecies]*e[ispecies]/rho[ispecies];
		dPdrho[ispecies]=T+P*T*e[ispecies]/(rho[ispecies]*(E[ispecies]-e2tot)+e[ispecies]*e[ispecies]);
		//printf("E[%d]=%g\n",ispecies,E[ispecies]);
	}
	//printf("epsilon=%g, P=%g,EE=%g\n",etot,P,EE);
	dPde=P*T/(e2tot-EE);
	double sigma2=0.0;
	sigma2+=p2+dPde*dPde*e2tot;
	denstot=0.0;
	for(ispecies=0;ispecies<Nspecies;ispecies++){
		mi=species->m[ispecies];
		sigma2+=dPdrho[ispecies]*dPdrho[ispecies]*rho[ispecies];
		sigma2-=2.0*T*rho[ispecies]*dPdrho[ispecies];
		sigma2+=2.0*e[ispecies]*dPdrho[ispecies]*dPde;
		sigma2-=2.0*dPde*(e2[ispecies]-mi*mi*rho[ispecies])/3.0;
		intrinsic->density[ispecies]=rho[ispecies];
		alt_zetaovertau-=T*dPdrho[ispecies]*rho[ispecies];
		denstot+=intrinsic->density[ispecies];
	}
	printf("s=%g, denstot=%g, s/denstot=%g\n",intrinsic->sdens,denstot,intrinsic->sdens/denstot);
	//printf("check P=%g=?%g, dPdrho[0]=%g, dPde=%g\n",P,denstot*T,dPdrho[0],dPde);
	delete [] E;
	delete [] dPdrho;
	delete [] rho;
	delete [] e;
	delete [] e2;
	intrinsic->P=P;
	intrinsic->bsigma=sqrt(sigma2);
	intrinsic->epsilon=etot;
	alt_zetaovertau-=T*dPde*(P+etot);
	double altalt_zetaovertau=(p2/T)-P-P*dPde;
	double zeta=(sigma2/T)/(2.5*HBARC*denstot);
	//printf("p2/PT=%g, dPde=%g\n",p2/(P*T),dPde);
	//printf("T=%g, sigma2=%g, sigma2/PT=%g =? %g =? %g, 4pi*zeta/s=%g\n",T,sigma2,sigma2/(P*T),alt_zetaovertau/(P*T),altalt_zetaovertau/P,4.0*PI*zeta/stot);
	return sigma2;
}

#endif
