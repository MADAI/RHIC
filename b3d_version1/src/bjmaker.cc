#ifndef __bjmaker_cc__
#define __bjmaker_cc__
#define __JOSH_FORMAT__

#include "b3d.h"
using namespace std;

void CBjMaker::Init(){
	double pi,ei,sigma2i,dedti,degen;
	int ires;
	randy=b3d->randy;
	nwrite=0;
	nsample=b3d->NSAMPLE;
	reslist=b3d->reslist;
	nres=reslist->NResonances;
	etamax=parameter::getD(b3d->parmap,"B3D_ETAMAX",1.0);
	etagauss=parameter::getD(b3d->parmap,"B3D_BJMAKER_ETAGAUSS",2.0);
	gaussian=parameter::getB(b3d->parmap,"B3D_BJMAKER_GAUSSIAN",false);
	balance=parameter::getB(b3d->parmap,"B3D_BJMAKER_BALANCE",true);
	tau=parameter::getD(b3d->parmap,"B3D_BJMAKER_TAU",1.0);
	Rx=parameter::getD(b3d->parmap,"B3D_BJMAKER_Rx",1.5);
	Ry=parameter::getD(b3d->parmap,"B3D_BJMAKER_Ry",3.0);
	T=parameter::getD(b3d->parmap,"B3D_BJMAKER_T",160.0);
	//printf("epsilon_H=%g\n",intrinsic->epsilon);
	density=new double[nres];
	//ID=new int[nres];
	epsilon=P=denstot=0.0;
	CResInfo *resinfo=reslist->GfirstResInfoptr;
	for(ires=0;ires<nres;ires++){
		density[ires]=0.0;
		//ID[ires]=resinfo->code;
		if(resinfo->code!=22){
			degen=2.0*resinfo->spin+1.0;
			//printf("ires=%d, ID=%d, degen=%g, mass=%g, T=%g\n",ires,resinfo->code,degen,resinfo->mass,T);
			b3d->freegascalc_onespecies(resinfo->mass,T,pi,ei,density[ires],sigma2i,dedti);
			//printf("mass=%g, T=%g, pi=%g, ei=%g, di=%g\n",resinfo->mass,T,pi,ei,density[ires]);
			density[ires]*=degen;
			denstot+=density[ires];
			P+=degen*pi;
			epsilon+=degen*ei;
		}
		resinfo=resinfo->nextResInfoptr;
	}
}

int CBjMaker::MakeEvent(){
	int nparts;
	b3d->Reset();
	nwrite=0;
	if(gaussian){
		nparts=b3d->NPARTSMAX;
		if(balance){
			GenerateParticles_Gaussian_Balance(nparts);	
		}
		else GenerateParticles_Gaussian(nparts);
	}
	else{
		nparts=GenerateParticles();
	}
	return nparts;
}

// \lambda_{ij}=\delta T_{ij}/(\epsilon+P) \frac{1}{\lambda_{\rm fact}}

int CBjMaker::GenerateParticles(){
	double volume=PI*Rx*Ry*tau*2.0*etamax;
	double x,y,m,mt,eta;
	double rapidity,p[4];
	double remainder;
	int Ni,i,nparts=0,ires;
	CResInfo *resinfo=reslist->GfirstResInfoptr;
	for(ires=0;ires<nres;ires++){
		m=resinfo->mass;
		remainder=(density[ires]*volume*double(nsample))-floor(density[ires]*volume*double(nsample));
		Ni=lrint(floor(density[ires]*volume*nsample));
		if(randy->ran()<remainder) Ni+=1;
		for(i=0;i<Ni;i++){
			randy->gauss2(&x,&y);
				//do{
				//x=(1.0-2.0*randy->ran());
				//y=(1.0-2.0*randy->ran());
				//}while(x*x+y*y>1.0);
			x*=Rx;
			y*=Ry;
			if(b3d->BJORKEN) eta=etamax*(1.0-2.0*randy->ran());
			else eta=etagauss*randy->gauss();
			randy->generate_boltzmann(m,T,p);
			mt=sqrt(m*m+p[1]*p[1]+p[2]*p[2]);
			rapidity=eta+asinh(p[3]/mt);
			b3d->partarray[nparts]->Init(resinfo->code,x,y,tau,eta,p[1],p[2],resinfo->mass,rapidity);
			nparts+=1;
		}
		resinfo=resinfo->nextResInfoptr;
	}
	return nparts;
}

void CBjMaker::GenerateParticles_Gaussian(int nparts){
	double volume=PI*Rx*Ry*tau*2.0*etamax;
	double x,y,m,mt,eta;
	double rapidity,p[4];
	double ranguy,denssum;
	int Ni,i,ipart,ires;
	CResInfo *resinfo;
	for(ipart=0;ipart<nparts;ipart++){
		ranguy=randy->ran()*denstot;
		ires=0;
		denssum=0.0;
		denssum=density[0];
		ires=0;
		resinfo=reslist->GfirstResInfoptr;
		while(ranguy>denssum){
			resinfo=resinfo->nextResInfoptr;
			ires+=1;
			denssum+=density[ires];
		}
		m=resinfo->mass;
		randy->gauss2(&x,&y);
			//do{
			//x=(1.0-2.0*randy->ran());
			//y=(1.0-2.0*randy->ran());
			//}while(x*x+y*y>1.0);
		x*=Rx;
		y*=Ry;
		if(b3d->BJORKEN) eta=etamax*(1.0-2.0*randy->ran());
		else{
			do{
				eta=etagauss*randy->gauss();
			} while(fabs(eta)>5.5);
		}
		randy->generate_boltzmann(m,T,p);
		mt=sqrt(m*m+p[1]*p[1]+p[2]*p[2]);
		rapidity=eta+asinh(p[3]/mt);
		b3d->partarray[ipart]->Init(resinfo->code,x,y,tau,eta,p[1],p[2],m,rapidity);
	}
}

void CBjMaker::GenerateParticles_Gaussian_Balance(int nparts){
	double volume=PI*Rx*Ry*tau*2.0*etamax;
	double x,y,m,mt,eta;
	double rapidity,p[4];
	double ranguy,denssum;
	int Ni,i,ipart,ires,ii,ID;
	CResInfo *resinfo;
	for(ipart=0;ipart<nparts/2;ipart++){
		ranguy=randy->ran()*denstot;
		ires=0;
		denssum=density[0];
		ires=0;
		resinfo=reslist->GfirstResInfoptr;
		while(ranguy>denssum){
			resinfo=resinfo->nextResInfoptr;
			denssum+=density[ires];
			ires+=1;
		}
		m=resinfo->mass;
		do{
			x=(1.0-2.0*randy->ran());
			y=(1.0-2.0*randy->ran());
		}while(x*x+y*y>1.0);
		x*=Rx;
		y*=Ry;
		if(b3d->BJORKEN) eta=etamax*(1.0-2.0*randy->ran());
		else{
			do{
				eta=etagauss*randy->gauss();
			} while(fabs(eta)>5.5);
		}
		ID=resinfo->code;
		for(ii=0;ii<2;ii++){
			if(ii==1 && abs(resinfo->code)==211) ID=-ID;
			randy->generate_boltzmann(m,T,p);
			mt=sqrt(m*m+p[1]*p[1]+p[2]*p[2]);
			rapidity=eta+asinh(p[3]/mt);
			b3d->partarray[2*ipart+ii]->Init(ID,x,y,tau,eta,p[1],p[2],m,rapidity);
		}
	}
	printf("CBjMaker::GenerateParticles_Gaussian_Balance, %d particles generated\n",nparts);
}
//if(ii=1 && abs(resinfo->code)==211) ID=-ID;

#endif
