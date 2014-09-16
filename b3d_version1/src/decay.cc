#ifndef __DECAY_CC__
#define __DECAY_CC__
#include "b3d.h"

using namespace std;

void CB3D::Decay(CPart *&mother,int &nbodies, CPart **&daughter){
	const double HBARC=197.326;
	int ibody,alpha;
	double u[4],mass[6],mtot;
	CPart *dptr;

	double *p[6],kprime[4],qprime[4],ptot[4];
	double q,weight,wmax,sthet,cthet,phi;
	double p3mag,kprimemax,p3max,ppmax,kprimemax2,kprimemag2,qprimemax,qprimemax2,qprimemag2,ppmag;
	double e1prime,e2prime,e3prime,e4prime,e1max,e2max,e3max,e4max;
	double e12,u12[4],pp[4];

	mass[0]=mother->GetMass();
	p[0]=mother->p;
	
	/* Create daughter objects */
	mtot=0.0;
	for(ibody=0;ibody<nbodies;ibody++){
		mass[ibody+1]=daughter[ibody]->resinfo->mass;
		mtot+=mass[ibody+1];
		p[ibody+1]=daughter[ibody]->p;
	}
	if(mtot>mass[0]){
		printf("CB3D::Decay, This mass can't decay, mothermass=%g, mtot of products=%g\n",mass[0],mtot);
		exit(1);
	}

	/* TWO-BODY DECAYS */
	if(nbodies==2){
		cthet=1.0-2.0*randy->ran();
		sthet=sqrt(1.0-cthet*cthet);
		phi=2.0*PI*randy->ran();
		q=sqrt(Misc::triangle(mass[0],mass[1],mass[2]));
		p[1][3]=q*cthet;
		p[1][1]=q*sthet*cos(phi);
		p[1][2]=q*sthet*sin(phi);
		p[2][3]=-p[1][3];
		p[2][2]=-p[1][2];
		p[2][1]=-p[1][1];
		p[1][0]=sqrt(mass[1]*mass[1]+p[1][1]*p[1][1]+p[1][2]*p[1][2]+p[1][3]*p[1][3]);
		p[2][0]=sqrt(mass[2]*mass[2]+p[2][1]*p[2][1]+p[2][2]*p[2][2]+p[2][3]*p[2][3]);
	}
	/* THREE-BODY DECAYS */
	else if(nbodies==3){
		kprimemax2=Misc::triangle(mass[0]-mass[3],mass[1],mass[2]);
		kprimemax=sqrt(kprimemax2);
		p3max=sqrt(Misc::triangle(mass[0],mass[1]+mass[2],mass[3]));
		e1max=sqrt(pow(mass[1],2)+p3max*p3max);
		e2max=sqrt(pow(mass[2],2)+p3max*p3max);
		e3max=sqrt(pow(mass[3],2)+p3max*p3max);
		//wmax=p3max*(e1max*e2max/(mass[1]*mass[2]))*(mass[1]+mass[2])/(e1max+e2max);
		wmax=p3max*pow(e1max+e2max,2)*e3max/(mass[1]+mass[2]);
		do{
			TRY_AGAIN:
			do{
				kprime[1]=kprimemax*(2.0*randy->ran()-1.0);
				kprime[2]=kprimemax*(2.0*randy->ran()-1.0);
				kprime[3]=kprimemax*(2.0*randy->ran()-1.0);
				kprimemag2=kprime[1]*kprime[1]+
					kprime[2]*kprime[2]+kprime[3]*kprime[3];
			} while(kprimemag2>kprimemax2);
			e1prime=sqrt(kprimemag2+mass[1]*mass[1]);
			e2prime=sqrt(kprimemag2+mass[2]*mass[2]);
			if(e1prime+e2prime+mass[3]>mass[0]) goto TRY_AGAIN;
			p3mag=sqrt(Misc::triangle(mass[0],e1prime+e2prime,mass[3]));
			cthet=1.0-2.0*randy->ran();
			sthet=sqrt(1.0-cthet*cthet);
			phi=2.0*PI*randy->ran();
			p[3][3]=p3mag*cthet;
			p[3][1]=p3mag*sthet*cos(phi);
			p[3][2]=p3mag*sthet*sin(phi);
			p[3][0]=sqrt(p3mag*p3mag+mass[3]*mass[3]);
			e12=sqrt(pow(e1prime+e2prime,2)+p3mag*p3mag);
			for(alpha=1;alpha<4;alpha++) u12[alpha]=-p[3][alpha]/(e1prime+e2prime);
			u12[0]=sqrt(1.0+u12[1]*u12[1]+u12[2]*u12[2]+u12[3]*u12[3]);
			kprime[0]=e1prime;
			Misc::lorentz(u12,kprime,p[1]);
			kprime[0]=e2prime;
			for(alpha=1;alpha<=3;alpha++) kprime[alpha]=-kprime[alpha];
			Misc::lorentz(u12,kprime,p[2]);
			//weight=p3mag*(p[1][0]*p[2][0]/(e1prime*e2prime))*((e1prime+e2prime)/(p[1][0]+p[2][0]));
			weight=p3mag*pow(p[1][0]+p[2][0],2)*p[3][0]/(e1prime+e2prime);
			if(weight/wmax>1.0){
				printf("3-body decay with weight/wmax>1, =%g\n",weight/wmax);
				printf("weight=%g, wmax=%g\n",weight,wmax);
				exit(1);
			}
			//printf("3-body w/wmax=%g\n",weight/wmax);
		} while(randy->ran()>weight/wmax);
	}
	/* FOUR-BODY DECAYS */
	else if(nbodies==4){
		kprimemax2=Misc::triangle(mass[0]-mass[3]-mass[4],mass[1],mass[2]);
		kprimemax=sqrt(kprimemax2);
		qprimemax2=Misc::triangle(mass[0]-mass[1]-mass[2],mass[3],mass[4]);
		qprimemax=sqrt(qprimemax2);
		
		ppmax=sqrt(Misc::triangle(mass[0],mass[1]+mass[2],mass[3]+mass[4]));
		e1max=sqrt(pow(mass[1],2)+ppmax*ppmax);
		e2max=sqrt(pow(mass[2],2)+ppmax*ppmax);
		e3max=sqrt(pow(mass[3],2)+ppmax*ppmax);
		e4max=sqrt(pow(mass[4],2)+ppmax*ppmax);
		wmax=ppmax*pow(e1max+e2max,2)*pow(e3max+e4max,2)/((mass[1]+mass[2])*(mass[3]+mass[4]));
		
		do{
			TRY_AGAIN_4:
			do{
				kprime[1]=kprimemax*(2.0*randy->ran()-1.0);
				kprime[2]=kprimemax*(2.0*randy->ran()-1.0);
				kprime[3]=kprimemax*(2.0*randy->ran()-1.0);
				kprimemag2=kprime[1]*kprime[1]+kprime[2]*kprime[2]+kprime[3]*kprime[3];
			} while(kprimemag2>kprimemax2);
			e1prime=sqrt(kprimemag2+mass[1]*mass[1]);
			e2prime=sqrt(kprimemag2+mass[2]*mass[2]);
			do{
				qprime[1]=qprimemax*(2.0*randy->ran()-1.0);
				qprime[2]=qprimemax*(2.0*randy->ran()-1.0);
				qprime[3]=qprimemax*(2.0*randy->ran()-1.0);
				qprimemag2=qprime[1]*qprime[1]+qprime[2]*qprime[2]+qprime[3]*qprime[3];
			} while(qprimemag2>qprimemax2);
			e3prime=sqrt(qprimemag2+mass[3]*mass[3]);
			e4prime=sqrt(qprimemag2+mass[4]*mass[4]);
			
			if(e1prime+e2prime+e3prime+e4prime>mass[0]) goto TRY_AGAIN_4;

			ppmag=Misc::triangle(mass[0],e1prime+e2prime,e3prime+e4prime);
			if(ppmag>0){
				ppmag=sqrt(ppmag);
				cthet=1.0-2.0*randy->ran();
				sthet=sqrt(1.0-cthet*cthet);
				phi=2.0*PI*randy->ran();
				pp[3]=ppmag*cthet;
				pp[1]=ppmag*sthet*cos(phi);
				pp[2]=ppmag*sthet*sin(phi);

				pp[0]=sqrt(ppmag*ppmag+pow(e1prime+e2prime,2));
				for(alpha=0;alpha<4;alpha++) u[alpha]=pp[alpha]/(e1prime+e2prime);
				kprime[0]=sqrt(mass[1]*mass[1]+kprimemag2);
				Misc::lorentz(u,kprime,p[1]);
				kprime[0]=sqrt(mass[2]*mass[2]+kprimemag2);
				for(alpha=1;alpha<4;alpha++) kprime[alpha]=-kprime[alpha];
				Misc::lorentz(u,kprime,p[2]);

				for(alpha=1;alpha<4;alpha++) pp[alpha]=-pp[alpha];
				pp[0]=sqrt(ppmag*ppmag+pow(e3prime+e4prime,2));
				for(alpha=0;alpha<4;alpha++) u[alpha]=pp[alpha]/(e3prime+e4prime);
				qprime[0]=sqrt(mass[3]*mass[3]+qprimemag2);
				Misc::lorentz(u,qprime,p[3]);
				qprime[0]=sqrt(mass[4]*mass[4]+qprimemag2);
				for(alpha=1;alpha<4;alpha++) qprime[alpha]=-qprime[alpha];
				Misc::lorentz(u,qprime,p[4]);

				weight=ppmag*pow(p[1][0]+p[2][0],2)*pow(p[3][0]+p[4][0],2)/((e1prime+e2prime)*(e3prime+e4prime));
				if(weight/wmax>1.0){
					printf("weight=%g, wmax=%g\n",weight,wmax);
					printf("four body decay weight > 1.0, =%g\n",weight/wmax);
					exit(1);
				}
			}
			else weight=0.0;
			if(weight/wmax>1.0){
				printf("4-body decay with weight/wmax>1, =%g\n",weight/wmax);
				exit(1);
			}
			/*
			printf("4-body w/wmax=%g\n",weight/wmax);
			printf("ppmag=%g, ppmax=%g\n",ppmag,ppmax);
			for(alpha=0;alpha<4;alpha++) ptot[alpha]=p[1][alpha]+p[2][alpha]+p[3][alpha]+p[4][alpha];
			printf("p1=(%g,%g,%g,%g), e1prime=%g\n",p[1][0],p[1][1],p[1][2],p[1][3],e1prime);
			printf("p2=(%g,%g,%g,%g), e2prime=%g\n",p[2][0],p[2][1],p[2][2],p[2][3],e2prime);
			printf("p3=(%g,%g,%g,%g), e3prime=%g\n",p[3][0],p[3][1],p[3][2],p[3][3],e3prime);
			printf("p4=(%g,%g,%g,%g), e4prime=%g\n",p[4][0],p[4][1],p[4][2],p[4][3],e4prime);
			printf("m_mother=%g, ptot=(%g,%g,%g,%g)\n",mass[0],ptot[0],ptot[1],ptot[2],ptot[3]);
			*/
		} while(randy->ran()>weight/wmax);
	}

	/* Boost the new particles */
	for(alpha=0;alpha<4;alpha++) u[alpha]=mother->p[alpha]/mother->GetMass();
	double pprime[4];
	for(ibody=0;ibody<nbodies;ibody++){
		dptr=daughter[ibody];
		Misc::lorentz(u,p[ibody+1],pprime);
		for(alpha=0;alpha<4;alpha++) dptr->p[alpha]=pprime[alpha];
		dptr->SetY();
		dptr->tau0=mother->tau0;
		for(alpha=0;alpha<4;alpha++) dptr->r[alpha]=mother->r[alpha];
		dptr->eta=mother->eta;
		if(dptr->eta>fabs(ETAMAX) && BJORKEN){
			printf("in decay.cc, eta out of bounds for daughter, =%g\n",dptr->eta);
			exit(1);
		}
		dptr->mass=dptr->resinfo->mass;
		dptr->Setp0();
		dptr->tau_lastint=daughter[ibody]->tau0;
		dptr->tau0=mother->tau0;
	}
}

#endif
