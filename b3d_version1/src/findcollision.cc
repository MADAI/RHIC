#ifndef __FINDCOLLISION_CC__
#define __FINDCOLLISION_CC__

#include "b3d.h"
using namespace std;

//this assumes part1 is inside (lower r0) than part2
bool CB3D::FindCollision(CPart *part1,CPart *part2,double &taucoll){
	double denom;
	bool collide=false;
	CB3DCell *cell1=part1->cell;
	CB3DCell *cell2=part2->cell;


	if(part1->active==false || part2->active==false || part1->actionmother==part2->actionmother || cell1==NULL || cell2==NULL){
		printf("FindCollision:: Why am I here?\n");
		part1->Print();
		part2->Print();
		exit(1);
		return false;
	}
	double q[4],P[4],r[4],p1dotp2=0.0,p1dotr=0.0,p2dotr=0.0,p1squared=0.0,p2squared=0.0,rsquared=0.0;
	double tau1,tau2,eta1,y1,mt,t1,t2,z1,z2,x1,x2,y2,u[4];
	double *p1=part1->p,*p2=part2->p,*r1=part1->r,*r2=part2->r;
	double pibsquared;
	const int g[4]={1,-1,-1,-1};
	int alpha;
	bool flip=false;
	if(BJORKEN && ((cell1->ieta==0 && cell2->ieta==2*NETA-1) || (cell1->ieta==2*NETA-1 && cell2->ieta==0))){
		flip=true;
		p1=new double[4];
		r1=new double[4];
		for(alpha=0;alpha<4;alpha++){
			p1[alpha]=part1->p[alpha]; r1[alpha]=part1->r[alpha];
		}
		tau1=part1->tau0;
		if(cell1->ieta==0){
			eta1=part1->eta+2.0*ETAMAX; y1=part1->y+2.0*ETAMAX;
		}
		else{
			eta1=part1->eta-2.0*ETAMAX; y1=part1->y-2.0*ETAMAX;
		}
		r1[0]=tau1*cosh(eta1);
		r1[3]=tau1*sinh(eta1);
		mt=part1->resinfo->mass;
		mt=sqrt(mt*mt+p1[1]*p1[1]+p1[2]*p1[2]);
		p1[0]=mt*cosh(y1);
		p1[3]=mt*sinh(y1);
	}

	for(alpha=0;alpha<4;alpha++){
		r[alpha]=r1[alpha]-r2[alpha];
		rsquared+=g[alpha]*r[alpha]*r[alpha];
		p1dotp2+=g[alpha]*p1[alpha]*p2[alpha];
		p1dotr+=g[alpha]*p1[alpha]*r[alpha];
		p2dotr+=g[alpha]*p2[alpha]*r[alpha];
		p1squared+=g[alpha]*p1[alpha]*p1[alpha];
		p2squared+=g[alpha]*p2[alpha]*p2[alpha];		
	}
	denom=p1dotp2*p1dotp2-p1squared*p2squared;
	t1=p1[0]*(p1dotr*p2squared-p2dotr*p1dotp2)/denom;
	t2=-p2[0]*(p2dotr*p1squared-p1dotr*p1dotp2)/denom;
		//double tt=t1; t1=t2; t2=tt;
		//printf("t1=%g, r1[0]=%g, t2=%g, r2[0]=%g\n",t1,r1[0],t2,r2[0]);
	if(t1+r1[0]>0 && t2+r2[0]>0){
		z1=r1[3]+(p1[3]/p1[0])*t1;
		z2=r2[3]+(p2[3]/p2[0])*t2;
		x1=r1[1]+(p1[1]/p1[0])*t1;
		x2=r2[1]+(p2[1]/p2[0])*t2;
		y1=r1[2]+(p1[2]/p1[0])*t1;
		y2=r2[2]+(p2[2]/p2[0])*t2;
		t1+=r1[0];
		t2+=r2[0];
		if(fabs(z1)<t1 && fabs(z2)<t2){
			pibsquared=PI*(-rsquared+(2.0*p1dotp2*p1dotr*p2dotr-p1dotr*p1dotr*p2squared-p2dotr*p2dotr*p1squared)/denom);
			/*
			if(fabs(pibsquared-GetPiBsquared(part1,part2))>1.0E-4){
			printf("failed to get pibsquared, pibsquared=%g != %\n",pibsquared,-GetPiBsquared(part1,part2));
			}
			*/
			if(pibsquared<SIGMAMAX){
				tau1=sqrt(t1*t1-z1*z1);
				tau2=sqrt(t2*t2-z2*z2);
				taucoll=0.5*(tau1+tau2);
		//printf("tau1=%g,tau2=%g, taucoll=%g\n",tau1,tau2,taucoll);
				if(taucoll>tau && taucoll<part1->tauexit && taucoll<part2->tauexit && taucoll<TAUCOLLMAX){
					collide=true;
					//printf("taucoll=%13.7e\n",taucoll);
				}
				//if((part1->cell->ieta==0 && part2->cell->ieta==NETA-1) || (part2->cell->ieta==0 && part1->cell->ieta==NETA-1)){
					//printf("CB3D:FindCollision, cells different, taucoll=%g\n",taucoll);
					//part1->cell->Print(); part2->cell->Print();
				//}
			}
		}
	}
	if(flip){
		delete [] r1; delete [] p1;
	}


	if(collide==true){
		//CB3DCell *cell1=part1->cell;
		//CB3DCell *cell2=part2->cell;
		//printf("found a collision!!!!\n");
		//Misc::Pause();
		if(taucoll<part1->tauexit && taucoll<part2->tauexit){
			AddAction_Collision(part1,part2,taucoll);
		}
	}
	
	return collide;
}

#endif
