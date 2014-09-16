#ifndef __FINDCOLLISION_CC__
#define __FINDCOLLISION_CC__

#include "b3d.h"

bool CB3D::FindCollision(CPart *part1,CPart *part2,double &taucoll){
	if(!part1->reality && !part2->reality)
		return false;
	double denom,sigmamax;
	CResInfo *resinfo1,*resinfo2;
	bool collide=false;
	CB3DCell *cell1=part1->cell;
	CB3DCell *cell2=part2->cell;
	double p1dotp2=0.0,p1dotr=0.0,p2dotr=0.0,m1squared,m2squared,rsquared=0.0;
	double tau1,tau2,eta1,y1,mt,t1,t2,z1,z2,x1,x2,y2;
	FourVector u,p1,p2,r1,r2,q,P,r;
	double pibsquared;
	const int g[4]={1,-1,-1,-1};
	int alpha;
	for(alpha=0;alpha<4;alpha++){
		r1[alpha]=part1->r[alpha];
		r2[alpha]=part2->r[alpha];
		p1[alpha]=part1->p[alpha];
		p2[alpha]=part2->p[alpha];
	}
	m1squared=part1->msquared;
	m2squared=part2->msquared;
	bool flip=false; // check for collisions across cyclic boundary
	if(BJORKEN && ((cell1->ieta==0 && cell2->ieta==2*NETA-1) || (cell1->ieta==2*NETA-1 && cell2->ieta==0))){
		flip=true;
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
	}
		
	denom=p1dotp2*p1dotp2-m1squared*m2squared;
	pibsquared=PI*(-rsquared+(2.0*p1dotp2*p1dotr*p2dotr-p1dotr*p1dotr*m2squared-p2dotr*p2dotr*m1squared)/denom);
	resinfo1=part1->resinfo;
	resinfo2=part2->resinfo;
	sigmamax=SIGMADEFAULT+reslist->SigmaMaxArray[resinfo1->ires][resinfo2->ires]+SIGMAINELASTIC;
	if(BARYON_ANNIHILATION && resinfo1->baryon*resinfo2->baryon<0)
		sigmamax+=30.0; // cuts off annihilation Xsection at 300 mb
	if(pibsquared<sigmamax/double(NSAMPLE)){
		
		t1=p1[0]*(p1dotr*m2squared-p2dotr*p1dotp2)/denom;
		t2=-p2[0]*(p2dotr*m1squared-p1dotr*p1dotp2)/denom;
		
		if(t1+r1[0]>0 && t2+r2[0]>0 && (t1+t2)>0.0){
			z1=r1[3]+(p1[3]/p1[0])*t1;
			z2=r2[3]+(p2[3]/p2[0])*t2;
			x1=r1[1]+(p1[1]/p1[0])*t1;
			x2=r2[1]+(p2[1]/p2[0])*t2;
			y1=r1[2]+(p1[2]/p1[0])*t1;
			y2=r2[2]+(p2[2]/p2[0])*t2;
			t1+=r1[0];
			t2+=r2[0];
			if(fabs(z1)<t1 && fabs(z2)<t2){
				tau1=sqrt(t1*t1-z1*z1);
				tau2=sqrt(t2*t2-z2*z2);
				taucoll=0.5*(tau1+tau2);
				if(taucoll>tau && taucoll<part1->tauexit && taucoll<part2->tauexit && taucoll<TAUCOLLMAX){
					AddAction_Collision(part1,part2,taucoll,pibsquared);
					collide=true;
				}
			}
		}
	}
	return collide;
}

#endif
