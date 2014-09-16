#ifndef __ANNIHILATE_CC__
#define __ANNIHILATE_CC__
#include "b3d.h"
using namespace std;

int CB3D::Annihilate(CPart *part1,CPart *part2){
	CPart *daughter[5],*dptr;
	CResInfo daughterresinfo[5];
	int netq,netb,nets,nK0bar,nK0,nKplus,nKminus,npi0,npiplus,npiminus,npions,nkaons,npaircheck,qpions;
	double *pa,*pb,pc[4],ma,mb,q,cthet,sthet,phi;
	double u[4];
	int idaughter,ndaughters=0,iK,ipair,alpha;
	double mt,Minv,Ptot[4];
	double MM,P[4]={0.0},PP[4],T;
	const double g[4]={1.0,-1.0,-1.0,-1.0};
	CPartMap::iterator ppos;
	CB3DCell *newcell;

	netq = part1->resinfo->charge+part2->resinfo->charge;
	netb = part1->resinfo->baryon+part2->resinfo->baryon;
	nets = part1->resinfo->strange+part2->resinfo->strange;
	nkaons=abs(nets);
	RECHARGE:
	nK0bar=nK0=nKplus=nKminus=npiplus=npiminus=npi0=0;
	for(iK=0;iK<nkaons;iK++){
		if(randy->ran()>0.5){
			if(nets>0)
				nK0+=1;
			else
				nK0bar+=1;
		}
		else{
			if(nets>0)
				nKplus+=1;
			else
				nKminus+=1;
		}
	}
	qpions=netq-nKplus+nKminus;
	if(qpions>0)
		npiplus+=qpions;
	else
		npiminus-=qpions;
	npions=5-nkaons;
	npi0=npions-npiplus-npiminus;
	if(netq!= nKplus+npiplus-nKminus-npiminus){
		printf("charges don't add up\n");
		exit(1);
	}
	if(nets!= nKplus+nK0-nKminus-nK0bar){
		printf("charges don't add up\n");
		exit(1);
	}

	if(npi0<0){
		printf("yikes: npi0<0, =%d\n",npi0);
		goto RECHARGE;
	}
	npaircheck=0;
	if(npi0>=2)
		npaircheck=1;
	if(npi0>=4)
		npaircheck=2;
	for(ipair=0;ipair<npaircheck;ipair++){
		if(randy->ran()<0.66666666667){
			npiplus+=1;
			npiminus+=1;
			npi0-=2;
		}
	}
	//printf("npi=(%d,%d,%d), nK=(%d,%d,%d,%d)\n",npi0,npiplus,npiminus,nK0,nK0bar,nKplus,nKminus);
	if(npi0+npiplus+npiminus+nKplus+nKminus+nK0+nK0bar != 5){
		printf("annihilation doesn't go to 5 particles\n");
		exit(1);
	}
	Minv=0.0;
	ndaughters=npi0+npiplus+npiminus+nKplus+nKminus+nK0+nK0bar;
	ppos=DeadPartMap.begin();
	for(idaughter=0;idaughter<ndaughters;idaughter++){
		daughter[idaughter]=ppos->second;
		daughter[idaughter]->actionmap.clear();
		++ppos;
	}
	for(alpha=0;alpha<4;alpha++){
		P[alpha]=part1->p[alpha]+part2->p[alpha];
		Minv+=g[alpha]*P[alpha]*P[alpha];
	}
	Minv=sqrt(Minv);
	T=Minv; // T will be KE of emitted particles
	idaughter=0;
	while(nK0bar>0){
		reslist->GetResInfoptr(311,daughter[idaughter]->resinfo);
		nK0bar-=1;
		T-=daughter[idaughter]->resinfo->mass;
		idaughter+=1;
	}
	while(nK0>0){
		reslist->GetResInfoptr(-311,daughter[idaughter]->resinfo);
		nK0-=1;
		T-=daughter[idaughter]->resinfo->mass;
		idaughter+=1;
	}
	while(nKplus>0){
		reslist->GetResInfoptr(321,daughter[idaughter]->resinfo);
		nKplus-=1;
		T-=daughter[idaughter]->resinfo->mass;
		idaughter+=1;
	}while(nKminus>0){
		reslist->GetResInfoptr(-321,daughter[idaughter]->resinfo);
		nKminus-=1;
		T-=daughter[idaughter]->resinfo->mass;
		idaughter+=1;
	}
	while(npi0>0){
		reslist->GetResInfoptr(111,daughter[idaughter]->resinfo);
		T-=daughter[idaughter]->resinfo->mass;
		npi0-=1;
		idaughter+=1;
	}
	while(npiplus>0){
		reslist->GetResInfoptr(211,daughter[idaughter]->resinfo);
		npiplus-=1;
		T-=daughter[idaughter]->resinfo->mass;
		idaughter+=1;
	}
	while(npiminus>0){
		reslist->GetResInfoptr(-211,daughter[idaughter]->resinfo);
		npiminus-=1;
		T-=daughter[idaughter]->resinfo->mass;
		idaughter+=1;
	}
	T=0.4*T/double(npions+nkaons); // Pick temperature of 0.4*KE/particles
	//printf("T=%g, idaughter=%d=?%d\n",T,idaughter,ndaughters);
	//printf("BEFORE: Minv=%g, P=(%g,%g,%g,%g)\n",Minv,P[0],P[1],P[2],P[3]);
	
	// Do a 2-body decay for the last 2 particles
DO_OVER:
	//printf("DO_OVER\n");
	for(alpha=0;alpha<4;alpha++)
		PP[alpha]=0;
	for(idaughter=0;idaughter<ndaughters-2;idaughter++){
		dptr=daughter[idaughter];
		randy->generate_boltzmann(dptr->resinfo->mass,T,dptr->p);
		for(alpha=0;alpha<4;alpha++){
			PP[alpha]-=dptr->p[alpha];
		}
	}
	PP[0]=Minv+PP[0];
	MM=0;
	for(alpha=0;alpha<4;alpha++)
		MM+=g[alpha]*PP[alpha]*PP[alpha];
	ma=daughter[ndaughters-2]->resinfo->mass;
	mb=daughter[ndaughters-1]->resinfo->mass;
	if(MM<(ma+mb)*(ma+mb))
		goto DO_OVER;
	MM=sqrt(MM);
	if(MM<ma+mb)
		goto DO_OVER;
	//printf("PP=(%g,%g,%g,%g)\n",PP[0],PP[1],PP[2],PP[3]);
	for(alpha=0;alpha<4;alpha++)
		u[alpha]=PP[alpha]/MM;
	pa=daughter[ndaughters-2]->p;
	pb=daughter[ndaughters-1]->p;
	cthet=1.0-2.0*randy->ran();
	sthet=sqrt(1.0-cthet*cthet);
	phi=2.0*PI*randy->ran();
	q=sqrt(Misc::triangle(MM,ma,mb));
	//printf("ma=%g,mb=%g,q=%g\n",ma,mb,q);
	pa[3]=q*cthet;
	pa[1]=q*sthet*cos(phi);
	pa[2]=q*sthet*sin(phi);
	pb[3]=-pa[3];
	pb[2]=-pa[2];
	pb[1]=-pa[1];
	pa[0]=sqrt(ma*ma+q*q);
	pb[0]=sqrt(mb*mb+q*q);
	//printf("pa=(%g,%g,%g,%g)\n",pa[0],pa[1],pa[2],pa[3]);
	//printf("pb=(%g,%g,%g,%g)\n",pb[0],pb[1],pb[2],pb[3]);
	Misc::Boost(u,pa,pc);
	for(alpha=0;alpha<4;alpha++)
		pa[alpha]=pc[alpha];
	Misc::Boost(u,pb,pc);
	for(alpha=0;alpha<4;alpha++)
		pb[alpha]=pc[alpha];
	//printf("pa=(%g,%g,%g,%g)\n",pa[0],pa[1],pa[2],pa[3]);
	//printf("pb=(%g,%g,%g,%g)\n",pb[0],pb[1],pb[2],pb[3]);
	//printf("pa+pb=(%g,%g,%g,%g)\n",pa[0]+pb[0],pa[1]+pb[1],pa[2]+pb[2],pa[3]+pb[3]);

	for(alpha=0;alpha<4;alpha++){
		u[alpha]=P[alpha]/Minv;
		P[alpha]=0.0;
	}
	for(idaughter=0;idaughter<ndaughters;idaughter++){
		pa=daughter[idaughter]->p;
		Misc::Boost(u,pa,pc);
		for(alpha=0;alpha<4;alpha++){
			pa[alpha]=pc[alpha];
			P[alpha]+=pa[alpha];
		}
	}
	Minv=0.0;
	for(alpha=0;alpha<4;alpha++)
		Minv+=g[alpha]*P[alpha]*P[alpha];
	Minv=sqrt(Minv);
	//printf("AFTER: Minv=%g, P=(%g,%g,%g,%g)\n",Minv,P[0],P[1],P[2],P[3]);
	//printf("-------------------------------------------------\n");

	double rbar[4],etabar;

	rbar[1]=0.5*(part1->r[1]+part2->r[1]);
	rbar[2]=0.5*(part1->r[2]+part2->r[2]);
	etabar=0.5*(part1->eta+part2->eta);
	rbar[0]=tau*cosh(etabar);
	rbar[3]=tau*sinh(etabar);
	for(idaughter=0;idaughter<ndaughters;idaughter++){
		dptr=daughter[idaughter];
		dptr->actionmap.clear();
		dptr->tau0=tau;
		dptr->eta=etabar;
		for(alpha=0;alpha<4;alpha++)
			dptr->r[alpha]=rbar[alpha];
		dptr->y=atanh(dptr->p[3]/dptr->p[0]);
		dptr->tau_lastint=dptr->tau0;
		dptr->actionmother=nactions;
		dptr->active=true;
		dptr->CyclicReset();
		newcell=dptr->FindCell();
		dptr->ChangeCell(newcell);
		if(newcell==NULL){
			dptr->ChangeMap(&FinalPartMap);
		}
		else{
			dptr->ChangeMap(&PartMap);
		}
		if(fabs(dptr->eta)>ETAMAX){
			printf("out of range\n");
			exit(1);
		}
	}
	part1->Kill();
	part2->Kill();
	for(idaughter=0;idaughter<ndaughters;idaughter++){
		daughter[idaughter]->FindActions();
	}

	return ndaughters;
	
}

#endif
