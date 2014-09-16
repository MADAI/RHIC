#ifndef __ANNIHILATE_CC__
#define __ANNIHILATE_CC__
#include "b3d.h"

int CB3D::Annihilate(CPart *part1,CPart *part2,int &ndaughters,array<CPart*,5> &daughter){
	CPart *dptr;
	int netq,netb,nets,nK0bar,nK0,nKplus,nKminus,npi0,npiplus,npiminus,npions,nkaons,npaircheck,qpions;
	FourVector *pa,*pb,pc,u,Ptot;
	double ma,mb,q,cthet,sthet,phi;
	bool bjtranslate=false;
	int idaughter,iK,ipair,alpha;
	double mt,Minv;
	double MM,P[4]={0.0},PP[4],T;
	const double g[4]={1.0,-1.0,-1.0,-1.0};
	CPartMap::iterator ppos;
	CB3DCell *newcell;
	ndaughters=0;
	if(BJORKEN && fabs(part1->eta)>ETAMAX)
		bjtranslate=true;

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
		// printf("Annihilate:: yikes: npi0<0, =%d\n",npi0);
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
	ndaughters=npi0+npiplus+npiminus+nKplus+nKminus+nK0+nK0bar;
	if(ndaughters != 5){
		printf("annihilation doesn't go to 5 particles\n");
		exit(1);
	}
	Minv=0.0;
	for(alpha=0;alpha<4;alpha++){
		P[alpha]=part1->p[alpha]+part2->p[alpha];
		Minv+=g[alpha]*P[alpha]*P[alpha];
	}
	Minv=sqrt(Minv);
	T=Minv; // T will be KE of emitted particles
	idaughter=0;
	while(nK0bar>0){
		daughter[idaughter]->resinfo=reslist->GetResInfoPtr(311);
		nK0bar-=1;
		T-=daughter[idaughter]->resinfo->mass;
		idaughter+=1;
	}
	while(nK0>0){
		daughter[idaughter]->resinfo=reslist->GetResInfoPtr(-311);
		nK0-=1;
		T-=daughter[idaughter]->resinfo->mass;
		idaughter+=1;
	}
	while(nKplus>0){
		daughter[idaughter]->resinfo=reslist->GetResInfoPtr(321);
		nKplus-=1;
		T-=daughter[idaughter]->resinfo->mass;
		idaughter+=1;
	}while(nKminus>0){
		daughter[idaughter]->resinfo=reslist->GetResInfoPtr(-321);
		nKminus-=1;
		T-=daughter[idaughter]->resinfo->mass;
		idaughter+=1;
	}
	while(npi0>0){
		daughter[idaughter]->resinfo=reslist->GetResInfoPtr(111);
		T-=daughter[idaughter]->resinfo->mass;
		npi0-=1;
		idaughter+=1;
	}
	while(npiplus>0){
		daughter[idaughter]->resinfo=reslist->GetResInfoPtr(211);
		npiplus-=1;
		T-=daughter[idaughter]->resinfo->mass;
		idaughter+=1;
	}
	while(npiminus>0){
		daughter[idaughter]->resinfo=reslist->GetResInfoPtr(-211);
		npiminus-=1;
		T-=daughter[idaughter]->resinfo->mass;
		idaughter+=1;
	}
	T=0.5*T/double(npions+nkaons); // Pick temperature of 0.5*KE/particles
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
	// PP is momentum of remaining two particles
	ma=daughter[ndaughters-2]->resinfo->mass;
	mb=daughter[ndaughters-1]->resinfo->mass;
	PP[0]=Minv+PP[0];
	if(PP[0]<ma+mb)
		goto DO_OVER;
	MM=0;
	for(alpha=0;alpha<4;alpha++)
		MM+=g[alpha]*PP[alpha]*PP[alpha];
	if(MM<(ma+mb)*(ma+mb))
		goto DO_OVER;
	MM=sqrt(MM);
	for(alpha=0;alpha<4;alpha++)
		u[alpha]=PP[alpha]/MM;
	pa=&daughter[ndaughters-2]->p;
	pb=&daughter[ndaughters-1]->p;
	cthet=1.0-2.0*randy->ran();
	sthet=sqrt(1.0-cthet*cthet);
	phi=2.0*PI*randy->ran();
	q=sqrt(Misc::triangle(MM,ma,mb));
	(*pa)[3]=q*cthet;
	(*pa)[1]=q*sthet*cos(phi);
	(*pa)[2]=q*sthet*sin(phi);
	(*pb)[3]=-(*pa)[3];
	(*pb)[2]=-(*pa)[2];
	(*pb)[1]=-(*pa)[1];
	(*pa)[0]=sqrt(ma*ma+q*q);
	(*pb)[0]=sqrt(mb*mb+q*q);
	Misc::Boost(u,*pa,pc);
	for(alpha=0;alpha<4;alpha++)
		(*pa)[alpha]=pc[alpha];
	Misc::Boost(u,*pb,pc);
	for(alpha=0;alpha<4;alpha++)
		(*pb)[alpha]=pc[alpha];

	//printf("Before: P=(%g,%g,%g,%g)\n",P[0],P[1],P[2],P[3]);
	for(alpha=0;alpha<4;alpha++){
		u[alpha]=P[alpha]/Minv;
		P[alpha]=0.0;
	}
	for(idaughter=0;idaughter<ndaughters;idaughter++){
		pa=&daughter[idaughter]->p;
		Misc::Boost(u,*pa,pc);
		for(alpha=0;alpha<4;alpha++){
			(*pa)[alpha]=pc[alpha];
			P[alpha]+=(*pa)[alpha];
		}
	}
	//printf("After: P=(%g,%g,%g,%g)\n",P[0],P[1],P[2],P[3]);
	Minv=0.0;
	for(alpha=0;alpha<4;alpha++)
		Minv+=g[alpha]*P[alpha]*P[alpha];
	Minv=sqrt(Minv);
	double rbar[4],etabar;

	rbar[1]=0.5*(part1->r[1]+part2->r[1]);
	rbar[2]=0.5*(part1->r[2]+part2->r[2]);
	etabar=0.5*(part1->eta+part2->eta);
	rbar[0]=tau*cosh(etabar);
	rbar[3]=tau*sinh(etabar);
	for(idaughter=0;idaughter<ndaughters;idaughter++){
		dptr=daughter[idaughter];
		dptr->tau0=tau;
		dptr->eta=etabar;
		for(alpha=0;alpha<4;alpha++)
			dptr->r[alpha]=rbar[alpha];
		dptr->active=true;
		dptr->SetMass();
		dptr->SetY();
		if(bjtranslate)
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
		if(dptr->p[0]<0.0){
			printf("dptr->p[0]=%g\n",dptr->p[0]);
			exit(1);
		}
	}
	return ndaughters;
}

double CB3D::GetAnnihilationSigma(CPart *part1,CPart *part2,double &vrel){
	double *p1=part1->p,*p2=part2->p;
	const double g[4]={1,-1,-1,-1};
	double Plab,sigma,p1dotp2,triangle,sigma_annihilation,rstrange;
	int alpha;
	part1->SetMass(); part2->SetMass();
	double m1squared=part1->msquared,m2squared=part2->msquared;
	p1dotp2=0.0;
	for(alpha=0;alpha<4;alpha++){
		p1dotp2+=part1->p[alpha]*part2->p[alpha]*g[alpha];
	}
	//Plab=sqrt((p1dotp2*p1dotp2/(part2->msquared))-part1->msquared);
	triangle=p1dotp2*p1dotp2-m1squared*m2squared;
	Plab=0.5*(m1squared+m2squared)*triangle/(m1squared*m2squared);
	Plab=sqrt(Plab);
	sigma_annihilation=6.7*pow(Plab/1000.0,-0.7)/double(NSAMPLE);
	rstrange=0.5*sqrt(sigma_annihilation);
	rstrange*=pow(ANNIHILATION_SREDUCTION,abs(part1->resinfo->strange))+pow(ANNIHILATION_SREDUCTION,abs(part2->resinfo->strange));
	sigma_annihilation=rstrange*rstrange;
	vrel=sqrt(triangle)/(p1[0]*p2[0]);
	//printf("sigma=%g,vrel=%g, p1=(%g,%g,%g,%g), p2=(%g,%g,%g,%g)\n",sigma_annihilation,vrel,p1[0],p1[1],p1[2],p1[3],p2[0],p2[1],p2[2],p2[3]);
	return sigma_annihilation;
}

#endif
