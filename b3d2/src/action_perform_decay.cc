#ifndef __ACTION_PERFORM_DECAY_CC__
#define __ACTION_PERFORM_DECAY_CC__
#include "b3d.h"

void CAction::PerformDecay(){
	CPart *mother,*dptr,*dptrb;
	CPartMap::iterator ppos;
	int ibody,nbodies,alpha;
	double mtot,mt,etamax=b3d->ETAMAX,mothermass;
	double deleta,Pcheck;
	FourVector Ptot;
	ppos=partmap.begin();
	mother=ppos->second;
	b3d->GetDeadParts(product);

	mothermass=mother->GetMass();
	if(mother->cell!=NULL && mother->cell!=mother->FindCell() && tau<b3d->TAUCOLLMAX){
		printf("Cells don't match for decaying mother\n");
		mother->cell->Print();
		mother->FindCell()->Print();
		mother->Print();
		exit(1);
	}
	if(tau>b3d->TAUCOLLMAX || mother->cell==NULL){
		while(b3d->BJORKEN && (mother->eta<-etamax || mother->eta>etamax)){
			if(mother->eta<-etamax) deleta=2.0*etamax*ceil((-etamax-mother->eta)/(2.0*etamax));
			if(mother->eta>etamax) deleta=-2.0*etamax*ceil((mother->eta-etamax)/(2.0*etamax));
			mother->eta+=deleta;
			mother->y+=deleta;
			mt=mother->GetMT();
			mother->p[0]=mt*cosh(mother->y);
			mother->p[3]=mt*sinh(mother->y);
			mother->r[0]=tau*cosh(mother->eta);
			mother->r[3]=tau*sinh(mother->eta);
		}
	}
	int ntry=0;
	do{
		mtot=0.0;
		if(ntry<25)
			mother->resinfo->DecayGetResInfoPtr(nbodies,daughterresinfo);
		else{
			mother->resinfo->DecayGetResInfoPtr_minmass(nbodies,daughterresinfo);
		}
		for(ibody=0;ibody<nbodies;ibody++){
			mtot+=daughterresinfo[ibody]->mass;
		}
		if(ntry>25){
			printf("FATAL: action_perform_decay, ntry too big, mothermass=%g\n",mother->GetMass());
			mother->Print();
			exit(1);
		}
		ntry++;
	}while(mtot>mothermass);
	for(ibody=0;ibody<nbodies;ibody++){
		product[ibody]->resinfo=daughterresinfo[ibody];
	}
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for(alpha=0;alpha<4;alpha++)
		Ptot[alpha]=mother->p[alpha];
	
	b3d->Decay(mother,nbodies,product);
	
	for(alpha=0;alpha<4;alpha++){
		for(ibody=0;ibody<nbodies;ibody++){
			Ptot[alpha]-=product[ibody]->p[alpha];
		}
		Pcheck+=fabs(Ptot[alpha]);
	}
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	mother->actionmother=b3d->nactions;
	int jbody;
	for(ibody=0;ibody<nbodies;ibody++){
		dptr=product[ibody];
		dptr->active=true;
		dptr->weight=mother->weight;
		dptr->reality=mother->reality;
		dptr->balanceID=mother->balanceID;
		if(b3d->BALANCE_DECAY){
			if(mother->reality){
				dptrb=b3d->GetDeadPart();
				dptrb->Copy(dptr);
				dptrb->balanceID=b3d->ibalmax+1;
				b3d->SplitPart(dptr,dptrb);
				dptrb->nscatt=0;
				dptrb->tau_lastint=tau;
				dptrb->actionmother=b3d->nactions;
				if(dptrb->currentmap!=mother->currentmap)
					dptrb->ChangeMap(mother->currentmap);
				dptrb->ChangeCell(dptrb->FindCell());
				dptrb->reality=false;
				dptrb->weight=1;
				dptrb->FindActions();
			}
		}
		if(dptr->reality){
			dptr->nscatt=0;
		}
		else{
			dptr->nscatt=mother->nscatt;
		}
		dptr->tau_lastint=tau;
		dptr->actionmother=b3d->nactions;
		dptr->ChangeCell(dptr->FindCell());
		if(dptr->currentmap!=mother->currentmap)
			dptr->ChangeMap(mother->currentmap);
		dptr->FindActions();
	}
	if(b3d->BALANCE_DECAY && mother->reality)
		b3d->ibalmax+=1;
	mother->Kill();
	b3d->ndecay+=1;
}
#endif
