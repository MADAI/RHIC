#ifndef __ACTION_PERFORM_DECAY_CC__
#define __ACTION_PERFORM_DECAY_CC__
#include "b3d.h"

void CAction::PerformDecay(){
	CPart *mother;
	CPart **daughter=new CPart *[5],*dptr;
	CResInfo **daughterresinfo=new CResInfo *[5];
	CPartMap::iterator ppos,p1,plast,pend;
	int ibody,nbodies;
	double mtot,mt,etamax=b3d->ETAMAX,mothermass;
	double deleta;
	bool rekey=false;
	CPartMap::iterator pos,nextpos;
	if(fabs(b3d->tau-tau)>1.0E-3){
		printf("taus do not match in PerformDecay, b3d->tau=%g, action->tau=%g\n",b3d->tau,tau);
		exit(1);
	}
	b3d->ndecay+=1;
	ppos=partmap.begin();
	mother=ppos->second;
	mothermass=mother->GetMass();
	if(mothermass<mother->resinfo->minmass){
		printf("FATAL: decaying mother has mass below minimum\n");
		cout << "minimum mass: " <<  mother->resinfo->minmass << "\tmass: " << mothermass << endl;
		mother->Print();
		mother->resinfo->Print();
		exit(1);
	}
	if(!mother->active){
		printf("FATAL: particle should be active to decay!!!!!!!!!!\n");
		Print();
		exit(1);
	}
	mother->Propagate(tau);
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
		if(ntry<25) mother->resinfo->DecayGetResInfoptr(nbodies,daughterresinfo);
		else mother->resinfo->DecayGetResInfoptr_minmass(nbodies,daughterresinfo);
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
	
	ppos=b3d->DeadPartMap.begin();
	for(ibody=0;ibody<nbodies;ibody++){
		daughter[ibody]=ppos->second;
		dptr=daughter[ibody];
		dptr->resinfo=daughterresinfo[ibody];
		++ppos;
	}
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	b3d->Decay(mother,nbodies,daughter);
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if(nbodies>4){
		printf("FATAL: Action.Perform(), In decay, nbodies=%d\n",nbodies);
		exit(1);
	}
	if(b3d->DeadPartMap.size()<nbodies){
		printf("PERFORM_DECAY: MUST INCREASE NPARTMAX\n");
		exit(1);
	}
	CB3DCell *newcell;
	for(ibody=0;ibody<nbodies;ibody++){
		dptr=daughter[ibody];
		dptr->active=true;
		dptr->actionmother=b3d->nactions;
		dptr->actionmap.clear();
		newcell=dptr->FindCell();
		dptr->ChangeCell(newcell);
		if(newcell==NULL){
			dptr->ChangeMap(&b3d->FinalPartMap);
		}
		else{
			dptr->ChangeMap(&b3d->PartMap);
		}
	}
	mother->Kill();
	for(ibody=0;ibody<nbodies;ibody++){
		daughter[ibody]->FindActions();
	}

	b3d->nactions+=1;
	b3d->ndecay+=1;

	delete [] daughterresinfo;
	delete [] daughter;
}
#endif
