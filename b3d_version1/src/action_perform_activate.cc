#ifndef __ACTION_PERFORM_ACTIVATE_CC__
#define __ACTION_PERFORM_ACTIVATE_CC__
#include "b3d.h"

void CAction::PerformActivate(){
	CPart *part;
	CPartMap::iterator ppos,pend;
	ppos=partmap.begin();
	part=ppos->second;
	part->Propagate(tau);
	double mt;
	part->active=true;
	part->tau0=tau;
	part->tau_lastint=tau;
	part->actionmother=b3d->nactions;
	b3d->nactions++;
	if(part->currentmap!=&(b3d->PartMap)){
		printf("FATAL: particles to be activated should be in PartMap\n");
		part->Print();
		exit(1);
	}
	part->KillActions();
	if(part->currentmap!=&(b3d->PartMap)){
		printf("A FATAL: particles to be activated should be in PartMap\n");
		part->Print();
		exit(1);
	}
	part->ChangeMap(&(b3d->PartMap));
	part->CyclicReset();
	part->ChangeCell(part->FindCell());
	part->FindActions();
	b3d->nactivate+=1;
}
#endif
