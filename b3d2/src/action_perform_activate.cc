#ifndef __ACTION_PERFORM_ACTIVATE_CC__
#define __ACTION_PERFORM_ACTIVATE_CC__
#include "b3d.h"

void CAction::PerformActivate(){
	CPart *part;
	CPartMap::iterator ppos,pend;
	ppos=partmap.begin();
	part=ppos->second;
	
	if(b3d->COLLISIONS){
		part->ChangeMap(&(b3d->PartMap));
		part->ChangeCell(part->FindCell());
		part->CyclicReset();
	}
	else{
		part->ChangeMap(&(b3d->FinalPartMap));
		part->ChangeCell(NULL);
	}

	part->tau_lastint=tau;
	part->active=true;
	part->actionmother=b3d->nactions;
	part->FindActions();
	b3d->nactivate+=1;
}
#endif
