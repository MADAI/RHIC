#ifndef __ACTION_PERFORM_COLLIDE_CC__
#define __ACTION_PERFORM_COLLIDE_CC__
#include "b3d.h"

void CAction::PerformExitCell(){
	CPart *part;
	CPartMap::iterator ppos;
	double mt;
	
	ppos=partmap.begin();
	part=ppos->second;

	ppos=part->GetPos(&(part->cell->partmap));
	
	FourVector *r=&part->r;
	FourVector *p=&part->p;

	if(b3d->BJORKEN && part->cell->creflection!=NULL && part->nextcell==part->cell->creflection && fabs(fabs(part->eta)-b3d->ETAMAX)<1.0E-6){
		if(part->y<0){
			part->y+=2.0*b3d->ETAMAX;
			part->eta=b3d->ETAMAX;
		}
		else{
			part->y-=2.0*b3d->ETAMAX;
			part->eta=-b3d->ETAMAX;
		}
		(*r)[0]=tau*cosh(part->eta);
		(*r)[3]=tau*sinh(part->eta);
		mt=sqrt((*p)[0]*(*p)[0]-(*p)[3]*(*p)[3]);
		(*p)[3]=mt*sinh(part->y);
		(*p)[0]=mt*cosh(part->y);
	}
	part->ChangeCell(part->nextcell);
	if(part->cell==NULL){
		part->ChangeMap(&(b3d->FinalPartMap));
	}
	part->FindActions();
	b3d->nexit+=1;
}

#endif
