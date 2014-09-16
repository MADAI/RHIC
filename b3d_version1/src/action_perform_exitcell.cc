#ifndef __ACTION_PERFORM_COLLIDE_CC__
#define __ACTION_PERFORM_COLLIDE_CC__
#include "b3d.h"

void CAction::PerformExitCell(){
	CPart *part;
	CPartMap::iterator ppos;
	double mt;
	
	if(partmap.size()!=1){
		printf("FATAL: wrong number of particles in partmap for ExitCell, =%d\n",int(partmap.size()));
		exit(1);
	}
	ppos=partmap.begin();
	part=ppos->second;

	ppos=part->GetPos(&(part->cell->partmap));
	if(ppos==part->cell->partmap.end()){
		printf("YIKES, particle not in cell, tau=%g\n",tau);
		part->Print();
		part->cell->Print();
		exit(1);
	}


	part->Propagate(tau);
	
	double *r=part->r;
	double *p=part->p;
	
	/*
	 if(part->nextcell==NULL)
	 printf("EXITING (part %d): tau=%g, old cell: (%d,%d,%d), new cell = NULL\n",part->listid,tau, part->cell->ix,part->cell->iy,part->cell->ieta);
	 else
	 printf("EXITING (part %d): tau=%g, old cell: (%d,%d,%d), new cell: (%d,%d,%d)\n",part->listid,tau, part->cell->ix,part->cell->iy,part->cell->ieta, part->nextcell->ix,part->nextcell->iy,part->nextcell->ieta);
  */

	if(b3d->BJORKEN && part->cell->creflection!=NULL && part->nextcell==part->cell->creflection && fabs(fabs(part->eta)-b3d->ETAMAX)<1.0E-6){
		if(part->y<-b3d->ETAMAX){
			part->y+=2.0*b3d->ETAMAX;
			part->eta=b3d->ETAMAX;
		}
		else if(part->y>b3d->ETAMAX){
			part->y-=2.0*b3d->ETAMAX;
			part->eta=-b3d->ETAMAX;
		}
		else{
			printf("trying to reflect particle not moving fast enough!");
			Print();
			exit(1);
		}
		r[0]=tau*cosh(part->eta);
		r[3]=tau*sinh(part->eta);
		mt=sqrt(p[0]*p[0]-p[3]*p[3]);
		p[3]=mt*sinh(part->y);
		p[0]=mt*cosh(part->y);
	}
	part->KillActions();
	part->ChangeCell(part->nextcell);
	if(part->cell==NULL){
		part->ChangeMap(&(b3d->FinalPartMap));
		//printf("tau=%g, eta=%g, x=%g, y=%g\n",tau,part->eta,part->r[1],part->r[2]);
	}
	part->KillActions();
	part->FindActions();
	part->actionmother=b3d->nactions;
	b3d->nexit+=1;
	b3d->nactions+=1;
}

#endif
