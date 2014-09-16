#ifndef __ACTION_PERFORM_COLLIDE_CC__
#define __ACTION_PERFORM_COLLIDE_CC__
#include "b3d.h"

void CAction::PerformCollide(){
	int colltype,iproduct,nproducts,wproduct,nrealscatt=0,nfakescatt=0;
	bool productreality=true;
	CPart *part1,*part2,*part,*splitpart;
	CPartMap::iterator ppos;
	CB3DCell *cell;

	ppos=partmap.begin();
	part1=ppos->second;
	++ppos;
	part2=ppos->second;
	part1->actionmother=b3d->nactions;
	part2->actionmother=b3d->nactions;
	
	b3d->GetDeadParts(product);
	colltype=b3d->Collide(part1,part2,nproducts,product,pibsquared);
	if(colltype==0 || nproducts==0){
		b3d->npass+=1;
		return;
	}
	if(colltype==1){
		b3d->nmerge+=1;
	}
	if(colltype==2)
		b3d->nscatter+=1;
	if(colltype==3)
		b3d->ninelastic+=1;
	if(colltype==4)
		b3d->nannihilate+=1;
	
	wproduct=part1->weight*part2->weight;
	if(part1->reality && part2->reality){
		part1->Kill();
		part2->Kill();
	}
	else{
		productreality=false;
		if(part2->reality){
			splitpart=part1;
			part1=part2;
			part2=splitpart;
		}
		nfakescatt=part2->nscatt;
		nrealscatt=part1->nscatt;
		b3d->SplitPart(part1,part2);
		part1->tau_lastint=tau;
		part2->weight=-part2->weight;
		part2->tau_lastint=tau;
		part1->FindActions();
		part2->nscatt=nfakescatt+1;
		part2->FindActions();
	}
		
	for(iproduct=0;iproduct<nproducts;iproduct++){
		part=product[iproduct];
		part->reality=productreality;
		if(productreality)
			part->nscatt=0;
		else
			part->nscatt=nfakescatt+1;
		part->SetMass();
		part->active=true;
		part->weight=wproduct;
		part->tau_lastint=tau;
		part->actionmother=b3d->nactions;
		cell=part->FindCell();
		part->ChangeCell(cell);
		if(part->cell!=NULL){
			if(part->currentmap!=&b3d->PartMap)
				part->ChangeMap(&b3d->PartMap);
		}
		else{
			if(part->currentmap!=&b3d->FinalPartMap)
				part->ChangeMap(&b3d->FinalPartMap);
		}
		part->FindActions();
	}
}

#endif
