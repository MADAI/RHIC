#ifndef __ACTION_PERFORM_COLLIDE_CC__
#define __ACTION_PERFORM_COLLIDE_CC__
#include "b3d.h"

void CAction::PerformCollide(){
	if(b3d->BALANCE_CALC || b3d->BALANCE_DECAY){
		PerformCollide_BALANCE();
		return;
	}
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
	if((part1->reality && part2->reality)){
		part1->Kill();
		part2->Kill();
	}
	else{
		productreality=false;
		if(part2->reality){ //make sure part 2 is the fake part and part 1 is the real part
			splitpart=part1;
			part1=part2;
			part2=splitpart;
		}
		nfakescatt=part2->nscatt;
		nrealscatt=part1->nscatt;
		b3d->SplitPart(part1,part2);

		part1->tau_lastint=tau;
		part1->FindActions();

		part2->weight=-part1->weight;
		part2->tau_lastint=tau;
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

void CAction::PerformCollide_BALANCE(){
	int colltype,iproduct,nproducts,wproduct,nrealscatt=0,nfakescatt=0;
	bool productreality=true;
	CPart *part1,*part1a,*part2,*part2a,*part,*splitpart;
	CPartMap::iterator ppos;
	CB3DCell *cell;

	ppos=partmap.begin();
	part1=ppos->second;
	++ppos;
	part2=ppos->second;
	int balanceID=0;
	if(part1->reality==false){
		balanceID=part1->balanceID;
		productreality=false;
	}
	if(part2->reality==false){
		balanceID=part2->balanceID;
		productreality=false;
	}
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
	if(colltype==2){
		part1->Kill();
		part2->Kill();
	}
	else{
		if(part1->reality&&part2->reality){
			b3d->ibalmax+=1;
			balanceID=b3d->ibalmax;
		}
		if(part1->reality){
			part1a=b3d->GetDeadPart();
			b3d->SplitPart(part1,part1a);
			
			part1->tau_lastint=tau;
			part1->nscatt=0;
			part1->balanceID=0;
			part1->FindActions();
			
			part1a->weight=-part1->weight;
			part1a->tau_lastint=tau;
			part1a->nscatt=nfakescatt+1;
			part1a->balanceID=balanceID;
			part1a->FindActions();
		}
		if(part2->reality){
			part2a=b3d->GetDeadPart();
			b3d->SplitPart(part2,part2a);
			
			part2->tau_lastint=tau;
			part2->nscatt=0;
			part2->balanceID=0;
			part2->FindActions();
			
			part2a->weight=-part2->weight;
			part2a->tau_lastint=tau;
			part2a->nscatt=nfakescatt+1;
			part2a->balanceID=balanceID;
			part2a->FindActions();
		}
	}
	for(iproduct=0;iproduct<nproducts;iproduct++){
		part=product[iproduct];
		part->reality=productreality;
		part->balanceID=balanceID;
		if((b3d->BALANCE_CALC && colltype==2) || (b3d->BALANCE_DECAY && colltype==2)){
			if(iproduct==0 && part1->reality==true){
				part->reality=true;
				part->balanceID=0;
			}
			if(iproduct==1 && part2->reality==true){
				part->reality=true;
				part->balanceID=0;
			}
		}
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