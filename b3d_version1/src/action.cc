#ifndef __ACTION_CC__
#define __ACTION_CC__

#include "b3d.h"

CB3D *CAction::b3d=NULL;
CAction::CAction(){
}

// type=0(creation) 1(decay) 2(collision) 3(VizWrite) 4(DensCalc)

CActionMap::iterator CAction::GetPos(CActionMap *emap){
	pair<CActionMap::iterator,CActionMap::iterator> epospair;
	CActionMap::iterator epos;
	epospair=emap->equal_range(key);
	epos=epospair.first;
	while(epos->second!=this && epos!=epospair.second){
		++epos;
	}
	if(epos->second!=this){
		return emap->end();
	}
	else return epos;
}

void CAction::MoveToActionMap(){
	CActionMap::iterator epos,eepos;
	if(currentmap==&b3d->ActionMap){
		printf("trying to move action to ActionMap even though action is already in ActionMap\n");
		printf("wrong current map\n");
		exit(1);
	}
	if(currentmap==&b3d->DeadActionMap){
		epos=GetPos(currentmap);
		if(epos!=currentmap->end()){
			b3d->DeadActionMap.erase(epos);
			partmap.clear();
		}
		else{
			printf("can't find epos!!!\n");
			exit(1);
		}
		key=tau;
		AddToMap(&b3d->ActionMap);
	}
}

void CAction::Kill(){
	if(currentmap==&b3d->ActionMap){
		CActionMap::iterator eepos,epos=GetPos(currentmap);
		currentmap->erase(epos);
		key=double(b3d->NACTIONSMAX+b3d->nactionkills);
		b3d->nactionkills+=1;
		AddToMap(b3d->DeadActionMap.end(),&b3d->DeadActionMap);
		CPart *part;
		CPartMap::iterator ppos=partmap.begin();
		while(ppos!=partmap.end()){
			part=ppos->second;
			epos=part->actionmap.begin();
			while(epos!=part->actionmap.end()){
				eepos=epos;
				++eepos;
				if(epos->second==this)
					part->actionmap.erase(epos);
				epos=eepos;
			}
			++ppos;
		}
	}
}

void CAction::AddToMap(CActionMap *newmap){
	currentmap=newmap;
	newmap->insert(CActionPair(key,this));
}

void CAction::AddToMap(CActionMap::iterator guess,CActionMap *newmap){
	currentmap=newmap;
	newmap->insert(guess,CActionPair(key,this));
}

void CAction::AddPart(CPart *part){
	partmap.insert(CPartPair(part->key,part));
}

void CAction::Print(){
	printf("___________ type=%d, tau=%g, nparts=%d ___________\n",type,tau,int(partmap.size()));
	CPartMap::iterator ppos;
	CPart *part;
	for(ppos=partmap.begin();ppos!=partmap.end();++ppos){
		part=ppos->second;
		part->Print();
	}
	printf("_________________________________________________\n");
}

void CAction::CheckPartList(){
	CPart *part;
	CPartMap::iterator ppos,ppos2;
	CPartMap *pmap;
	ppos=partmap.begin();
	while(ppos!=partmap.end()){
		part=ppos->second;
		ppos2=part->GetPos(&(b3d->PartMap));
		if(ppos2==b3d->PartMap.end()){
			printf("____________ CAction::CheckPartList FATAL, action type=%d ________________\n",type);
			printf("iterator not in expected pmap\n");
			part->Print();
			//exit(1);
		}
		++ppos;
	}
}

void CAction::PerformDensCalc(){
	int itau;
	itau=lrint(floor((tau-0.001)/b3d->DENSWRITE_DELTAU));
	if(itau>=b3d->DENSWRITE_NTAU){
		printf("trying to perform CAction::DensCalc() for itau>=DENSWRITE_NTAU, =%d",itau);
		exit(1);
	}
	CB3DCell *cell;
	int ix,iy,ieta;
	for(ix=0;ix<2*b3d->NXY;ix++){
		for(iy=0;iy<2*b3d->NXY;iy++){
			for(ieta=0;ieta<2*b3d->NETA;ieta++){
				cell=b3d->cell[ix][iy][ieta];
				cell->dens[itau]+=cell->partmap.size();
			}
		}
	}
}

#endif
