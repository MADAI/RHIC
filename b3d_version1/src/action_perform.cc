#ifndef __ACTION_PERFORM_CC__
#define __ACTION_PERFORM_CC__
#include "b3d.h"

void CAction::Perform(){
	CPartMap::iterator ppos;
	int ttype=type;
	//printf("performing action type %d, tau=%g, action listid=%d, nactions=%d\n",type,tau,listid,int(b3d->nactions));

/*
	ppos=partmap.begin();
	while(ppos!=partmap.end()){
		if(ppos->second->listid==302){
			printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
			printf("beginning action action type=%d, tau=%g, action listid=%d\n",type,tau,listid);
			ppos->second->Print();
			if(ppos->second->cell!=NULL)
				ppos->second->cell->Print();
			//Misc::Pause();
		}
		++ppos;
	}
	*/

	if(currentmap!=&(b3d->ActionMap)){
		printf("FATAL: trying to perform dead action\n");
		exit(1);
	}
	Kill();

	b3d->tau=tau;

	if(tau+1.0E-4<b3d->tau){
		printf("FATAL:: action earlier than tau!!!!, b3d->tau=%15.10e, action tau=%15.10e\n",b3d->tau,tau);
		exit(1);
	}
	//printf("performing action: type=%d, tau=%g, listid=%d\n",type,tau,listid);
	if(ttype==0) PerformActivate();
	else if(ttype==1) PerformDecay();
	else if(ttype==2) PerformCollide();
#ifdef VIZWRITE
	else if(ttype==3) PerformVizWrite();
#endif
	else if(ttype==4) PerformDensCalc();
	else if(ttype==6) PerformExitCell();
	else{
		printf("FATAL: action type = %d is unknown, exiting\n",type);
		exit(1);
	}
	/*
	ppos=partmap.begin();
	while(ppos!=partmap.end()){
		if(ppos->second->listid==302){
			printf("finished action action type=%d, tau=%g, action listid=%d\n",type,tau,listid);
			ppos->second->Print();
			if(ppos->second->cell!=NULL)
				ppos->second->cell->Print();
			//Misc::Pause();
			printf("YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY\n");
		}
		++ppos;
	}
	*/
	//printf("action finished, listid=%d, %d more to go\n",listid,int(b3d->ActionMap.size()));
}


#endif
