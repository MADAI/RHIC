#ifndef __ACTION_PERFORM_COLLIDE_CC__
#define __ACTION_PERFORM_COLLIDE_CC__
#include "b3d.h"

void CAction::PerformCollide(){
	double scompare,Mtest;
	double tautest;
	int ncolls,colltype;
	CPart *part1,*part2,*part;
	CPartMap::iterator ppos1,ppos2,pbegin,pend,ppos;
	double sigma,einitial,efinal;
	
	if(partmap.size()!=2){
		printf("FATAL: wrong number of particles in partmap for collision, =%d\n",int(partmap.size()));
		exit(1);
	}
	ppos1=partmap.begin();
	ppos2=ppos1; ++ppos2;
	part1=ppos1->second;
	part2=ppos2->second;
	
	part1->Propagate(b3d->tau);
	part2->Propagate(b3d->tau);
	
	if(part1->active && part2->active){
		colltype=b3d->Collide(part1,part2);
  }
	else{
		colltype=-1;
		printf("Trying to collide with dead dude\n");
		exit(1);
	}

	if(colltype==0)
		b3d->npass++;
	if(colltype==1)
		b3d->nmerge+=1;
	if(colltype==2)
		b3d->nscatter+=1;
	if(colltype == 3)
		b3d->ninelastic+=1;
	if(colltype==4)
		b3d->nannihilate+=1;
	b3d->nactions+=1;
}

#endif
