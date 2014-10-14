#ifndef __PART_CC__
#define __PART_CC__

#include "b3d.h"
CB3D *CPart::b3d=NULL;

CPart::CPart(){
	cell=NULL;
	currentmap=NULL;
	weight=1;
	tau0=0.0;
	reality=true;
	currentmap=NULL;
	cell=NULL;
	active=false;
	balanceID=0;
}

CPart::CPart(int keyset){
	weight=1;
	reality=true;
	key=keyset;
	tau0=0.0;
	tauexit=0.0;
	tau_lastint=0.0;
	taudecay=0.0;
	currentmap=&b3d->DeadPartMap;
	cell=NULL;
	actionmap.clear();
	active=false;
	resinfo=b3d->reslist->GetResInfoPtr(211);
	for(int alpha=0;alpha<4;alpha++){
		r[alpha]=p[alpha]=0.0;
	}
	y=eta=0.0;
	nscatt=0;
	reality=true;
	weight=1;
	balanceID=0;
	b3d->DeadPartMap.insert(CPartPair(key,this));
	b3d->npartstot+=1;
}

CPart::~CPart(){
}

void CPart::Copy(CPart *part){  //copies all info except actionmap
	int alpha;
	tau0=part->tau0;
	tau_lastint=part->tau_lastint;
	tauexit=part->tauexit;
	taudecay=part->taudecay;
	y=part->y;
	eta=part->eta;
	for(alpha=0;alpha<4;alpha++){
		p[alpha]=part->p[alpha];
		r[alpha]=part->r[alpha];
	}
	msquared=part->msquared;
	resinfo=part->resinfo;
	//reality=part->reality;
	//weight=part->weight;
}

void CPart::CopyPositionInfo(CPart *part){  //copies all info except actionmap
	int alpha;
	tau0=part->tau0;
	eta=part->eta;
	for(alpha=0;alpha<4;alpha++){
		r[alpha]=part->r[alpha];
	}
}

void CPart::CopyMomentumInfo(CPart *part){  //copies all info except actionmap
	int alpha;
	y=part->y;
	msquared=part->msquared;
	for(alpha=0;alpha<4;alpha++){
		p[alpha]=part->p[alpha];
	}
}

void CPart::InitBalance(int IDset,double rxset,double ryset,double tauset,double etaset,double pxset,double pyset,double mset,double rapidityset,int weightset,bool realset,int balanceIDset){
	realset=false;
	Init(IDset,rxset,ryset,tauset,etaset,pxset,pyset,mset,rapidityset,weightset,realset);
	balanceID=balanceIDset;
}

void CPart::Init(int IDset,double rxset,double ryset,double tauset,double etaset,double pxset,double pyset,double mset,double rapidityset,int weightset,bool realset){
	double et;
	CResInfo *resinfoptr;
	int ID;
	balanceID=0;
	resinfo=b3d->reslist->GetResInfoPtr(IDset);
	ID=resinfo->code;
	if(ID!=IDset){
		printf("ID mismatch, ID=%d, resinfo->codeID=%d\n",IDset,ID);
	}
	p[1]=pxset; p[2]=pyset; msquared=mset*mset; y=rapidityset;
	r[1]=rxset; r[2]=ryset; tau0=tauset; eta=etaset;
	tau_lastint=tau0-1.0E-6;
	nscatt=0;
	weight=weightset;
	reality=realset;
	if(fabs(eta)>b3d->ETAMAX){
		printf("in part->init, eta out of bounds, =%g\n",eta);
	}
	resinfoptr=b3d->reslist->GetResInfoPtr(ID);
	if(resinfoptr->decay==false){
		msquared=resinfoptr->mass*resinfoptr->mass;
	}
	r[3]=tau0*sinh(eta);
	r[0]=tau0*cosh(eta);
	et=sqrt(p[1]*p[1]+p[2]*p[2]+msquared);
	p[3]=et*sinh(y);
	Setp0();
	if(tau0<0.0){
		printf("FATAL: tau0<0, tau0^2=%g\n",tau0);
		Print();
		exit(1);
	}
	if(b3d->BJORKEN && fabs(eta)>b3d->ETAMAX){
		CyclicReset();
		printf("performed cyclic reset in CPart::Init()\n");
	}
	active=false;
	if(b3d->COLLISIONS==true)
		ChangeMap(&(b3d->PartMap));
	else
		ChangeMap(&(b3d->FinalPartMap));
	b3d->AddAction_Activate(this);
	actionmother=b3d->nactions;
}

void CPart::CyclicReset(){
	double eta_offset,etamax=b3d->ETAMAX;
	double mt;
	while(fabs(eta)>etamax){
		eta_offset=2.0*etamax*floor(0.5+0.5*eta/etamax);
		eta-=eta_offset;
		y-=eta_offset;
		mt=sqrt(p[0]*p[0]-p[3]*p[3]);
		p[3]=mt*sinh(double(y));
		p[0]=sqrt(mt*mt+p[3]*p[3]);
		r[3]=tau0*sinh(double(eta));
		r[0]=tau0*cosh(double(eta));
	}
}

void CPart::Print(){
	printf("________________ PART INFO FOR key=%d _____________________________\n",key);
	printf("Minv^2=%g, ID=%d\n",p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3],resinfo->code);
	printf("ID=%d, m_onshell=%g, M=%g, tau0=%g=?%g, tauexit=%g\n r=(%g,%g,%g,%g) eta=%g=?%g\n", 
		resinfo->code,resinfo->mass,sqrt(msquared),double(tau0),sqrt(r[0]*r[0]-r[3]*r[3]),tauexit,r[0],r[1],r[2],r[3],eta,GetEta(tau0));
	printf("weight=%d, reality=%d, key=%d, actionmother=%d, active=%d, balanceID=%d\n",
		weight,int(reality),key,actionmother,int(active),balanceID);
	printf("p=(%15.9e,%15.9e,%15.9e,%15.9e), y=%g =? %g\n",p[0],p[1],p[2],p[3],y,atanh(p[3]/p[0]));
	string currentmapname="IN CELL";
	if(currentmap==&(b3d->PartMap)) currentmapname="PartMap";
	if(currentmap==&(b3d->DeadPartMap)) currentmapname="DeadPartMap";
	if(currentmap==&(b3d->FinalPartMap)) currentmapname="FinalPartMap";
	printf("currentmap=%s\n",currentmapname.c_str());
	if(cell==NULL) printf("CELL=NULL\n");
	if(nextcell==NULL) printf("NEXTCELL=NULL\n");
	if(cell!=NULL) printf("Cell No: ix=%d, iy=%d, ieta=%d\n",cell->ix,cell->iy,cell->ieta);
	printf("________________________________________________________________________\n");
}

void CPart::ChangeMap(CPartMap *newmap){
	DeleteFromCurrentMap();
	AddToMap(newmap);
}

CPartMap::iterator CPart::DeleteFromCurrentMap(){
	CPartMap::iterator neighbor;
	CPartMap::iterator ppos=GetPos(currentmap);
	neighbor=ppos;
	neighbor++;
	if(ppos==currentmap->end()){
		printf("FATAL: In CPart::DeleteFromCurrentMap, can't find ppos!!!\n");
		Print();
		printf("currentmap has length %d\n",int(currentmap->size()));
		exit(1);
	}
	else currentmap->erase(ppos);
	currentmap=NULL;
	return neighbor;
}

CPartMap::iterator CPart::DeleteFromMap(CPartMap *partmap){
	CPartMap::iterator neighbor;
	CPartMap::iterator ppos=GetPos(partmap);
	neighbor=ppos;
	if(ppos==partmap->end()){
		printf("FATAL: In CPart::DeleteFromMap, can't find ppos!!!\n");
		Print();
		exit(1);
	}
	else{
		neighbor++;
		partmap->erase(ppos);
	}
	partmap=NULL;
	return neighbor;
}

void CPart::AddToMap(CPartMap *newmap){
	newmap->insert(CPartPair(key,this));
	if(newmap==&b3d->PartMap || newmap==&b3d->DeadPartMap || newmap==&b3d->FinalPartMap)
		currentmap=newmap;
}

void CPart::AddToMap(CPartMap::iterator guess,CPartMap *newmap){
	newmap->insert(guess,CPartPair(key,this));
	currentmap=newmap;
}

void CPart::SubtractAction(CAction *action){
	CActionMap::iterator epos=action->GetPos(&actionmap);
	if(epos!=actionmap.end()) actionmap.erase(epos);
}

void CPart::AddAction(CAction *action){
	actionmap.insert(CActionPair(action->key,action));
}

void CPart::Propagate(double tau){
	if(b3d->BJORKEN && abs(eta)>b3d->ETAMAX){
		printf("eta screwy before propagation\n");
		printf("eta=%g\n",eta);
	}
	double t0,vz,gamma,gammav;
	CPartMap *pmap=currentmap;
	CPartMap::iterator neighbor;
	if(active==true){
		eta=GetEta(tau);//y-asinh((tau0/tau)*sinh(y-eta));
		tau0=tau;
		gamma=cosh(eta);
		gammav=sinh(eta);
		t0=r[0];
		r[0]=tau0*gamma;
		r[3]=tau0*gammav;
		r[1]+=(p[1]/p[0])*(r[0]-t0);
		r[2]+=(p[2]/p[0])*(r[0]-t0);
	}
	else{
		r[0]=tau*cosh(eta);
		r[3]=tau*sinh(eta);
		tau0=tau;
	}
	if(currentmap==&(b3d->PartMap) && tau<b3d->TAUCOLLMAX && fabs(eta)>0.000001+b3d->ETAMAX){
		printf("eta out of bounds after propagation\n");
		Print();
		exit(1);
	}
}

/*
CPartMap::iterator CPart::GetPos(CPartMap *pmap){
	pair<CPartMap::iterator,CPartMap::iterator> ppospair;
	CPartMap::iterator ppos,pend;
	long long int count;
	count=pmap->count(key);
	ppospair=pmap->equal_range(key);
	ppos=ppospair.first;
	pend=ppospair.second;
	while(ppos->second!=this && ppos!=pend){
		++ppos;
	}
	if(ppos->second!=this){
		//printf("failed in CPart::GetPos()\n");
		//printf("key=%d, count=%d\n",key,count);
		return pmap->end();
	}
	else return ppos;
}
*/

CPartMap::iterator CPart::GetPos(CPartMap *pmap){
	CPartMap::iterator ppos=pmap->find(key);
	return ppos;
}

void CPart::CheckMap(CPartMap *expectedpartmap){
	if(currentmap!=expectedpartmap){
		printf("FATAL: XXXXXXXXX particle not in expected map XXXXXXXXX\n");
		if(currentmap==&(b3d->DeadPartMap)){
			printf("particle in DeadPartMap\n");
		}
		if(currentmap==&(b3d->PartMap)){
			printf("particlein PartMap\n");
		}
		Print();
		exit(1);
	}
}

double CPart::GetMass(){
	if(resinfo->code==22)
		return 0.0;
	else
		return sqrt(msquared);
}

void CPart::Setp0(){
	p[0]=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]+msquared);
}

void CPart::SetY(){
	y=asinh(p[3]/GetMT());
}

void CPart::SetMass(){
	msquared=p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3];
}

double CPart::GetEta(double tau){
	double dy,deta,dtau0,dtau;
	if(active){
		dy=y;
		deta=eta;
		dtau0=tau0;
		dtau=tau;
		deta=dy-asinh((dtau0/dtau)*sinh(dy-deta));
		return deta;
	}
	else return eta;
}

double CPart::GetMT(){
	if(p[0]<fabs(p[3])){
		printf("CPart::GetMT, catastrophe\n");
		Print();
		exit(1);
	}
	return sqrt(p[0]*p[0]-p[3]*p[3]);
}

void CPart::KillActions(){
	CActionMap::iterator ep,epp;
	CAction *action;
	ep=actionmap.begin();
	while(ep!=actionmap.end()){
		action=ep->second;
		epp=ep; ++epp;
		action->Kill();
		ep=epp;
	}
	actionmap.clear();
}

void CPart::Kill(){
	KillActions();
	if(cell!=NULL){
		RemoveFromCell();
		cell=NULL;
	}
	DeleteFromCurrentMap();
	eta=0.0;
	nscatt=0;
	tau0=tau_lastint=tauexit=-1.0;
	AddToMap(b3d->DeadPartMap.begin(),&(b3d->DeadPartMap));
	weight=1;
	active=false;
}

void CPart::BjorkenTranslate(){
	if(eta<-b3d->ETAMAX || eta>b3d->ETAMAX){
		printf("eta out of bounds before translation\n");
		Print();
		exit(1);
	}
	double mt;
	if(cell->ieta==0){
		eta+=2.0*b3d->ETAMAX;
		y+=2.0*b3d->ETAMAX;
	}
	else{
		eta-=2.0*b3d->ETAMAX;
		y-=2.0*b3d->ETAMAX;
	}
	mt=sqrt(p[0]*p[0]-p[3]*p[3]);
	p[3]=mt*sinh(y);
	p[0]=mt*cosh(y);
	r[0]=tau0*cosh(eta);
	r[3]=tau0*sinh(eta);
	SetMass();
}

void CPart::BjorkenUnTranslate(){
	double mt;
	if(eta>b3d->ETAMAX){
		eta-=2.0*b3d->ETAMAX;
		y-=2.0*b3d->ETAMAX;
	}
	else{
		eta+=2.0*b3d->ETAMAX;
		y+=2.0*b3d->ETAMAX;
	}
	mt=sqrt(p[0]*p[0]-p[3]*p[3]);
	p[3]=mt*sinh(y);
	p[0]=mt*cosh(y);
	r[0]=tau0*cosh(eta);
	r[3]=tau0*sinh(eta);
	SetMass();
}

void CPart::FindCollisions(){
	int ix,iy,ieta;
	double taucoll;
	CPart *part2,*part1=this;
	CPartMap::iterator ppos;
	CB3DCell *cell2;
	for(ix=0;ix<3;ix++){
		for(iy=0;iy<3;iy++){
			for(ieta=0;ieta<3;ieta++){
				cell2=cell->neighbor[ix][iy][ieta];
				if(cell2!=NULL){
					ppos=cell2->partmap.begin();
					while(ppos!=cell2->partmap.end()){
						part2=ppos->second;
						if(part1!=part2 && part1->actionmother!=part2->actionmother){
							b3d->FindCollision(part1,part2,taucoll);
						}
						++ppos;
					}
				}
			}
		}
	}
}

CB3DCell *CPart::FindCell(){
	if(tau0>b3d->TAUCOLLMAX || !b3d->COLLISIONS){
		return NULL;
	}
	int ieta,ix,iy;
	double deta=b3d->ETAMAX/double(b3d->NETA);
	ieta=lrint(floor((eta+b3d->ETAMAX)/deta));
	if(ieta<0 ||ieta>=2*b3d->NETA){
		return NULL;
	}
	double dx=b3d->XYMAX/double(b3d->NXY);
	double dy=dx;
	ix=lrint(floor((r[1]+b3d->XYMAX)/dx));
	if(ix<0 || ix>=2*b3d->NXY){
		return NULL;
	}
	iy=lrint(floor((r[2]+b3d->XYMAX)/dy));
	if(iy<0 || iy>=2*b3d->NXY){
		return NULL;
	}
	return b3d->cell[ix][iy][ieta];
}

void CPart::FindDecay(){
	CAction *action;
	double t,gamma,vz,newt,newz;
	t=HBARC/resinfo->width;
	gamma=p[0]/sqrt(msquared);
	t=-t*gamma*log(b3d->randy->ran());
	vz=p[3]/p[0];
	newt=r[0]+t;
	newz=r[3]+vz*t;
	taudecay=sqrt(newt*newt-newz*newz);
	if(taudecay<tauexit || tauexit>b3d->TAUCOLLMAX || cell==NULL){
		b3d->AddAction_Decay(this,taudecay);
	}
}

void CPart::FindCellExit(){
	if(active){
		double t,taux,tauy,taueta,z;
		double etamax=cell->etamax,etamin=cell->etamin;
		nextcell=NULL;
		tauexit=1.0E50;
		taueta=taux=tauy=tauexit;
		double vx=p[1]/p[0];
		double vy=p[2]/p[0];
		double vz=p[3]/p[0];

		if(vx>0)
			t=(cell->xmax-r[1])/vx;
		else
			t=(cell->xmin-r[1])/vx;
		t=t+r[0];
		z=r[3]+vz*(t-r[0]);
		taux=sqrt(t*t-z*z);

		if(vy>0)
			t=(cell->ymax-r[2])/vy;
		else
			t=(cell->ymin-r[2])/vy;
		t=t+r[0];
		z=r[3]+vz*(t-r[0]);
		tauy=sqrt(t*t-z*z);

		if(y<cell->etamin)
			taueta=tau0*sinh(y-eta)/sinh(y-etamin);
		else if(y>cell->etamax)
			taueta=tau0*sinh(y-eta)/sinh(y-etamax);

		if(taux<tauy && taux<taueta){
			tauexit=taux;
			if(vx<0)
				nextcell=cell->neighbor[0][1][1];
			else
				nextcell=cell->neighbor[2][1][1];
		}
		else if (tauy<taueta){
			tauexit=tauy;
			if(vy<0)
				nextcell=cell->neighbor[1][0][1];
			else
				nextcell=cell->neighbor[1][2][1];
		}
		else{
			tauexit=taueta;
			if(y<etamin)
				nextcell=cell->neighbor[1][1][0];
			else
				nextcell=cell->neighbor[1][1][2];
		}

		if(tauexit<b3d->TAUCOLLMAX)
			b3d->AddAction_ExitCell(this);
	}
}

void CPart::FindActions(){
	if(active!=true){
		printf("CPart::FindActions(), trying to Reset Inactive particle\n");
		exit(1);
	}
	KillActions();

	if(msquared<resinfo->minmass*resinfo->minmass){
		msquared=resinfo->minmass*resinfo->minmass;
		Setp0();
	}

	if(cell!=NULL){
		if(b3d->COLLISIONS && tau0<b3d->TAUCOLLMAX && (nscatt<b3d->NSCATT_MAX || reality)){
			FindCellExit();
			FindCollisions();
		}
		else{
			ChangeCell(NULL);
			ChangeMap(&b3d->FinalPartMap);
		}
	}
	if(resinfo->decay)
		FindDecay();
}

double CPart::GetPseudoRapidity(){
	double pmag,eta_ps;
	pmag=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
	eta_ps=atanh(p[3]/pmag);
	return eta_ps;
}

void CPart::Boost(FourVector &u){
	BoostP(u);
	BoostR(u);
}

void CPart::BoostP(FourVector &u){
	int alpha;
	FourVector pprime;
	Misc::Boost(u,p,pprime);
	for(alpha=0;alpha<4;alpha++)
		p[alpha]=pprime[alpha];
	y=atanh(p[3]/p[0]);
}

void CPart::BoostR(FourVector &u){
	int alpha;
	FourVector rprime;
	Misc::Boost(u,r,rprime);
	for(alpha=0;alpha<4;alpha++)
		r[alpha]=rprime[alpha];
	eta=atanh(r[3]/r[0]);
	tau0=sqrt(r[0]*r[0]-r[3]*r[3]);
}

void CPart::RemoveFromCell(){
	if(cell!=NULL){
		CPartMap::iterator ppos=GetPos(&(cell->partmap));
		if(ppos==cell->partmap.end()){
			printf("FATAL: In CPart::RemoveFromCell, can't find ppos!!!\n");
			Print();
			printf("cell partmap has length %d\n",int(cell->partmap.size()));
			exit(1);
		}
		else{
			cell->partmap.erase(ppos);
		}
	}
}

void CPart::CheckCell(){
	if(cell!=NULL){
		CPartMap::iterator ppos=GetPos(&(cell->partmap));
		if(ppos==cell->partmap.end()){
			printf("FATAL: In CPart::RemoveFromCell, can't find ppos!!!\n");
			Print();
			printf("cell partmap has length %d\n",int(cell->partmap.size()));
			exit(1);
		}
	}
}

void CPart::ChangeCell(CB3DCell *newcell){
	if(newcell!=cell){
		if(cell!=NULL)
			RemoveFromCell();
		if(newcell!=NULL){
			newcell->partmap.insert(CPartPair(key,this));
		}
		cell=newcell;
	}
}

void CPart::GetHBTPars(double &t,double &rout,double &rside,double &rlong){
	const double tcompare=15.0;
	double pt,ptsquared,et;
	rlong=tau0*sinh(eta-y);
	t=tau0*cosh(eta-y);
	ptsquared=p[1]*p[1]+p[2]*p[2];
	pt=sqrt(ptsquared);
	et=sqrt(ptsquared+msquared);
	rout=(p[1]*r[1]+p[2]*r[2])/pt;
	rout=rout-(pt/et)*(t-tcompare);
	rside=(p[1]*r[2]-p[2]*r[1])/pt;
}

//void CB3DBinaryPartInfo::Print(){
//	printf("ID=%d, tau=%g, x=%g, y=%g, eta=%g\n",ID,tau,x,y,eta);
//	printf("px=%g,py=%g, rapidity=%g, weight=%g, reality=%d\n",px,py,rapidity,weight,int(reality));
//}
#endif
