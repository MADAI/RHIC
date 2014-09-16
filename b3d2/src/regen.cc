#include "regen.h"
CB3D *CRegenerate::b3d=NULL;

CRegenerate::CRegenerate(){
	CResInfoMap::iterator rpos;
	CResInfoMap *resmap=&(b3d->reslist->resmap);
	CResInfo *resinfo;
	randy=b3d->randy;
	int iT;
	for(rpos=resmap->begin();rpos!=resmap->end();++rpos){
		resinfo=rpos->second;
		if(resinfo->baryon==1){
			baryonmap.insert(CResInfoPair(resinfo->code,resinfo));
		}
	}
	Tmin=50.0;
	delT=2.0;
	sigmavmax=5.8/b3d->NSAMPLE;
	NT=lrint((180-Tmin)/delT);
	BDens.resize(NT+1);
	BDensMap.resize(NT+1);
	for(iT=0;iT<=NT;iT++){
		CalculateBDens(iT);
	}
}

void CRegenerate::CalculateBDens(int iT){
	double epsilon,P,dens0,sigma2,dedt;
	double T=Tmin+delT*iT;
	int S;
	CResInfo *resinfo;
	CResInfoMap *resmap=&(b3d->reslist->resmap);
	CDensMap::iterator dpos;
	CResInfoMap::iterator rpos;
	for(S=0;S<=3;S++){
		BDens[iT][S]=0.0;
	}
	for(rpos=baryonmap.begin();rpos!=baryonmap.end();++rpos){
		resinfo=rpos->second;
		if(resinfo->baryon==1){
			S=abs(resinfo->strange);
			b3d->reslist->freegascalc_onespecies(resinfo->mass,T,epsilon,P,dens0,sigma2,dedt);
			BDens[iT][S]+=dens0*(2.0*resinfo->spin+1.0);
			BDensMap[iT][S].insert(CDensPair(BDens[iT][S],resinfo));
		}
	}
}

bool CRegenerate::CheckForRegeneration(CB3DCell *cell,CMuTInfo *muTinfo){
	CPart *part;
	CResInfo *resinfo1,*resinfo2;
	bool regencheck=false;
	int nproducts,iproduct;
	/*
	double Tpi,mupi,TK,muK;
	FourVector upi,uK;
	FourVector p1,p2;
	muTinfo->MuTCalc();
	Tpi=muTinfo->Tpi;
	mupi=muTinfo->mupi;
	TK=muTinfo->TK;
	muK=muTinfo->muK;
	int itau=lrint(b3d->tau/CMuTInfo::DELTAU)-1;
	*/
	GetdNMax(muTinfo);
	MC_Nbar+=dNMaxNet;
	while(MC_Nbar>MC_Ntarget){
		MC_Ntarget-=log(b3d->randy->ran());
		if(GetBBbarResInfo(muTinfo,resinfo1,resinfo2)){
			b3d->GetDeadParts(part1,part2);
			part1->resinfo=resinfo1;
			part2->resinfo=resinfo2;
			GetP1P2SigmaV(muTinfo);
			if(sigmavrel>sigmavmax){
				printf("sigmavrel=%g, is > than sigmavmax=%g\n",sigmavrel,sigmavmax);
			}
			if(b3d->randy->ran()<sigmavrel/sigmavmax){
				regencheck=true;
				b3d->nregenerate+=1;
				b3d->GetDeadParts(product);
				FillOutPartPRInfo(part1,cell,muTinfo);
				FillOutPartPRInfo(part2,cell,muTinfo);
				b3d->Annihilate(part1,part2,nproducts,product);
				
				part1->reality=part2->reality=true;
				part1->active=part2->active=true;
				part1->weight=part2->weight=1;
				part1->tau_lastint=part2->tau_lastint=b3d->tau;
				part1->actionmother=part2->actionmother=b3d->nactionstot;
				part1->FindActions();
				part2->FindActions();
				
				for(iproduct=0;iproduct<nproducts;iproduct++){
					part=product[iproduct];
					FillOutPartPRInfo(part,cell,muTinfo);
					part->reality=false;
					part->nscatt=1;
					part->active=true;
					part->weight=-1;
					part->tau_lastint=b3d->tau;
					part->actionmother=b3d->nactionstot;
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
		}
	}
	return regencheck;
}

void CRegenerate::GetdNMax(CMuTInfo *muTinfo){
	double xpi,xK,dNNet,T,Tpi,TK,mupi,muK,tau,volume,pidens,Kdens;
	int smax,itau,norm;
	tau=b3d->tau;
	volume=tau*2.0*b3d->ETAMAX*b3d->DXY*b3d->DXY;
	norm=CMuTInfo::NETEVENTS*b3d->NSAMPLE;
	volume=volume;
	itau=lrint(tau/CMuTInfo::DELTAU)-1;
	pidens=muTinfo->Npi[itau]/(volume*norm);
	Kdens=muTinfo->NK[itau]/(volume*norm);
	if(pidens>0.025 && muTinfo->Npi[itau]>5){
		muTinfo->MuTCalc();
		int S1,S2,iT,smax;
		Tpi=muTinfo->Tpi;
		TK=muTinfo->TK;
		mupi=muTinfo->mupi;
		muK=muTinfo->muK;
		xpi=exp(muTinfo->mupi/muTinfo->Tpi);
		xK=exp(muTinfo->muK/muTinfo->TK);
		dNMaxNet=0.0;
		smax=3;
		if(muTinfo->NK[itau]<=5 || TK<Tmin)
			smax=0;
		if(Tpi>Tmin){
			for(S1=0;S1<=smax;S1++){
				for(S2=0;S2<=S1;S2++){
					NK=abs(S1-S2);
					T=0.2*(Tpi*(5-NK)+TK*NK);
					iT=lrint((T-Tmin)/delT);
					if(iT>NT)
						iT=NT;
					if(iT>0 && iT<=NT){
						dNMax[S1][S2]=pow(xpi,5-NK)*pow(xK,NK)*BDens[iT][S1]*BDens[iT][S2]*volume*b3d->NSAMPLE;
						dNMaxNet+=dNMax[S1][S2];
						if(S1!=S2){
							dNMax[S2][S1]=dNMax[S1][S2];
							dNMaxNet+=dNMax[S2][S1];
						}
					}
					else{
						dNMax[S1][S2]=0.0;
						if(S1!=S2)
							dNMax[S2][S1]=0.0;
					}
				}
			}
			//if(dNMaxNet>0.01)
				//printf("dNMaxNet=%g, pidens=%g, xpi=%g, Tpi=%g, Kdens=%g, xK=%g, TK=%g, muK=%g, itau=%d\n",
			//dNMaxNet,pidens,xpi,Tpi,Kdens,xK,TK,muK,itau);
		}
	}
	else{
		dNMaxNet=0;
	}
	
}

bool CRegenerate::GetBBbarResInfo(CMuTInfo *muTinfo,CResInfo *&resinfo1,CResInfo *&resinfo2){
	CDensMap::iterator dpos;
	CResInfo *resinfo;
	double key,wtot=0.0,wsum,w[16];
	int S1,S2,iS,NK,iT,pid,itau=lrint(b3d->tau/CMuTInfo::DELTAU)-1;
	double Tpi=muTinfo->Tpi,mupi=muTinfo->mupi,TK=muTinfo->TK,muK=muTinfo->muK;
	double T;
	if(muTinfo->NK[itau]>1){
		for(iS=0;iS<16;iS++){
			w[iS]=0.0;
			S1=(iS/4);
			S2=iS-4*S1;
			NK=abs(S1-S2);
			T=0.2*(TK*NK+Tpi*(5-NK));
			iT=lrint((T-Tmin)/delT);
			if(iT>=0 && iT<=NT){
				w[iS]=exp((muK/TK)*NK+(mupi/Tpi)*(5-NK));
				w[iS]*=BDens[iT][S1]*BDens[iT][S2];
				wtot+=w[iS];
			}
		}
		iS=0;
		wsum=0.0;
		while(b3d->randy->ran()>wsum){
			wsum+=w[iS]/wtot;
			iS+=1;
		}
	}
	else{
		T=Tpi;
		iT=lrint((T-Tmin)/delT);
		wtot=1.0;
		iS=0;
	}
	if(iT>0 && iT<=NT && wtot>1.0E-13){
		S1=(iS/4);
		S2=iS-4*S1;
		key=BDens[iT][S1]*(b3d->randy->ran());
		dpos=(BDensMap[iT][S1]).lower_bound(key);
		resinfo1=dpos->second;
		key=BDens[iT][S2]*(b3d->randy->ran());
		dpos=(BDensMap[iT][S2]).lower_bound(key);
		pid=(dpos->second)->code;
		resinfo2=b3d->reslist->GetResInfoPtr(-pid);
		return true;
	}
	else
		return false;
}

//cell info only used for choosing x,y coordinates
void CRegenerate::FillOutPartPRInfo(CPart *part,CB3DCell *cell,CMuTInfo *muTinfo){
	FourVector u;
	part->tau0=b3d->tau;
	part->r[1]=cell->xmin+b3d->randy->ran()*(cell->xmax-cell->xmin);
	part->r[2]=cell->ymin+b3d->randy->ran()*(cell->ymax-cell->ymin);
	part->eta=-b3d->ETAMAX+2.0*b3d->ETAMAX*b3d->randy->ran();
	u[1]=0.2*(muTinfo->uxpi*(5-NK)+muTinfo->uxK*NK);
	u[2]=0.2*(muTinfo->uypi*(5-NK)+muTinfo->uyK*NK);
	u[3]=sinh(part->eta);
	u[0]=sqrt(1.0+u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
	part->SetMass();
	part->BoostP(u);
	part->SetY();
	part->r[3]=part1->tau0*sinh(part1->eta);
	part->ChangeCell(part->FindCell());
}

void CRegenerate::GetP1P2SigmaV(CMuTInfo *muTinfo){
	CRandom *randy=b3d->randy;
	double T,vrel;
	NK=abs(part1->resinfo->strange-part2->resinfo->strange);
	T=0.2*(muTinfo->Tpi*(5-NK)+muTinfo->TK*NK);
	randy->generate_boltzmann(part1->resinfo->mass,T,part1->p);
	randy->generate_boltzmann(part2->resinfo->mass,T,part2->p);
	double sigma=b3d->GetAnnihilationSigma(part1,part2,vrel);
	sigmavrel=sigma*vrel;
	if(sigmavrel>sigmavmax){
		printf("CRegenerate::GetP1P2SigmaV, sigmavrel=%g, sigmavmax=%g\n",sigmavrel,sigmavmax);
		exit(1);
	}
}