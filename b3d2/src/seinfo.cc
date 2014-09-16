#include "seinfo.h"

CSEInfo::CSEInfo(CB3D *b3dset){
	int itau;
	b3d=b3dset;
	DELTAU=1.0;
	NTAU=lrint(parameter::getD(b3d->parmap,"B3D_TAUCOLLMAX",20));
	TAU0=parameter::getD(b3d->parmap,"SEINFO_TAU0",5);
	ETAOVERS=0.0;
	NETEVENTS=0;
	R=10.0;
	epsilon.resize(NTAU);
	Tzz.resize(NTAU);
	Pbar.resize(NTAU);
	nhadrons.resize(NTAU);
	K0.resize(NTAU);
	F0.resize(NTAU);
	Zero();
}

void CSEInfo::Zero(){
	int itau;
	NETEVENTS=0;
	for(itau=0;itau<NTAU;itau++){
		Tzz[itau]=Pbar[itau]=epsilon[itau]=nhadrons[itau]=K0[itau]=F0[itau]=0.0;
	}
}

void CSEInfo::SECalc(){
	double e,pmag2,pz,et,r;
	CPartMap::iterator ppos;
	CPart *part;
	int itau;
	itau=lrint(b3d->tau/DELTAU)-1;
	for(ppos=b3d->PartMap.begin();ppos!=b3d->PartMap.end();ppos++){
		part=ppos->second;
		part->Propagate(b3d->tau);
		r=sqrt(part->r[1]*part->r[1]+part->r[2]*part->r[2]);
		if(r<R){
			et=sqrt(part->msquared+part->p[1]*part->p[1]+part->p[2]*part->p[2]);
			pz=et*sinh(part->y-part->eta);
			e=sqrt(et*et+pz*pz);
			pmag2=part->p[1]*part->p[1]+part->p[2]*part->p[2];
			pmag2+=pz*pz;
			Tzz[itau]+=pz*pz/e;
			Pbar[itau]+=pmag2/(3.0*e);
			epsilon[itau]+=e;
			K0[itau]+=pmag2*(1.0-0.2*pmag2/(e*e))/(3.0*e);
			F0[itau]+=pmag2*pmag2/(15.0*e*e);
			nhadrons[itau]+=1;
		}
	}
}

void CSEInfo::Print(){
	int itau,itau0=lrint(TAU0/DELTAU)-1;
	double etaf,eta;
	double Pf=b3d->sampler->Pf,nf=b3d->sampler->nhadronsf,ef=b3d->sampler->epsilonf,Tf=b3d->sampler->Tf;
	double sf=(Pf+ef)/Tf;
	double alpha,dpizzoveralpha_dt,p0,p1,p2;
	CResList *reslist=b3d->reslist;
	double tau,volume;
	char filename[150];
	sprintf(filename,"model_output/%s/%s/tij.dat",b3d->run_name.c_str(),b3d->qualifier.c_str());
	FILE *fptr=fopen(filename,"w");
	printf("#tau  epsilon   Pbar    Tzz   nhadrons  eta      K0      F0     -pizz/alpha\n");
	fprintf(fptr,"#tau  epsilon   Pbar    Tzz   nhadrons  eta      K0      F0     -pizz/alpha\n");
	for(itau=itau0;itau<NTAU;itau++){
		tau=(itau+1)*DELTAU;
		volume=2.0*b3d->ETAMAX*tau*PI*R*R*NETEVENTS;
		eta=0.75*tau*(Pbar[itau]-Tzz[itau])/volume;
		alpha=(TAU0/tau)*sqrt(F0[itau]/F0[itau0]);
		printf("%5.2f %7.3f %7.3f %7.3f %7.4f %7.4f %7.4f %7.2f %7.4f\n",
		tau,epsilon[itau]/volume,Pbar[itau]/volume,Tzz[itau]/volume,nhadrons[itau]/volume,eta,K0[itau]/volume,F0[itau]/volume,(Pbar[itau]-Tzz[itau])/(alpha*volume));
		fprintf(fptr,"%5.2f %7.3f %7.3f %7.3f %7.4f %7.4f %7.4f %7.2f %7.4f\n",
		tau,epsilon[itau]/volume,Pbar[itau]/volume,Tzz[itau]/volume,nhadrons[itau]/volume,eta,K0[itau]/volume,F0[itau]/volume,(Pbar[itau]-Tzz[itau])/(alpha*volume));
		if(itau==itau0)
			p0=(Pbar[itau]-Tzz[itau])/(alpha*volume);
		if(itau==itau0+1)
			p1=(Pbar[itau]-Tzz[itau])/(alpha*volume);
		if(itau==itau0+2)
			p2=(Pbar[itau]-Tzz[itau])/(alpha*volume);
		printf("tau=%g, p=%g, pNS=%g, alpha=%g\n",tau,(Pbar[itau]-Tzz[itau])/(alpha*volume),(4.0/3.0)*eta/(tau*alpha),alpha);
	}
	dpizzoveralpha_dt=(2.0*p1-1.5*p0-0.5*p2)/DELTAU;
	printf("&&&&&&&&&eta/s=%g D(pizz/alpha)/Dt=%g &&&&&&&&&&&&\n",ETAOVERS,dpizzoveralpha_dt);
	fprintf(fptr,"#&&&&&&&&&eta/s=%g D(pizz/alpha)/Dt=%g &&&&&&&&&&&&\n",ETAOVERS,dpizzoveralpha_dt);
	fclose(fptr);
}


