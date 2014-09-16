#ifndef __RESONANCES_CC__
#define __RESONANCES_CC__
#include "b3d.h"

CB3D *CResList::b3d=NULL;

CResList::CResList(){
	if(b3d!=NULL){
		parmap=&(b3d->parmap);
		ReadResInfo();
	}
}

CResList::~CResList(){
	CResInfo *resinfo;
	CResInfoMap::iterator rpos=resmap.begin();
	while(rpos!=resmap.end()){
		resinfo=rpos->second;
		resmap.erase(rpos);
		delete resinfo;
		rpos=resmap.begin();
	}
}

CResList::CResList(parameterMap* parmap_in){
	parmap=parmap_in;
	ReadResInfo();
}

CRandom *CResInfo::ranptr=new CRandom(-1234);

CResInfo::CResInfo(){
	minmass=0.0;
	branchlist.clear();
}

CMerge::CMerge(CResInfo *resinfo_in,double branching_in, int L_in){
	resinfo=resinfo_in;
	branching=branching_in;
	L = L_in;
	next=NULL;
}

void CResInfo::DecayGetResInfoPtr(int &nbodies,array<CResInfo *,5> &daughterresinfo){
	double randy,bsum;
	int ibody,ibranch;
	CBranchInfo *bptr;
	bptr=NULL;
	bsum=0.0;
	randy=ranptr->ran();
	ibranch=0;
	do{
		bptr=branchlist[ibranch];
		bsum+=bptr->branching;
		ibranch++;
		if(bsum>1.00000001){
			printf("FATAL: In DecayGetResInfo: bsum too large, = %g\n",bsum);
			exit(1);
		}
	}while(bsum<randy);
	nbodies=bptr->resinfoptr.size();
	for(ibody=0;ibody<nbodies;ibody++){
		daughterresinfo[ibody]=bptr->resinfoptr[ibody];
	}

}

bool CResInfo::CheckForDaughters(int codecheck){
	//checks to see if any decay daughters match code, for code=0, checking for charged parts
	int ibody,nbodies,ibranch;
	bool exists=false;
	CResInfo *daughter;
	CBranchInfo *bptr;
	CBranchList::iterator bpos;
	if(codecheck!=0){
		if(code==codecheck){
			exists=true;
			return exists;
		}
		if(decay){
			ibranch=0;
			do{
				bptr=branchlist[ibranch];
				ibranch++;
				nbodies=bptr->resinfoptr.size();
				for(ibody=0;ibody<nbodies;ibody++){
					daughter=bptr->resinfoptr[ibody];
					if(daughter->code==codecheck){
						exists=true;
						return exists;
					}
					if(daughter->CheckForDaughters(codecheck)){
						exists=true;
						return exists;
					}
				}
				bpos++;
			}while(bpos!=branchlist.end());
		}
	}
	else{
		if(charge!=0){
			exists=true;
			return exists;
		}
		if(decay){
			ibranch=0;
			do{
				bptr=branchlist[ibranch];
				ibranch++;
				nbodies=bptr->resinfoptr.size();
				for(ibody=0;ibody<nbodies;ibody++){
					daughter=bptr->resinfoptr[ibody];
					if(daughter->charge!=0){
						exists=true;
						return exists;
					}
					if(daughter->CheckForDaughters(codecheck)){
						exists=true;
						return exists;
					}
				}
				++bpos;
			}while(bpos!=branchlist.end());
		}
	}
	return exists;
}

void CResInfo::DecayGetResInfoPtr_minmass(int &nbodies,array<CResInfo *,5> &daughterresinfo){
	nbodies=bptr_minmass->resinfoptr.size();
	for(int ibody=0;ibody<nbodies;ibody++){
		daughterresinfo[ibody]=bptr_minmass->resinfoptr[ibody];
	}
}

CBranchInfo::CBranchInfo(){
}

void CResInfo::Print(){
	printf("+++++++ ID=%d, M=%g, M_min=%g, %s +++++++++\n",code,mass,minmass,name.c_str());
	printf("Gamma=%g, Spin=%g, Decay=%d\n",width,spin,int(decay));
	printf("Q=%d, B=%d, S=%d, G_parity=%d\n",charge,baryon,strange,G_Parity);
}

void CResList::freegascalc_onespecies_offshell(CResInfo *resinfo,double T,double &epsilon,double &P,double &dens,double &sigma2,double &dedt){
	double width=resinfo->width;
	double mass=resinfo->mass;
	double minmass=resinfo->minmass;
	double wtot=0.0,w;
	double m0,dm=0.25*width;
	double dP,depsilon,ddens,dsigma2,ddedt;
	epsilon=P=dens=sigma2=dedt=0.0;
	if(width/T<0.01){
		freegascalc_onespecies(mass,T,epsilon,P,dens,sigma2,dedt);
	}
	else{
		m0=mass-4.0*width;
		if(m0<minmass)
			m0=minmass;
		while(m0<mass+4.0*width){
			Misc::Pause();
		}
	}
}

void CResList::freegascalc_onespecies(double m,double T,double &epsilon,double &P,double &dens,double &sigma2,double &dedt){
	const double prefactor=1.0/(2.0*PI*PI*pow(HBARC,3));
	double k0,k1,z,k0prime,k1prime,m2,m3,m4,t2,t3,I1,I2,Iomega;
	m2=m*m;
	m3=m2*m;
	m4=m2*m2;
	t2=T*T;
	t3=t2*T;
	z=m/T;
	if(z>1000.0){
		P=epsilon=dens=dedt=0.0;
		printf("z is huge=%g, m=%g, t=%g\n",z,m,T);
	}
	else{
		if(z<0.0){
			printf("___z=%g,m=%g,T=%g ___\n",z,m,T);
			exit(1);
		}
		k0=Bessel::K0(z);
		k1=Bessel::K1(z);
		P=prefactor*(m2*t2*k0+2.0*m*t3*k1);
		epsilon=prefactor*(3.0*m2*t2*k0+(m3*T+6.0*m*t3)*k1);
		dens=P/T;
		k0prime=-k1;
		k1prime=-k0-k1/z;
		dedt=prefactor*(6.0*m2*T*k0+(m3+18.0*m*t2)*k1-3.0*m3*k0prime-((m4/T)+6.0*m2*T)*k1prime);
		Iomega=exp(-z)/(30.0*PI*PI*HBARC*HBARC*HBARC);
		I1=pow(m,1.5)*pow(T,3.5)*7.5*sqrt(2.0*PI);
		I2=24.0*pow(T,5);
		sigma2=Iomega*(I1+I2+0.5*sqrt(I1*I2));  // this is an approximation (+/-2%) to messy integral
	}
}

void CResList::ReadResInfo(){
	CMerge *merge;
	int mothercode,code,decay,strange,charge,baryon,NResonances;
	double mass,mothermass,spin,width,bsum,netm,qR2;
	int ires,jres,ires1,ires2,iresflip,ichannel,nchannels,ibody,nbodies,length, LDecay, i_inel;
	int netq,netb,nets, netg, G_Parity;
	string name, filename;
	CResInfo *resinfoptr=NULL,*oldresinfoptr=NULL, *resinfoptr_1 = NULL;
	CBranchInfo *bptr=NULL,*oldbptr=NULL;
	FILE *resinfofile;
	FILE * decayinfofile;
	char dummy[200],cname[200];
	
	filename=parameter::getS(*parmap,"B3D_RESONANCES_INFO_FILE",string("resinfo/resonances_standardhadrons.dat"));
	printf("will read res info from %s\n",filename.c_str());
	resinfofile=fopen(filename.c_str(),"r");
	fgets(dummy,200,resinfofile);
	fgets(dummy,200,resinfofile);
	fgets(dummy,200,resinfofile);
	fscanf(resinfofile,"%d",&NResonances);
	fgets(dummy,200,resinfofile);
	printf("NResonances=%d\n",NResonances);
	MergeArray=new CMerge **[NResonances];
	SigmaMaxArray=new double *[NResonances];
	for(ires=0;ires<NResonances;ires++){
		MergeArray[ires]=new CMerge *[NResonances];
		SigmaMaxArray[ires]=new double[NResonances];
		for(jres=0;jres<NResonances;jres++){
			MergeArray[ires][jres]=NULL;
			SigmaMaxArray[ires][jres]=0.0;
		}
	}
	for(ires=0;ires<NResonances;ires++){
		resinfoptr=new CResInfo();
		fscanf(resinfofile,"%d %lf %d %d %d %lf %d %d %lf", &resinfoptr->code,&resinfoptr->mass,&resinfoptr->charge,&resinfoptr->baryon, &resinfoptr->strange,&resinfoptr->spin,&resinfoptr->G_Parity,&decay,&resinfoptr->width);
		fgets(cname,100,resinfofile);
		cname[int(strlen(cname))-1]='\0';
		resinfoptr->name=cname;
		resinfoptr->decay=bool(decay);
		resinfoptr->ires=ires;
		resinfoptr->branchlist.clear();
		resmap.insert(CResInfoPair(resinfoptr->code,resinfoptr));
	}			
	fclose(resinfofile);

	filename=parameter::getS(*parmap,"B3D_RESONANCES_DECAYS_FILE",string("resinfo/decays_pdg_weak.dat"));
	printf("will read decay info from %s\n",filename.c_str());
	decayinfofile=fopen(filename.c_str(),"r");
	while(fscanf(decayinfofile,"%d %lf",&mothercode,&mothermass) && !feof(decayinfofile)){
		fgets(dummy,200,decayinfofile);
		fscanf(decayinfofile,"%d %d",&mothercode,&nchannels);
		resinfoptr=GetResInfoPtr(mothercode);
		resinfoptr->minmass=1.0E10;
		bsum=0.0;
		for(ichannel=0;ichannel<nchannels;ichannel++){
			bptr=new CBranchInfo();
			bptr->resinfoptr.clear();
			resinfoptr->branchlist.push_back(bptr);
			fscanf(decayinfofile,"%d",&nbodies);
			netq=-resinfoptr->charge;
			netb=-resinfoptr->baryon;
			nets=-resinfoptr->strange;
			netm=0.0;
			for(ibody=0;ibody<nbodies;ibody++){
				fscanf(decayinfofile,"%d",&code);
				bptr->resinfoptr.push_back(GetResInfoPtr(code));
				netq+=bptr->resinfoptr[ibody]->charge;
				netb+=bptr->resinfoptr[ibody]->baryon;
				nets+=bptr->resinfoptr[ibody]->strange;
				netm+=bptr->resinfoptr[ibody]->mass;
			}	
			//if the total mass is smaller than the minimum required mass, replace it
			if(netm<resinfoptr->minmass){
				resinfoptr->minmass=netm;
				resinfoptr->bptr_minmass=bptr;
			}
			
			//total charge and baryon number should be conserved, and shouldn't be larger than single strangeness
			if(netq!=0 || netb!=0 || abs(nets)>1){
				printf("Charge conservation failure while reading decay info,\nnetq=%d, netb=%d, nets=%d\n",netq,netb,nets);
				printf("MOTHER (ichannel=%d, nbodies=%d):\n",ichannel,nbodies);
				resinfoptr->Print();
				printf("DAUGHTERS:\n");
				for(ibody=0;ibody<nbodies;ibody++)
					bptr->resinfoptr[ibody]->Print();
				if(netq!=0 || netb!=0)
					exit(1);
			}
			fscanf(decayinfofile,"%lf %d",&bptr->branching,&LDecay);
			//store two body decays only
			if(nbodies==2){
				ires1=bptr->resinfoptr[0]->ires;
				ires2=bptr->resinfoptr[1]->ires;
				if(ires1>ires2){
					iresflip=ires1; ires1=ires2; ires2=iresflip;
				}
				merge=MergeArray[ires1][ires2];
				if(merge==NULL){
					MergeArray[ires1][ires2]=new CMerge(resinfoptr,bptr->branching, LDecay);
				}
				else{
					while(merge->next!=NULL){
						merge=merge->next;
					}
					merge->next=new CMerge(resinfoptr,bptr->branching, LDecay);
				}
				if(resinfoptr->mass>bptr->resinfoptr[0]->mass+bptr->resinfoptr[1]->mass){
					qR2=Misc::triangle(resinfoptr->mass,
					bptr->resinfoptr[0]->mass,bptr->resinfoptr[1]->mass);
					SigmaMaxArray[ires1][ires2]+=
						(bptr->branching*4.0*PI*HBARC*HBARC)*(2.0*resinfoptr->spin+1.0)/
						((2.0*bptr->resinfoptr[0]->spin+1.0)*(2.0*bptr->resinfoptr[1]->spin+1.0)*qR2);
					if(ires2!=ires1)
						SigmaMaxArray[ires2][ires1]=SigmaMaxArray[ires1][ires2];
				}
			}
			bsum+=bptr->branching;
		}
	}
	fclose(decayinfofile);
}

CResInfo* CResList::GetResInfoPtr(int code){
	CResInfoMap::iterator rpos;
	rpos=resmap.find(code);
	return rpos->second;
}

void CResList::CalcEoS(double T0,double Tf,double delT){
	CResInfo *resinfoptr;
	CResInfoMap::iterator rpos;
	printf("#_____________________\n#  T       s         P        epsilon\n");
	double T,P,epsilon,s,m,degen;
	double pi,epsiloni,densi,sigma2i,dedti,si;
	for(T=T0;T<Tf+0.00000001;T+=delT){
		P=epsilon=s=0.0;
		for(rpos=resmap.begin();rpos!=resmap.end();rpos++){
			resinfoptr=rpos->second;
			if(resinfoptr->code!=22){
				degen=2.0*resinfoptr->spin+1.0;
				m=resinfoptr->mass;
				freegascalc_onespecies(m,T,epsiloni,pi,densi,sigma2i,dedti);
				P+=pi*degen;
				epsilon+=epsiloni*degen;
				s+=(pi+epsiloni)*degen/T;
			}
		}
		printf("%6.2f %15.10e %15.10e %15.10e\n",T,s,P,epsilon);
	}
}

void CResList::CalcEoS(double T,double &epsilon,double &P,double &nhadrons,
vector<double> &density,vector<double> &boseweight){
	CResInfo *resinfoptr;
	CResInfoMap::iterator rpos;
	double s,m,degen;
	double pi,epsiloni,densi,sigma2i,dedti,bosenorm=0.0;
	int ires=0,nres,ibose,nbose;
	char dummy[100];
	if(boseweight.size()==0){
		boseweight.resize(2);
	}
	P=epsilon=nhadrons=0.0;
	density.clear();
	nres=resmap.size();
	if(GetResInfoPtr(22)->code==22)
		nres-=1;
	density.resize(nres);
	for(rpos=resmap.begin();rpos!=resmap.end();rpos++){
		resinfoptr=rpos->second;
		if(resinfoptr->code!=22){
			degen=2.0*resinfoptr->spin+1.0;
			m=resinfoptr->mass;
			nbose=1;
			if(abs(resinfoptr->code)==211 || resinfoptr->code==111)
				nbose=boseweight.size()-1;
			if(nbose<1)
				nbose=1;
			density[ires]=0.0;
			boseweight[0]=0.0;
			for(ibose=1;ibose<=nbose;ibose++){
				freegascalc_onespecies(m,T/double(ibose),epsiloni,pi,densi,sigma2i,dedti);
				if(resinfoptr->code==211)
					boseweight[ibose]=boseweight[ibose-1]+densi*degen;
				P+=pi*degen;
				epsilon+=epsiloni*degen;
				density[ires]+=densi*degen;
				nhadrons+=density[ires];
			}
			ires+=1;
		}
	}
	for(ibose=1;ibose<=boseweight.size();ibose++)
		boseweight[ibose]/=boseweight[nbose];
	//printf("T=%g, epsilon=%g, P=%g, nhadrons=%g\n",T,epsilon,P,nhadrons);
}

void CResList::CalcEoS(double T,double &epsilon,double &P,double &nhadrons,double &cs2,vector<double> &density){
	CResInfo *resinfoptr;
	CResInfoMap::iterator rpos;
	double s,m,degen,dedt;
	double pi,epsiloni,densi,sigma2i,dedti;
	int ires=0,nres;
	char dummy[100];
	P=epsilon=nhadrons=dedt=0.0;
	density.clear();
	nres=resmap.size();
	if(GetResInfoPtr(22)->code==22)
		nres-=1;
	density.resize(nres);
	for(rpos=resmap.begin();rpos!=resmap.end();rpos++){
		resinfoptr=rpos->second;
		if(resinfoptr->code!=22){
			degen=2.0*resinfoptr->spin+1.0;
			m=resinfoptr->mass;
			freegascalc_onespecies(m,T,epsiloni,pi,densi,sigma2i,dedti);
			P+=pi*degen;
			epsilon+=epsiloni*degen;
			density[ires]=densi*degen;
			dedt+=dedti*degen;
			nhadrons+=density[ires];
			ires+=1;
		}
	}
	s=(P+epsilon)/T;
	cs2=s/dedt;
	//printf("T=%g, epsilon=%g, P=%g, nhadrons=%g\n",T,epsilon,P,nhadrons);
}

void CResList::CalcEoSandChi(double T){
	CResInfo *resinfoptr;
	CResInfoMap::iterator rpos;
	double P,epsilon,s,m,degen;
	double Qu,Qd,Qs,Q,S,B,chi[3][3],q[3];
	double pi,epsiloni,densi,sigma2i,dedti,si;
	char dummy[100];
	int a,b;
	for(a=0;a<3;a++){
		for(b=0;b<3;b++) chi[a][b]=0.0;
	}
	P=epsilon=s=0.0;
	for(rpos=resmap.begin();rpos!=resmap.end();rpos++){
		resinfoptr=rpos->second;
		if(resinfoptr->code!=22){
			degen=2.0*resinfoptr->spin+1.0;
			m=resinfoptr->mass;
			freegascalc_onespecies(m,T,epsiloni,pi,densi,sigma2i,dedti);
			P+=pi*degen;
			epsilon+=epsiloni*degen;
			s+=(pi+epsiloni)*degen/T;
			Q=resinfoptr->charge;
			B=resinfoptr->baryon;
			S=resinfoptr->strange;
			q[0]=(B+Q);
			q[1]=(2*B+S-Q);
			q[2]=-S;
			//sprintf(dummy,"%s",resinfoptr->name.c_str());
			//printf("%40s: u=%3.1f,  d=%3.1f,  s=%3.1f\n",dummy,q[0],q[1],q[2]);
			for(a=0;a<3;a++){
				for(b=0;b<3;b++) chi[a][b]+=densi*q[a]*q[b];
			}
		}
	}
	printf("T=%6.2f  s/T^3=%10.4e\n",T,s*HBARC*HBARC*HBARC/(T*T*T));
	printf("---- chi/s -----\n");
	for(a=0;a<3;a++){
		for(b=0;b<3;b++) printf("%7.4f ",chi[a][b]/s);
		printf("\n");
	}
	printf("-------------------\n");
}

void CResList::CalcEosandKubo(double T,double &epsilon,double &P,double &nhadrons,double &sdens,double &kubo){
	CResInfo *resinfoptr;
	CResInfoMap::iterator rpos;
	double m,degen,dedt;
	double pi,epsiloni,densi,sigma2i,dedti;
	double e,p,delp=2.0;
	int ires=0,nres;
	char dummy[100];
	vector<double> density;
	P=epsilon=nhadrons=dedt=kubo=0.0;
	density.clear();
	nres=resmap.size();
	if(GetResInfoPtr(22)->code==22)
		nres-=1;
	density.resize(nres);
	for(rpos=resmap.begin();rpos!=resmap.end();rpos++){
		resinfoptr=rpos->second;
		if(resinfoptr->code!=22){
			degen=2.0*resinfoptr->spin+1.0;
			m=resinfoptr->mass;
			freegascalc_onespecies(m,T,epsiloni,pi,densi,sigma2i,dedti);
			P+=pi*degen;
			epsilon+=epsiloni*degen;
			density[ires]=densi*degen;
			nhadrons+=density[ires];
			for(p=0.5*delp;p<20*sqrt(m*T);p+=delp){
				e=sqrt(m*m+p*p);
				kubo+=degen*(pow(p,6)/(e*e))*exp(-e/T)*delp;
			}
			ires+=1;
		}
	}
	kubo=kubo/(30.0*PI*PI*T*pow(HBARC,4));
	sdens=(P+epsilon)/T;
}

#endif
