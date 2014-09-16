#ifndef __RESONANCES_CC__
#define __RESONANCES_CC__
#include "b3d.h"

using namespace std;

CB3D *CResList::b3d=NULL;

CResList::CResList(){
	if(b3d!=NULL){
		parmap=&(b3d->parmap);
		ReadResInfo();
	}
}

CResList::CResList(parameterMap* parmap_in){
	parmap=parmap_in;
	ReadResInfo();
}

CRandom *CResInfo::ranptr=new CRandom(-1234);

CResInfo::CResInfo(){
	count=0;
	minmass=0.0;
	nextResInfoptr=NULL;
	firstbptr=NULL;
}

CMerge::CMerge(CResInfo *resinfo_in,double branching_in, int L_in){
	resinfo=resinfo_in;
	branching=branching_in;
	L = L_in;
	next=NULL;
}

void CResInfo::DecayGetResInfoptr(int &nbodies,CResInfo **&daughterresinfoptr){
	double randy,bsum;
	int ibody;
	CBranchInfo *bptr;
	bptr=NULL;
	bsum=0.0;
	randy=ranptr->ran();

	do{
		if(bptr==NULL){
			bptr=firstbptr;
		}
		else{
			bptr=bptr->nextbptr;
		}
		bsum+=bptr->branching;
		if(bsum-1.0>1.0E-6){
			cout << "FATAL: In DecayGetResInfo: bsum too large, = " << bsum << endl;
			exit(1);
		}

	}while(bsum<randy);

	nbodies=bptr->nbodies;
	for(ibody=0;ibody<nbodies;ibody++){
		daughterresinfoptr[ibody]=bptr->resinfoptr[ibody];
	}

}

bool CResInfo::CheckForDaughters(int codecheck){//checks to see if any decay daughters match code
	int ibody,nbodies;
	bool exists=false;
	CResInfo *daughter;
	CBranchInfo *bptr;
	if(codecheck!=0){
		if(code==codecheck){
			exists=true;
			return exists;
		}
		if(decay){
			bptr=firstbptr;
			do{
				nbodies=bptr->nbodies;
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
				bptr=bptr->nextbptr;
			}while(bptr!=NULL);
		}
	}
	else{
		if(charge!=0){
			exists=true;
			return exists;
		}
		if(decay){
			bptr=firstbptr;
			do{
				nbodies=bptr->nbodies;
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
				bptr=bptr->nextbptr;
			}while(bptr!=NULL);
		}
	}
	return exists;
}

void CResInfo::DecayGetResInfoptr_minmass(int &nbodies,CResInfo **&daughterresinfoptr){
	nbodies=bptr_minmass->nbodies;
	for(int ibody=0;ibody<nbodies;ibody++){
		daughterresinfoptr[ibody]=bptr_minmass->resinfoptr[ibody];
	}

}

CBranchInfo::CBranchInfo(){
	nextbptr=NULL;
}

void CResInfo::Print(){
	cout << "+++++++ " << name << " +++++++" << endl;
	cout << "Code= " << code;
	cout << ", Mass= " << mass;
	cout << ", Minmass= " << minmass;
	cout << ", Width= " << width;
	cout << ", Spin= " << spin;
	cout << ", Charge= " << charge;
	cout << ", Baryon= " << baryon;
	cout << ", Strange= " << strange;
	cout << ", G_Parity= ";
	if(G_Parity){
		cout <<  G_Parity;
	}else{
		cout << "NULL";
	}
	cout << ", Decay= " << decay << endl;
}

/*
This method reads in the resonance information stored in a given file. The information has been preformatted, therefore,
this simply reads in the uncommented preformatted data.
*/

void CResList::ReadResInfo(){
	CMerge *merge;
	int mothercode,code,decay,strange,charge,baryon;
	double mass,mothermass,spin,width,bsum,netm;
	int ires,ires1,ires2,iresflip,ichannel,nchannels,ibody,nbodies,length, LDecay, i_inel;
	int netq,netb,nets, netg, G_Parity;
	string dummy, name, filename;
	stringstream dummystream;
	CResInfo *resinfoptr=NULL,*oldresinfoptr=NULL, *resinfoptr_1 = NULL;
	CBranchInfo *bptr=NULL,*oldbptr=NULL;
	ifstream resinfofile,decayinfofile;
	stringstream sst;
	char dummychars[200];
	int foobar = 0;
	
	//read in the filename from the parameter map, or use a default string
	filename=parameter::getS(*parmap,"B3D_RESONANCES_INFO_FILE",string("resinfo/resonances_standardhadrons.dat"));
	cout << "will read res info from " << filename << endl;
	resinfofile.open(filename.c_str());
	if(resinfofile){
		//Comments may be prepended to the resonance info file, as long as they are led by a # character
		getline(resinfofile, dummy);
		while(dummy.at(0)=='#'){
			getline(resinfofile,dummy);
		}
		dummystream << dummy;

		dummystream >> NResonances; //the first line is the number of resonances
		printf("NResonances=%d\n",NResonances);
		MergeArray=new CMerge **[NResonances];
		for(ires=0;ires<NResonances;ires++){
			MergeArray[ires]=new CMerge *[NResonances];
			for(ires2=0;ires2<NResonances;ires2++){
				MergeArray[ires][ires2]=NULL;
			}
		}

		for(ires=0;ires<NResonances;ires++){
			resinfofile >> code >> mass >> charge >> baryon >> strange >> spin >> G_Parity 
			>> decay >> width;
			char cname[100];
			resinfofile.getline(cname,100);
			name=cname;
			resinfoptr=new CResInfo();
			if(oldresinfoptr==NULL){
				GfirstResInfoptr=resinfoptr;
			}
			else{
				oldresinfoptr->nextResInfoptr=resinfoptr;
			}
			resinfoptr->ires=ires;
			resinfoptr->code=code;
			resinfoptr->mass=mass;
			resinfoptr->charge=charge;
			resinfoptr->baryon=baryon;
			resinfoptr->strange=strange;
			resinfoptr->spin=spin;
			resinfoptr->decay=bool(decay);
			resinfoptr->width=width;    
			resinfoptr->name=name;

			if(G_Parity != 9){
				resinfoptr->G_Parity = G_Parity;
			}
			else{
				resinfoptr->G_Parity = 9;
			}
			oldresinfoptr=resinfoptr;
		}

		resinfofile.close();
	}else{
		cout << "Unable to open resonance info file: " << filename << endl;
		exit(-1);
	}
	//cout << "RESONANCES READ, WILL BEGIN READING DECAY INFO\n";

	filename=parameter::getS(*parmap,"B3D_RESONANCES_DECAYS_FILE",string("resinfo/decays_weak.dat"));
	cout << "will read decay info from " << filename << endl;
	decayinfofile.open(filename.c_str());
	if(decayinfofile) {
		while(decayinfofile >> mothercode >> mothermass){ //read in the mother's particle code and mass...
			getline(decayinfofile, dummy); // and ignore the rest
		decayinfofile >> mothercode >> nchannels;
		GetResInfoptr(mothercode,resinfoptr);
		resinfoptr->minmass=1.0E10;
		bsum=0.0;

		for(ichannel=0;ichannel<nchannels;ichannel++){
			bptr=new CBranchInfo();
			if(ichannel==0){
				resinfoptr->firstbptr=bptr;
			}
			else{
				oldbptr->nextbptr=bptr;
			}
			decayinfofile >> bptr->nbodies;
			netq=-resinfoptr->charge;
			netb=-resinfoptr->baryon;
			nets=-resinfoptr->strange;
			netm=0.0;

			for(ibody=0;ibody<bptr->nbodies;ibody++){
				decayinfofile >> code;
				GetResInfoptr(code,bptr->resinfoptr[ibody]);
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
				printf("MOTHER (ichannel=%d, nbodies=%d):\n",ichannel,bptr->nbodies);
				resinfoptr->Print();
				printf("DAUGHTERS:\n");
				for(ibody=0;ibody<bptr->nbodies;ibody++) bptr->resinfoptr[ibody]->Print();
					if(netq!=0 || netb!=0) exit(1);
			}
			decayinfofile >> bptr->branching;
			decayinfofile >> LDecay;
				//store two body decays only
			if(bptr->nbodies==2){
				ires1=bptr->resinfoptr[0]->ires;
				ires2=bptr->resinfoptr[1]->ires;
				if(ires1>ires2){
					iresflip=ires1; ires1=ires2; ires2=iresflip;
				}
				
				merge=MergeArray[ires1][ires2];
				foobar++;
				if(merge==NULL){
					MergeArray[ires1][ires2]=new CMerge(resinfoptr,bptr->branching, LDecay);
				}
				else{
					while(merge->next!=NULL){
						merge=merge->next;
					}
					merge->next=new CMerge(resinfoptr,bptr->branching, LDecay);
				}
			}
			
			bsum+=bptr->branching;
			oldbptr=bptr;
		}
	}
	decayinfofile.close();
}else{
	cout << "Unable to open decay information file " << filename << endl;
	exit(-1);
}
cout << "Done reading in decay and resonance info, read in " << foobar << " decays." << endl;
}

void CResList::GetResInfoptr(int code,CResInfo *&resinfoptr){
	resinfoptr=GfirstResInfoptr;
	while(resinfoptr->code!=code){
		resinfoptr=resinfoptr->nextResInfoptr;
		if(resinfoptr==NULL){
			cout << "FATAL: failure to identify code=" << code << endl;
			exit(1);
		}
	}
}

void CResList::PrintYields(){
	CResInfo *resinfoptr;
	cout << "_____________________\n     YIELDS\n";
	resinfoptr=GfirstResInfoptr;
	while(resinfoptr!=NULL){
		if(resinfoptr->count!=0) cout << resinfoptr->name << " " << resinfoptr->code << " " << resinfoptr->mass<< " : " << resinfoptr->count << endl;
		resinfoptr=resinfoptr->nextResInfoptr;
	}
}

void CResList::CalcEoS(double T0,double Tf,double delT){
	CResInfo *resinfoptr;
	cout << "#_____________________\n#  T       s         P        epsilon\n";
	double T,P,epsilon,s,m,degen;
	double pi,epsiloni,densi,sigma2i,dedti,si;
	for(T=T0;T<Tf+0.00000001;T+=delT){
		P=epsilon=s=0.0;
		resinfoptr=GfirstResInfoptr;
		while(resinfoptr!=NULL){
			if(resinfoptr->code!=22){
				degen=2.0*resinfoptr->spin+1.0;
				m=resinfoptr->mass;
				//printf("ID=%d, degen=%g, m=%g, q=%d, s=%d, b=%d\n",resinfoptr->code,degen,m,resinfoptr->charge,resinfoptr->strange,resinfoptr->baryon);
				freegascalc_onespecies(m,T,pi,epsiloni,densi,sigma2i,dedti);
				P+=pi*degen;
				epsilon+=epsiloni*degen;
				s+=(pi+epsiloni)*degen/T;
			}
			resinfoptr=resinfoptr->nextResInfoptr;
		}
		printf("%6.2f %15.10e %15.10e %15.10e\n",T,s,P,epsilon);
	}
}

void CResList::CalcEoSandChi(double T){
	CResInfo *resinfoptr;
	double P,epsilon,s,m,degen;
	double Qu,Qd,Qs,Q,S,B,chi[3][3],q[3];
	double pi,epsiloni,densi,sigma2i,dedti,si;
	char dummy[100];
	int a,b;
	for(a=0;a<3;a++){
		for(b=0;b<3;b++) chi[a][b]=0.0;
	}
P=epsilon=s=0.0;
resinfoptr=GfirstResInfoptr;
while(resinfoptr!=NULL){
	if(resinfoptr->code!=22){
		degen=2.0*resinfoptr->spin+1.0;
		m=resinfoptr->mass;
				//printf("ID=%d, degen=%g, m=%g, q=%d, s=%d, b=%d\n",resinfoptr->code,degen,m,resinfoptr->charge,resinfoptr->strange,resinfoptr->baryon);
		freegascalc_onespecies(m,T,pi,epsiloni,densi,sigma2i,dedti);
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
resinfoptr=resinfoptr->nextResInfoptr;
}
printf("T=%6.2f  s/T^3=%10.4e\n",T,s*HBARC*HBARC*HBARC/(T*T*T));
printf("---- chi/s -----\n");
for(a=0;a<3;a++){
	for(b=0;b<3;b++) printf("%7.4f ",chi[a][b]/s);
		printf("\n");
}
printf("-------------------\n");
}


void CResList::freegascalc_onespecies(double m,double t,double &p,double &e,double &dens,double &sigma2,double &dedt){
	const double prefactor=1.0/(2.0*PI*PI*pow(HBARC,3));
	double k0,k1,z,k0prime,k1prime,m2,m3,m4,t2,t3,I1,I2,Iomega;
	m2=m*m;
	m3=m2*m;
	m4=m2*m2;
	t2=t*t;
	t3=t2*t;
	z=m/t;
	if(z>1000.0){
		p=e=dens=dedt=0.0;
		printf("z is huge=%g, m=%g, t=%g\n",z,m,t);
	}
	else{
		if(z<0.0){
			printf("___z=%g,m=%g,T=%g ___\n",z,m,t);
			exit(1);
		}
		k0=Bessel::K0(z);
		k1=Bessel::K1(z);
		p=prefactor*(m2*t2*k0+2.0*m*t3*k1);
		e=prefactor*(3.0*m2*t2*k0+(m3*t+6.0*m*t3)*k1);
		dens=p/t;
		k0prime=-k1;
		k1prime=-k0-k1/z;
		dedt=prefactor*(6.0*m2*t*k0+(m3+18.0*m*t2)*k1-3.0*m3*k0prime-((m4/t)+6.0*m2*t)*k1prime);
		Iomega=exp(-m/t)/(30.0*PI*PI*HBARC*HBARC*HBARC);
		I1=pow(m,1.5)*pow(t,3.5)*7.5*sqrt(2.0*PI);
		I2=24.0*pow(t,5);
		sigma2=Iomega*(I1+I2+0.5*sqrt(I1*I2));  // this is an approximation (+/-2%) to messy integral
	}
}

#endif
