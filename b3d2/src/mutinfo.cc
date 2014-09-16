#include "mutinfo.h"

CB3D *CMuTInfo::b3d=NULL;
int CMuTInfo::NETEVENTS=0;
int CMuTInfo::NTAU=0;
double CMuTInfo::DELTAU=0.0;

CMuTInfo::CMuTInfo(){
	int itau;
	Pxpi.resize(NTAU);
	Pypi.resize(NTAU);
	Epi.resize(NTAU);
	Mpi.resize(NTAU);
	Npi.resize(NTAU);
	PxK.resize(NTAU);
	PyK.resize(NTAU);
	EK.resize(NTAU);
	MK.resize(NTAU);
	NK.resize(NTAU);
	Zero();
}

void CMuTInfo::Zero(){
	int itau;
	for(itau=0;itau<NTAU;itau++){
		Pxpi[itau]=Pypi[itau]=Epi[itau]=PxK[itau]=PyK[itau]=EK[itau]=0.0;
		Npi[itau]=NK[itau]=0;
	}
	muK=mupi=0.0;
	TK=Tpi=150.0;
}

void CMuTInfo::UpdateNMPE(CB3DCell *cell){
	double e;
	CPartMap::iterator ppos;
	CPart *part;
	int itau;
	itau=lrint(b3d->tau/DELTAU)-1;
	for(ppos=cell->partmap.begin();ppos!=cell->partmap.end();ppos++){
		part=ppos->second;
		if(part->resinfo->code==111 || abs(part->resinfo->code)==211){
			Npi[itau]+=1;
			Mpi[itau]+=part->resinfo->mass;
			Pxpi[itau]+=part->p[1];
			Pypi[itau]+=part->p[2];
			e=cosh(part->eta)*part->p[0]-sinh(part->eta)*part->p[3];
			Epi[itau]+=e;
		}
		else if(abs(part->resinfo->code)==321 || abs(part->resinfo->code)==311){
			NK[itau]+=1;
			MK[itau]+=part->resinfo->mass;
			PxK[itau]+=part->p[1];
			PyK[itau]+=part->p[2];
			e=cosh(part->eta)*part->p[0]-sinh(part->eta)*part->p[3];
			EK[itau]+=e;
		}
	}
}

void CMuTInfo::MuTCalc(){
	int ntry,itau;
	double Tguess,Etot,u0,e,rho0,rhotarget,etarget,dedT,drhodT,accuracy;
	double dT,de,ddedT,dens,p,sigma2,volume;
	double tau=b3d->tau;
	double mass;
	itau=lrint(b3d->tau/DELTAU)-1;
	FindMuTInfo_pi(itau);
	FindMuTInfo_K(itau);
}
	
void CMuTInfo::FindMuTUxUy(double tau,int N,double E,double M,double Px,double Py,double degen,double &T,double &mu,double &ux,double &uy){
	double ur,u0,e,dedT,Tguess,p,dens,sigma2,ddedT,rho0,Etot,rhotarget,etarget;
	double mass,volume,drhodT,accuracy,de,dT;
	int ntry;
	Tguess=T;
	if(Tguess!=Tguess || Tguess<20.0 || Tguess>180.0)
		Tguess=150.0;
	if(N==0){
		mu=T=ux=uy=0.0;
	}
	else if(N==1){
		ux=Px/E;
		uy=Py/E;
		u0=1.0/sqrt(1.0-ux*ux-uy*uy);
		ux*=u0;
		uy*=u0;
		mu=T=0.0;
	}
	else{
		ux=Px/E;
		uy=Py/E;
		u0=1.0/sqrt(1.0-ux*ux-uy*uy);
		ux*=u0;
		uy*=u0;
		Etot=sqrt(E*E-Px*Px-Py*Py);
		volume=u0*b3d->DXY*b3d->DXY*2.0*b3d->ETAMAX*tau;
		mass=M/N;
		etarget=Etot/N;
		rhotarget=N/(volume*NETEVENTS);
			
		if(((etarget-mass)/mass)<0.1){
			T=2.0*(etarget-mass)/3.0;
			if(T<5.0)
				T=5.0;
			b3d->reslist->freegascalc_onespecies(mass,T,de,p,dens,sigma2,ddedT);
			rho0=degen*dens;
		}
		else{
			ntry=0;
			do{
				ntry+=1;
				e=dedT=drhodT=rho0=0.0;
				accuracy=0.0;
				b3d->reslist->freegascalc_onespecies(mass,Tguess,de,p,dens,sigma2,ddedT);
				rho0=degen*dens;
				e=degen*de;
				dedT=degen*ddedT;
				drhodT=degen*de/(Tguess*Tguess);
				e=e/rho0;
				accuracy=fabs((e-etarget)/mass);
				dT=(etarget-e)/(dedT/rho0-e*drhodT/rho0);
				if(fabs(dT)>0.5*Tguess)
					dT=0.5*Tguess*dT/fabs(dT);
				Tguess+=dT;
			}while(accuracy>1.0E-6);
			T=Tguess;
		}
		mu=T*log(rhotarget/rho0);
	}
}

void CMuTInfo::FindMuTInfo_pi(int itau){
	double tau=(1+itau)*DELTAU;
	FindMuTUxUy(tau,Npi[itau],Epi[itau],Mpi[itau],Pxpi[itau],Pypi[itau],3,Tpi,mupi,uxpi,uypi);
}
void CMuTInfo::FindMuTInfo_K(int itau){
	double tau=(1+itau)*DELTAU;
	FindMuTUxUy(tau,NK[itau],EK[itau],MK[itau],PxK[itau],PyK[itau],4,TK,muK,uxK,uyK);
}

void CMuTInfo::Print(){
	int itau=lrint(b3d->tau/DELTAU)-1;
	printf("-------- MuT Info, itau=%d ----------\n",itau);
	printf("Npi=%d, Epi/N=%g, Mpi/Npi=%g, Pxpi/Npi=%g, Pypi/Npi=%g\n",
	Npi[itau],Epi[itau]/Npi[itau],Mpi[itau]/Npi[itau],Pxpi[itau]/Npi[itau],Pypi[itau]/Npi[itau]);
	printf("NK=%d, EK/N=%g, MK/NK=%g, PxK/NK=%g, PyK/NK=%g\n",
	NK[itau],EK[itau]/NK[itau],MK[itau]/NK[itau],PxK[itau]/NK[itau],PyK[itau]/NK[itau]);
	printf("Tpi=%g, mupi=%g, uxpi=%g, uypi=%g\n",Tpi,mupi,uxpi,uypi);
	printf("TK =%g, muK= %g, uxK =%g, uyK =%g\n",TK,muK,uxK,uyK);
}

/*
void CMuTInfo::MuCalc_PionsWithBose(){
int nbose=b3d->NBOSE;
double dalpha,alpha=0.0;
double delr2=CLocalInfo::DELR2;
int nrbins=CLocalInfo::NRBINS;
double r2max=nrbins*delr2;
double rmax=sqrt(r2max);
int ir2,netevents=CLocalInfo::NETEVENTS,ntry;
double Tguess,Etot,u0,rhotarget,etarget,dedT,ddedT,drhodT,accuracy,dT;
double dedalpha,drhodalpha,denom,dele,delrho,e,rho,x,de,dens,sigma2,p,volume;
double tau=b3d->tau,mass;
int n;

for(ir2=0;ir2<nrbins;ir2++){
if(N>0){
if(T>0)
Tguess=T;
else
Tguess=150.0;
ur=Pr/E;
ur=ur/sqrt(1.0-ur*ur);
if(N>1){
Etot=sqrt(E*E-Pr*Pr);
u0=sqrt(1.0+ur*ur);
volume=u0*PI*delr2*2.0*b3d->ETAMAX*tau;
etarget=Etot/(netevents*volume);
rhotarget=N/(netevents*volume);
etarget=etarget/rhotarget;
mass=M/N;
ntry=0;
do{
ntry+=1;
accuracy=0.0;
e=rho=dedT=drhodT=dedalpha=drhodalpha=0.0;
for(n=1;n<=nbose;n++){
b3d->reslist->freegascalc_onespecies(mass,Tguess/double(n),de,p,dens,sigma2,ddedT);
x=exp(n*alpha);
e+=degen*de*x;
rho+=degen*dens*x;
dedT+=degen*x*ddedT/double(n);
drhodT+=degen*x*n*de/(Tguess*Tguess);
dedalpha+=degen*de*x*n;
drhodalpha+=degen*dens*x*n;
}
dedT=(dedT/rho)-e*drhodT/(rho*rho);
dedalpha=(dedalpha/rho)-e*drhodalpha/(rho*rho);
e=e/rho;
denom=dedT*drhodalpha-dedalpha*drhodT;
dele=etarget-e;
delrho=rhotarget-rho;
if(ntry>30){
printf("MuCalc:: ntry=%d, ir1=%d, etarget=%g, dele=%g, rhotarget=%g, delrho=%g, T=%g, mu=%g\n",
ntry, ir2,etarget,dele,rhotarget,delrho,Tguess,alpha*Tguess);
//Misc::Pause();
}
accuracy=fabs(dele/etarget)+fabs(delrho/rhotarget);
dT=(drhodalpha*dele-dedalpha*delrho)/denom;
if(fabs(dT)>0.5*Tguess){
dT=0.5*dT*Tguess/fabs(dT);
}
Tguess+=dT;
dalpha=(dedT*delrho-drhodT*dele)/denom;
if(fabs(dalpha)>0.5){
dalpha=dalpha*0.5/fabs(dalpha);
}
alpha+=dalpha;
}while(accuracy>0.01*rhotarget);
T=Tguess;
mu=alpha*T;
}
}
}
}
*/
