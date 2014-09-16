#include "mucalc.h"

CB3D *CLocalSpeciesInfo::b3d=NULL;
CB3D *CLocalInfo::b3d=NULL;
int CLocalInfo::NRBINS=0;
int CLocalInfo::NETEVENTS=0;
int CLocalInfo::IMUCALC=0;
double CLocalInfo::DELR2=0.0;
bool CLocalInfo::printing=false;

CLocalSpeciesInfo::CLocalSpeciesInfo(string nameset){
	name=nameset;
	degen=0;
	int nrbins=CLocalInfo::NRBINS;
	T.resize(nrbins);
	mu.resize(nrbins);
	ur.resize(nrbins);
	N.resize(nrbins);
	E.resize(nrbins);
	Pr.resize(nrbins);
	M.resize(nrbins);
}

void CLocalSpeciesInfo::Zero(){
	int ir2;
	for(ir2=0;ir2<CLocalInfo::NRBINS;ir2++){
		T[ir2]=mu[ir2]=ur[ir2]=-1.0;
		Pr[ir2]=E[ir2]=M[ir2]=N[ir2]=0;
	}
}

CLocalInfo::CLocalInfo(int itauset){
	itau=itauset;
	
	pion=new CLocalSpeciesInfo("pion");
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(211,pion));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-211,pion));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(111,pion));
	
	kaon=new CLocalSpeciesInfo("kaon");
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(321,kaon));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-321,kaon));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(311,kaon));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-311,kaon));
	
	nucleon=new CLocalSpeciesInfo("nucleon");
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(2212,nucleon));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-2212,nucleon));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(2112,nucleon));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-2112,nucleon));

	delta=new CLocalSpeciesInfo("delta");
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(1114,delta));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-1114,delta));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(2114,delta));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-2114,delta));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(2214,delta));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-2214,delta));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(2224,delta));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-2224,delta));
		
	lambda=new CLocalSpeciesInfo("lambda");
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(3122,lambda));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-3122,lambda));
	
	sigma=new CLocalSpeciesInfo("sigma");
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(3222,sigma));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-3222,sigma));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(3212,sigma));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-3212,sigma));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(3112,sigma));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-3112,sigma));
	
	sigmastar=new CLocalSpeciesInfo("sigmastar");
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(3224,sigmastar));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-3224,sigmastar));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(3214,sigmastar));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-3214,sigmastar));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(3114,sigmastar));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-3114,sigmastar));
	
	xsi=new CLocalSpeciesInfo("xsi");
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(3312,xsi));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-3312,xsi));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(3322,xsi));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-3322,xsi));
	
	xsistar=new CLocalSpeciesInfo("xsistar");
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(3314,xsistar));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-3314,xsistar));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(3324,xsistar));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-3324,xsistar));
	
	omega=new CLocalSpeciesInfo("omega");
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(3334,omega));
	lsimap.insert(std::pair<int,CLocalSpeciesInfo*>(-3334,omega));
	
	CLSIMap::iterator lsipos;
	CResInfoMap::iterator rpos;
	CResInfoMap *resmap=&(b3d->reslist->resmap);
	CResInfo *resinfo;
	int pid;
	CLocalSpeciesInfo *lsi;
	for(rpos=resmap->begin();rpos!=resmap->end();rpos++){
		resinfo=rpos->second;
		pid=resinfo->code;
		lsipos=lsimap.find(pid);
		if(lsipos!=lsimap.end()){
			lsi=lsipos->second;
			lsi->degen+=2.0*resinfo->spin+1.0;
		}
	}
}

void CLocalInfo::Zero(){
	CLocalSpeciesInfo *lsi;
	CLSIMap::iterator lsipos;
	for(lsipos=lsimap.begin();lsipos!=lsimap.end();++lsipos){
		lsi=lsipos->second;
		lsi->Zero();
	}
}

void CLocalInfo::CheckMap(){
	CLSIMap::iterator lsipos;
	printf("Checking LSI Map for itau=%d\n",itau);
	for(lsipos=lsimap.begin();lsipos!=lsimap.end();++lsipos){
		printf("checking, key=%d, name=%s\n",lsipos->first,lsipos->second->name.c_str());
	}
	lsipos=lsimap.find(-211);
	printf("how about for key=-211, name=%s\n",lsipos->second->name.c_str());
}

void CLocalSpeciesInfo::MuCalc(){
	double delr2=CLocalInfo::DELR2;
	int nrbins=CLocalInfo::NRBINS;
	double r2max=nrbins*delr2;
	double rmax=sqrt(r2max);
	int ir2,netevents=CLocalInfo::NETEVENTS,ntry;
	double Tguess,Etot,u0,e,rho0,rhotarget,etarget,dedT,drhodT,accuracy;
	double dT,de,ddedT,dens,p,sigma2,volume;
	double tau=b3d->tau;
	double mass;
	
	for(ir2=0;ir2<nrbins;ir2++){
		if(T[ir2]>0)
			Tguess=T[ir2];
		else
			Tguess=150.0;
		if(N[ir2]>10){
			ur[ir2]=Pr[ir2]/E[ir2];
			ur[ir2]=ur[ir2]/sqrt(1.0-ur[ir2]*ur[ir2]);
		}
		if(N[ir2]>1){
			Etot=sqrt(E[ir2]*E[ir2]-Pr[ir2]*Pr[ir2]);
			u0=sqrt(1.0+ur[ir2]*ur[ir2]);
			volume=u0*PI*delr2*2.0*b3d->ETAMAX*tau;
			
			Etot=Etot/(volume*netevents);
			rhotarget=N[ir2]/(volume*netevents);
			etarget=Etot/rhotarget;
			mass=M[ir2]/N[ir2];
			if(etarget<mass){
				printf("etarget=%g, mass=%g, N[%d]=%d\n",etarget,mass,ir2,N[ir2]);
				exit(1);
			}
			if(((etarget-mass)/mass)<0.1){
				T[ir2]=2.0*(etarget-mass)/3.0;
				if(T[ir2]<5.0)
					T[ir2]=5.0;
				b3d->reslist->freegascalc_onespecies(mass,T[ir2],de,p,dens,sigma2,ddedT);
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
					accuracy=fabs(e-etarget);
					dT=(etarget-e)/(dedT/rho0-e*drhodT/rho0);
					if(fabs(dT)>0.5*Tguess)
						dT=0.5*Tguess*dT/fabs(dT);
					Tguess+=dT;
				}while(accuracy>0.1);
				T[ir2]=Tguess;
			}
			mu[ir2]=T[ir2]*log(rhotarget/rho0);
		}
	}
}

void CLocalSpeciesInfo::MuCalc_PionsWithBose(){
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
		if(N[ir2]>0){
			if(T[ir2]>0)
				Tguess=T[ir2];
			else
				Tguess=150.0;
			ur[ir2]=Pr[ir2]/E[ir2];
			ur[ir2]=ur[ir2]/sqrt(1.0-ur[ir2]*ur[ir2]);
			if(N[ir2]>1){
				Etot=sqrt(E[ir2]*E[ir2]-Pr[ir2]*Pr[ir2]);
				u0=sqrt(1.0+ur[ir2]*ur[ir2]);
				volume=u0*PI*delr2*2.0*b3d->ETAMAX*tau;
				etarget=Etot/(netevents*volume);
				rhotarget=N[ir2]/(netevents*volume);
				etarget=etarget/rhotarget;
				mass=M[ir2]/N[ir2];
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
				T[ir2]=Tguess;
				mu[ir2]=alpha*T[ir2];
			}
		}
	}
}

void CLocalInfo::MuCalc(){
	double r2,r,pr,e,eperp;
	double T,chi,alpha,s;
	double tau=b3d->tau;
	double etarget,starget,rhotarget,dalpha,dT,dchi,echeck;
	double ddens,accuracy,denom,dele,dels,delrho,de,x;
	double r2max=CLocalInfo::DELR2*CLocalInfo::NRBINS;
	CLocalSpeciesInfo *lsi;
	int ir2;
	CResInfo *resinfo;
	CPart *part;
	CPartMap::iterator ppos;
	multimap<int,CLocalSpeciesInfo *>::iterator lsipos;
	
	CPartMap *pmap=&(b3d->PartMap);
	ppos=pmap->begin();
	while(ppos!=pmap->end()){
		part=ppos->second;
		if(part->active){
			part->Propagate(tau);
			r2=part->r[1]*part->r[1]+part->r[2]*part->r[2];
			ir2=lrint(floor(NRBINS*r2/r2max));
			if(ir2<NRBINS){
				r=sqrt(r2);
				pr=(part->r[1]*part->p[1]+part->r[2]*part->p[2])/r;
				eperp=sqrt(part->msquared+part->p[1]*part->p[1]+part->p[2]*part->p[2]);
				e=eperp*cosh(part->y-part->eta);
				resinfo=part->resinfo;
				lsipos=lsimap.find(resinfo->code);
				if(lsipos!=lsimap.end()){
					lsi=lsipos->second;
					lsi->N[ir2]+=1;
					lsi->E[ir2]+=e;
					lsi->Pr[ir2]+=pr;
					lsi->M[ir2]+=part->GetMass();
				}
			}
		}
		ppos++;
	}
	
	pion->MuCalc();
	kaon->MuCalc();
	nucleon->MuCalc();
	delta->MuCalc();
	lambda->MuCalc();
	sigma->MuCalc();
	sigmastar->MuCalc();
	xsi->MuCalc();
	xsistar->MuCalc();
	omega->MuCalc();

	if(printing){
		pion->MuPrint();
		kaon->MuPrint();
		nucleon->MuPrint();
		delta->MuPrint();
		lambda->MuPrint();
		sigma->MuPrint();
		sigmastar->MuPrint();
		xsi->MuPrint();
		xsistar->MuPrint();
		omega->MuPrint();
	}
}

void CLocalSpeciesInfo::MuPrint(){
	double rho;
	FILE *fptr;
	char command[120];
	char filename[120];
	sprintf(command,"mkdir -p mucalc_results/%s",b3d->run_name.c_str());
	system(command);
	sprintf(filename,"mucalc_results/%s/%s_tau%g.dat",b3d->run_name.c_str(),name.c_str(),b3d->tau);
	fptr=fopen(filename,"w");
	fprintf(fptr,"#------- tau=%6.3f --------\n",b3d->tau);
	fprintf(fptr,"#r(fm)     T      mu    u_r      rho       N\n");
	for(int ir2=0;ir2<CLocalInfo::NRBINS;ir2++){
		rho=N[ir2]/(PI*CLocalInfo::DELR2*2.0*b3d->ETAMAX*b3d->tau*CLocalInfo::NETEVENTS);
		fprintf(fptr,"%6.3f %7.2f %7.2f %7.2f %9.7f %7d\n",
		sqrt(CLocalInfo::DELR2*(0.5+ir2)),T[ir2],mu[ir2],ur[ir2],rho,N[ir2]);
	}
	fclose(fptr);
}

double CLocalInfo::GetRegenFactor(CPart *part1,CPart *part2){
	int ir2,ib,NK;
	CLSIMap::iterator lsipos;
	CLocalSpeciesInfo *lsi,*lsi1,*lsi2;
	CResInfo *resinfo,*resinfo1,*resinfo2;
	double regenfactor,rsquared,E1,E2,Epi,Tpi,T1,T2,TK,muB1,muB2,mupi,muK,urpi,urK,ur1,ur2,u0;
	double Pr1,Pr2,Prpi;
	rsquared=0.25*(pow(part1->r[1]+part2->r[1],2)+pow(part2->r[2]+part2->r[2],2));
	ir2=lrint(floor(rsquared/DELR2));
	regenfactor=0.0;
	if(ir2<NRBINS){
		for(ib=1;ib<=2;ib++){
			if(ib==1)
				resinfo=part1->resinfo;
			else
				resinfo=part2->resinfo;
			lsipos=lsimap.find(resinfo->code);
			if(lsipos==lsimap.end()){
				if(resinfo->strange==0){
					if(resinfo->spin<1.0)
						lsi=nucleon;
					else
						lsi=delta;
				}
				else if(abs(resinfo->strange==1)){
					if(resinfo->spin<1.0){
						if(abs(resinfo->code)==3122)
							lsi=lambda;
						else
							lsi=sigma;
					}
					else{
						lsi=sigmastar;
					}
				}
				else if(abs(resinfo->strange==2)){
					if(resinfo->spin<1.0)
						lsi=xsi;
					else
						lsi=xsistar;
				}
				else
					lsi=omega;
			}
			else
				lsi=lsipos->second;
			
			if(ib==1){
				resinfo1=resinfo;
				lsi1=lsi;
			}
			if(ib==2){
				lsi2=lsi;
				resinfo2=resinfo;
			}
		}
		
		regenfactor=0.0;
		lsi=pion;
		if(lsi->N[ir2]>10 && lsi1->N[ir2]>10 && lsi2->N[ir2]>10){

			muB1=lsi1->mu[ir2];
			T1=lsi1->T[ir2];
			ur1=lsi1->ur[ir2];
			u0=sqrt(1.0+ur1*ur1);
			Pr1=part1->p[1]*part1->r[1]+part1->p[2]*part1->r[2];
			Pr1=Pr1/sqrt(part1->r[1]*part1->r[1]+part1->r[2]*part1->r[2]);
			E1=part1->p[0];
			E1=E1*u0-Pr1*ur1;
			
			muB2=lsi2->mu[ir2];
			T2=lsi2->T[ir2];
			ur2=lsi2->ur[ir2];
			u0=sqrt(1.0+ur2*ur2);
			Pr2=part2->p[1]*part2->r[1]+part2->p[2]*part2->r[2];
			Pr2=Pr2/sqrt(part2->r[1]*part2->r[1]+part2->r[2]*part2->r[2]);
			E2=part2->p[0];
			E2=E2*u0-Pr2*ur2;
			
			mupi=lsi->mu[ir2];
			Tpi=lsi->T[ir2];
			urpi=lsi->ur[ir2];
			u0=sqrt(1.0+urpi*urpi);
			Prpi=Pr1+Pr2;
			Epi=part1->p[0]+part2->p[0];
			Epi=Epi*u0-Prpi*urpi;
	
			NK=abs(resinfo1->strange+resinfo2->strange);
			if(NK==0){
				regenfactor=exp((5*mupi/Tpi)-(muB1/T1)-(muB2/T2)-Epi/Tpi+E1/T1+E2/T2);
				//regenfactor=exp((5*mupi/Tpi)-(muB1/T1)-(muB2/T2));
			}
			else{
				lsi=kaon;
				if(lsi->N[ir2]>10){
					TK=lsi->T[ir2];
					urK=lsi->ur[ir2];
					muK=lsi->mu[ir2];
					regenfactor=exp(((5-NK)*mupi/Tpi)+(NK*muK/TK)-(muB1/T1)-(muB2/T2)-Epi/Tpi+E1/T1+E2/T2);
					//regenfactor=exp(((5-NK)*mupi/Tpi)+(NK*muK/TK)-(muB1/T1)-(muB2/T2));
				}				
			}
			if(regenfactor>1.0){
				printf(">>>>> regenfactor=%g\n",regenfactor);
				printf("tau=%g, r=%g, PIDs=%d,%d\n",part1->tau0,sqrt((ir2+0.5)*CLocalInfo::DELR2),resinfo1->code,resinfo2->code);
				printf("mupi=%g, muK=%g, muB1=%g, muB2=%g, Tpi=%g, T1=%g, T2=%g, E1=%g, E2=%g,Epi=%g\n",mupi,NK*muK,muB1,muB2,Tpi,T1,T2,E1,E2,Epi);
				printf("ur1=%g, ur2=%g, urpi=%g, urK=%g\n",ur1,ur2,urpi,NK*urK);
				regenfactor=1.0;
			}
		}
	}
	return regenfactor;
}
