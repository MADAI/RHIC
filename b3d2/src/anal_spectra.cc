#ifndef __ANAL_SPECTRA_CC__
#define __ANAL_SPECTRA_CC__
#include "b3d.h"
//
double CB3D::CalcSpectra_PHENIX(){
	string anal_output_filename;
	FILE *anal_output;
	string detailsfilename,dirname,detailsdirname,command;
	bool WRITEDETAILS=parameter::getB(parmap,"WRITEDETAILS",false);
	bool WRITE_SPECTRA=false;
	FILE *detailsfile;
	const int NSPECIES=2,NPHI=12;
	double PHENIX_PTMAX[NSPECIES]={1200,1600},PHENIX_PTMIN[NSPECIES]={200,400};
	const int PHENIX_NPTBINS=20;
	double PHENIX_DELPT=100.0;
	double binnedPHENIX_spectra[NSPECIES][PHENIX_NPTBINS]={0.0};
	double Rx[PHENIX_NPTBINS][NPHI]={0.0},PHENIX_Ry[PHENIX_NPTBINS][NPHI]={0.0},PHENIX_Rz2[PHENIX_NPTBINS][NPHI]={0.0};
	double Rx2[PHENIX_NPTBINS][NPHI]={0.0},Ry2[PHENIX_NPTBINS][NPHI]={0.0},Rz2[PHENIX_NPTBINS][NPHI]={0.0};
	double wR_tot_PHENIX[PHENIX_NPTBINS][NPHI]={0.0},wR[PHENIX_NPTBINS][NPHI]={0.0};
	double Rout_PHENIX[PHENIX_NPTBINS][NPHI],Rlong_PHENIX[PHENIX_NPTBINS][NPHI],Rside_PHENIX[PHENIX_NPTBINS][NPHI];
	
	int ievent=1,nparts,nevents,ispecies,ID,ipt,iphi,n;
	double yield,pt,meanpt,spread,sigma,sigmapt,sigmaspread,etot=0.0,mass,pz,rapidity,et,phi;
	double DELRAPIDITY,pionspectra;
	double dmult,spectranorm,xsum,degen,phenixpionyield,weight;
	bool reality;
	int neventsmax=parameter::getI(parmap,"B3D_NEVENTSMAX",10);
	double rout,rlong,rside,tau;
	double w,wtot;
	double alpha,boseweight;
	CPartMap::iterator ppos;
	CPart *part;
	dirname="analysis/"+run_name+"/"+qualifier;
	command="mkdir -p "+dirname;
	system(command.c_str());
	anal_output_filename=dirname+"/results_phenix_spectra.dat";
	anal_output=fopen(anal_output_filename.c_str(),"w");

	parameter::set(parmap,string("B3D_RESONANCES_DECAYS_FILE"),string("progdata/madai/resinfo/decays_pdg_weak.dat"));
	parameter::set(parmap,string("B3D_RESONANCES_INFO_FILE"),string("progdata/madai/resinfo/resonances_pdg_weak.dat"));
	reslist=new CResList();
	reslist->ReadResInfo();
	DELRAPIDITY=2.0*parameter::getD(parmap,"B3D_ETAMAX",1.0);
	COLLISIONS=false;
	if(WRITEDETAILS){
		detailsdirname="analysis/"+run_name+"/"+qualifier+"/details";
		command="mkdir -p "+detailsdirname;
		system(command.c_str());
	}
	double meanpt_pion=0.0,meanpt_kaon=0.0,meanpt_proton=0.0,meanpt_omega=0.0;
	long long int npions=0,nkaons=0,nprotons=0,nomegas=0;
	nevents=0;
	do{
		KillAllParts();
		nparts=ReadOSCAR(nevents+1);
		if(nparts>0){
			nevents+=1;
			//printf("before, nparts=%d =? %d\n",nparts,int(FinalPartMap.size()));
			PerformAllActions(); // Decays unstable particles
			//printf("after decays, nparts=%d\n",int(FinalPartMap.size()));
			for(ppos=FinalPartMap.begin();ppos!=FinalPartMap.end();ppos++){
				part=ppos->second;
				pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
				et=sqrt(part->msquared+pt*pt);
				pz=et*sinh(part->y-part->eta);
				phi=acos(fabs(part->p[1]/pt));
				etot+=sqrt(et*et+part->p[3]*part->p[3]);
				ID=part->resinfo->code;
				if(ID!=111 && abs(ID)!=211 && abs(ID)!=2212 && abs(ID)!=2112 && ID!=22 && abs(ID)!=321 && abs(ID)!=311){
					printf("funny ID=%d\n",ID);
					exit(1);
				}
				weight=part->weight;
				reality=part->reality;
				ispecies=-1;
				if(ID==211 || ID==-211) ispecies=0;
				else if(ID==321 || ID==-321) ispecies=1;
				else if(ID==2112 || ID==2212 || ID==-2112 || ID==-2212) ispecies=2;
				else if(ID==3334 || ID==-3334) ispecies=3;
				//			else if(ID!=22 && ID!=-311){
				//				printf("why is particle here? ID=%d\n",ID);
				//			}
				//else if(abs(ID)!=311 &ma& abs(ID)!=3122 && abs(ID)!=22 && abs(ID)!=3222 
				// && abs(ID)!=3312 && abs(ID)!=3112 && abs(ID)!=3322 && abs(ID)!=3334) 
				//printf("ID=%d\n",ID);
				//else if(ID==3334 || ID==-3334) ispecies=3;
			
				if(ispecies==0){
					if(pt>PHENIX_PTMIN[0] && pt<PHENIX_PTMAX[0]){
						npions+=weight;
						meanpt_pion+=pt*weight;
					}
				}
				if(ispecies==1) {
					if(pt>PHENIX_PTMIN[1] && pt<PHENIX_PTMAX[1]){
						nkaons+=weight;
						meanpt_kaon+=pt*weight;
					}
				}
				if(ispecies==2){
						nprotons+=weight;
						meanpt_proton+=pt*weight;	
				}
				if(ispecies==3){
					nomegas+=weight;
					meanpt_omega+=pt*weight;
				}
			
				if(ispecies==0 || ispecies==1){
					ipt=lrint(floor(pt/PHENIX_DELPT));
					iphi=lrint(floor(2.0*NPHI*phi/PI));
					if(iphi>=NPHI || iphi<0){
						printf("iphi=%d???\n",iphi);
						exit(1);
					}
					if(ipt<PHENIX_NPTBINS){
						binnedPHENIX_spectra[ispecies][ipt]+=weight;
					}
					if(ispecies==0 && ipt<PHENIX_NPTBINS){
						part->GetHBTPars(tau,rout,rside,rlong);
						//printf("tau=%g, rout=%g, rside=%g, rlong=%g\n",tau,rout,rside,rlong);
						//					printf("rout=%g, rside=%g, rlong=%g, tau=%g\n",rout,rside,rlong,tau);
						tau=part->tau0;
						wR_tot_PHENIX[ipt][iphi]+=weight;
						w=exp(-pow(tau-20.0,2)/800.0);
						wR[ipt][iphi]+=w*weight;
						Rx2[ipt][iphi]+=w*rout*rout*weight;
						Ry2[ipt][iphi]+=w*rside*rside*weight;
						Rz2[ipt][iphi]+=w*rlong*rlong*weight;
						Rx[ipt][iphi]+=w*rout*weight;
					}
				}
			}
		}
	}while(!feof(oscarfile) && nevents<neventsmax);
	if(nevents!=neventsmax){
		printf("EVENT SHORTAGE?\n");
	}
	printf("nevents=%d, neventsmax=%d\n",nevents,neventsmax);
	printf("<pt>=%g for pions and %g for kaons\n",meanpt_pion/double(npions),meanpt_kaon/double(nkaons));
	printf("npions=%lld, nkaons=%lld, nprotons=%lld\n",npions,nkaons,nprotons);
	printf("read in info for %d events\n",nevents);
	/*
	for(ipt=0;ipt<PHENIX_NPTBINS;ipt++){
	printf("%6.2f ",(ipt+0.5)*50.0);
	for(ispecies=0;ispecies<3;ispecies++){
	printf("%15.6e ",binnedPHENIX_spectra[ispecies][ipt]);
	}
	printf("\n");
	}
	*/
	fclose(oscarfile);
	oscarfile=NULL;
	
	for(ispecies=0;ispecies<NSPECIES;ispecies++){
		if(ispecies==0)
			degen=2;
		else if(ispecies==1)
			degen=2;
		else if(ispecies==2)
			degen=4;
		else if (ispecies==3)
			degen=2;
		for(ipt=0;ipt<PHENIX_NPTBINS;ipt++){
			pt=(0.5+ipt)*PHENIX_DELPT;
			if(ispecies==0){
				binnedPHENIX_spectra[0][ipt]=0.0;
				pionspectra=0.0;
				for(iphi=0;iphi<NPHI;iphi++){
					Rx2[ipt][iphi]=Rx2[ipt][iphi]/wR[ipt][iphi];
					Ry2[ipt][iphi]=Ry2[ipt][iphi]/wR[ipt][iphi];
					Rz2[ipt][iphi]=Rz2[ipt][iphi]/wR[ipt][iphi];
					Rx[ipt][iphi]=Rx[ipt][iphi]/wR[ipt][iphi];
					//printf("ipt=%d, iphi=%d : Rx2/Ry2/Rz2=%g,%g,%g,wR=%g\n",ipt,iphi,Rx2[ipt][iphi],Ry2[ipt][iphi],Rz2[ipt][iphi],wR[ipt][iphi]);
					Rout_PHENIX[ipt][iphi]=sqrt(Rx2[ipt][iphi]-Rx[ipt][iphi]*Rx[ipt][iphi]);
					Rside_PHENIX[ipt][iphi]=sqrt(Ry2[ipt][iphi]);
					Rlong_PHENIX[ipt][iphi]=sqrt(Rz2[ipt][iphi]);
					//printf("ipt=%d, iphi=%d : Rout/side/long=%g,%g,%g\n",ipt,iphi,Rout_PHENIX[ipt][iphi],Rside_PHENIX[ipt][iphi],Rlong_PHENIX[ipt][iphi]);
					alpha=pow(2.0*PI,1.5)*HBARC*HBARC*HBARC/(Rout_PHENIX[ipt][iphi]*Rlong_PHENIX[ipt][iphi]*Rside_PHENIX[ipt][iphi]);
					alpha*=wR[ipt][iphi]/(sqrt(pt*pt+139.57*139.57));
					alpha*=double(NPHI)/(2.0*PI*PHENIX_DELPT*pt*DELRAPIDITY*degen*double(nevents*NSAMPLE));
					boseweight=0.0;
					for(n=1;n<10;n+=1){
						boseweight+=pow(double(n),-1.5)*pow(alpha,n-1);
					}
					boseweight=(1.0-wR[ipt][iphi]/wR_tot_PHENIX[ipt][iphi])+wR[ipt][iphi]*boseweight/wR_tot_PHENIX[ipt][iphi];
					//printf("ipt=%d, iphi=%d, boseweight=%g\n",ipt,iphi,boseweight);
					//boseweight=1.0;
					pionspectra=(wR[ipt][iphi]*boseweight+(wR_tot_PHENIX[ipt][iphi]-wR[ipt][iphi]));
					binnedPHENIX_spectra[0][ipt]+=pionspectra;
				}
				//printf("pt=%g, alpha=%g, spectra0=%g, spectra=%g, boseweight=%g\n",pt,alpha,
				//binnedPHENIX_spectra[0][ipt],pionspectra,pionspectra/binnedPHENIX_spectra[0][ipt]);
			}
		}
	
		yield=meanpt=0.0;
		if(WRITEDETAILS){
			if(ispecies==0) detailsfilename=detailsdirname+"/phenix_spectra_pion.dat";
			if(ispecies==1) detailsfilename=detailsdirname+"/phenix_spectra_kaon.dat";
			//printf("detailsfilename=%s\n",detailsfilename.c_str());
			detailsfile=fopen(detailsfilename.c_str(),"w");
			if(ispecies==0) fprintf(detailsfile,"# PION SPECTRA \n");
			if(ispecies==1) fprintf(detailsfile,"# KAON SPECTRA \n");
		}
		for(ipt=0;ipt<PHENIX_NPTBINS;ipt++){
			pt=PHENIX_DELPT*(0.5+ipt);
			dmult=binnedPHENIX_spectra[ispecies][ipt]/double(nevents);
			if(pt>PHENIX_PTMIN[ispecies] && pt<PHENIX_PTMAX[ispecies]){
				yield+=dmult;
				meanpt+=dmult*pt;
			}
			sigma=0.1*dmult*pt/1000.0;
			if(WRITE_SPECTRA){
				if(ispecies==0) fprintf(anal_output,"PHENIX_SPECTRA_PION_pt%d %g %g\n",int(lrint(pt)),dmult,sigma);
				if(ispecies==1) fprintf(anal_output,"PHENIX_SPECTRA_KAON_pt%d %g %g\n",int(lrint(pt)),dmult,sigma);
			}
			if(WRITEDETAILS && (ispecies==0 || ispecies==1)){
				fprintf(detailsfile,"%g %g\n",pt/1000.0,
				binnedPHENIX_spectra[ispecies][ipt]/(degen*PHENIX_DELPT*DELRAPIDITY*2*PI*pt/(1.E6)*double(nevents)));
			}	
		}
		meanpt=meanpt/yield;
		yield=yield/(DELRAPIDITY);
		//printf("yield=%g, meanpt=%g\n",yield,meanpt);
		
		if(ispecies==0){
			fprintf(anal_output,"PHENIX_SPECTRA_PION_YIELD %g %g\n",yield,0.06*yield);	
			fprintf(anal_output,"PHENIX_SPECTRA_PION_MEANPT %g %g\n",meanpt,0.06*meanpt);
		}
		if(ispecies==1){
			fprintf(anal_output,"PHENIX_SPECTRA_KAON_YIELD %g %g\n",yield,0.06*yield);	
			fprintf(anal_output,"PHENIX_SPECTRA_KAON_MEANPT %g %g\n",meanpt,0.06*meanpt);
		}
		if(WRITEDETAILS){
			fflush(detailsfile);
			fclose(detailsfile);
		}
	}
	delete reslist;
	fclose(anal_output);
	return phenixpionyield;
}

void CB3D::CalcSpectra_PHENIXppbar(){
	string anal_output_filename;
	FILE *anal_output;
	string detailsfilename,dirname,detailsdirname,command;
	bool WRITEDETAILS=parameter::getB(parmap,"WRITEDETAILS",false);
	bool WRITE_SPECTRA=false;
	FILE *detailsfile;
	CResInfo *resinfoptr;
	const int NPHI=12;
	double PHENIX_PTMAX=2000,PHENIX_PTMIN=600;
	const int PHENIX_NPTBINS=20;
	double PHENIX_DELPT=100.0;
	double binnedPHENIX_spectra[PHENIX_NPTBINS]={0.0};
	int neventsmax=parameter::getI(parmap,"B3D_NEVENTSMAX",10);
	int ievent=1,nparts,nevents,ispecies,ID,ipt,iphi,n;
	double yield,pt,meanpt,spread,sigma,sigmapt,sigmaspread,etot=0.0,mass,pz,rapidity,et,phi;
	double DELRAPIDITY,pionspectra;
	double dmult,spectranorm,xsum,degen,phenixpionyield,weight;
	bool reality;
	
	double rout,rlong,rside,tau;
	double w,wtot;
	double alpha,boseweight;
	CPartMap::iterator ppos;
	CPart *part;
	dirname="analysis/"+run_name+"/"+qualifier;
	command="mkdir -p "+dirname;
	system(command.c_str());
	anal_output_filename=dirname+"/results_phenix_ppbarspectra.dat";
	anal_output=fopen(anal_output_filename.c_str(),"w");
	
	parameter::set(parmap,string("B3D_RESONANCES_DECAYS_FILE"),string("progdata/madai/resinfo/decays_pdg_weak.dat"));
	parameter::set(parmap,string("B3D_RESONANCES_INFO_FILE"),string("progdata/madai/resinfo/resonances_pdg_weak.dat"));
	reslist=new CResList();
	reslist->ReadResInfo();
	resinfoptr=reslist->GetResInfoPtr(3122);
	resinfoptr->decay=0;
	resinfoptr=reslist->GetResInfoPtr(-3122);
	resinfoptr->decay=0;
	
	DELRAPIDITY=2.0*parameter::getD(parmap,"B3D_ETAMAX",1.0);
	COLLISIONS=false;
	if(WRITEDETAILS){
		detailsdirname="analysis/"+run_name+"/"+qualifier+"/details";
		command="mkdir -p "+detailsdirname;
		system(command.c_str());
	}
	double meanpt_proton=0.0;
	long long int nprotons=0;
	

	NSAMPLE=parameter::getI(parmap,"B3D_NSAMPLE",1);
	nevents=0;
	do{
		KillAllParts();
		nparts=ReadOSCAR(nevents+1);
		if(nparts>0){
			nevents+=1;
			//printf("before, nparts=%d =? %d\n",nparts,int(FinalPartMap.size()));
			
			PerformAllActions(); // Decays unstable particles
			//printf("after decays, nparts=%d\n",int(FinalPartMap.size()));			
			
			for(ppos=FinalPartMap.begin();ppos!=FinalPartMap.end();ppos++){
				part=ppos->second;
				pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
				et=sqrt(part->msquared+pt*pt);
				pz=et*sinh(part->y-part->eta);
				phi=acos(fabs(part->p[1]/pt));
				etot+=sqrt(et*et+part->p[3]*part->p[3]);
				ID=part->resinfo->code;
				weight=part->weight;
				reality=part->reality;
				ispecies=-1;
				if(ID==2112 || ID==2212 || ID==-2112 || ID==-2212) ispecies=2;
				if(ispecies==2){
					if(pt>PHENIX_PTMIN && pt<PHENIX_PTMAX){
						nprotons+=weight;
						meanpt_proton+=pt*weight;	
					}
					ipt=lrint(floor(pt/PHENIX_DELPT));
					if(ipt<PHENIX_NPTBINS){
						binnedPHENIX_spectra[ipt]+=weight;
					}
				}
			}
		}
	}while(!feof(oscarfile) && nevents<neventsmax);
	if(nevents!=neventsmax){
		printf("EVENT SHORTAGE?\n");
	}
	fclose(oscarfile);
	oscarfile=NULL;

	degen=4.0;
	spectranorm=0.0;
	yield=0.0;
	meanpt=0.0;
	for(ipt=0;ipt<PHENIX_NPTBINS;ipt++){
		pt=PHENIX_DELPT*(0.5+ipt);
		dmult=binnedPHENIX_spectra[ipt]/double(nevents);
		if(pt>PHENIX_PTMIN && pt<PHENIX_PTMAX){
			yield+=dmult;
			meanpt+=dmult*pt;
		}
		sigma=0.1*dmult*pt/1000.0;
		if(WRITE_SPECTRA)
			fprintf(anal_output,"PHENIX_SPECTRA_PPBAR_pt%d %g %g\n",int(lrint(pt)),dmult,sigma);
	}
	meanpt=meanpt/yield;
	printf("meanpt from array=%g\n",meanpt);
	yield=0.5*yield/DELRAPIDITY ; // dN/dy of p+pbar in acceptance, factor of 0.5 for neutrons
	fprintf(anal_output,"PHENIX_SPECTRA_PPBAR_YIELD %g %g\n",yield,0.06*yield);
	fprintf(anal_output,"PHENIX_SPECTRA_PPBAR_MEANPT %g %g\n",meanpt,0.06*meanpt);

	if(WRITEDETAILS){
		detailsfilename=detailsdirname+"/phenix_spectra_ppbar.dat";
		detailsfile=fopen(detailsfilename.c_str(),"w");
		fprintf(detailsfile,"# PPBAR SPECTRA \n");
		for(ipt=0;ipt<PHENIX_NPTBINS;ipt++){
			pt=PHENIX_DELPT*(0.5+ipt);
			fprintf(detailsfile,"%g %g\n",pt/1000.0,
			binnedPHENIX_spectra[ipt]/(degen*PHENIX_DELPT*DELRAPIDITY*2*PI*pt/(1.E6)*double(nevents)));
		}
		fflush(detailsfile);
		fclose(detailsfile);
	}
	delete reslist;
	fclose(anal_output);
}

double CB3D::CalcSpectra_ALICE(){
	string anal_output_filename;
	FILE *anal_output;
	string detailsfilename,dirname,detailsdirname,command;
	bool WRITEDETAILS=parameter::getB(parmap,"WRITEDETAILS",false);
	bool WRITE_SPECTRA=false;
	FILE *detailsfile;
	const int NSPECIES=3,NPHI=12;
	double ALICE_PTMAX[NSPECIES]={1200,1600,2000},ALICE_PTMIN[NSPECIES]={100,200,300};
	const int ALICE_NPTBINS=20;
	double ALICE_DELPT=100.0;
	double binnedALICE_spectra[NSPECIES][ALICE_NPTBINS]={0.0};
	double Rx[ALICE_NPTBINS][NPHI]={0.0},ALICE_Ry[ALICE_NPTBINS][NPHI]={0.0},ALICE_Rz2[ALICE_NPTBINS][NPHI]={0.0};
	double Rx2[ALICE_NPTBINS][NPHI]={0.0},Ry2[ALICE_NPTBINS][NPHI]={0.0},Rz2[ALICE_NPTBINS][NPHI]={0.0};
	double wR_tot_ALICE[ALICE_NPTBINS][NPHI]={0.0},wR[ALICE_NPTBINS][NPHI]={0.0};
	double Rout_ALICE[ALICE_NPTBINS][NPHI],Rlong_ALICE[ALICE_NPTBINS][NPHI],Rside_ALICE[ALICE_NPTBINS][NPHI];
	
	int ievent=1,nparts,nevents,ispecies,ID,ipt,iphi,n;
	double yield,pt,meanpt,spread,sigma,sigmapt,sigmaspread,etot=0.0,mass,pz,rapidity,et,phi;
	double DELRAPIDITY,pionspectra;
	double dmult,spectranorm,xsum,degen,ALICEpionyield,weight;
	bool reality;
	int neventsmax=parameter::getI(parmap,"B3D_NEVENTSMAX",10);
	double rout,rlong,rside,tau;
	double w,wtot;
	double alpha,boseweight;
	CPartMap::iterator ppos;
	CPart *part;
	dirname="analysis/"+run_name+"/"+qualifier;
	command="mkdir -p "+dirname;
	system(command.c_str());
	anal_output_filename=dirname+"/results_alice_spectra.dat";
	anal_output=fopen(anal_output_filename.c_str(),"w");

	parameter::set(parmap,string("B3D_RESONANCES_DECAYS_FILE"),string("progdata/madai/resinfo/decays_pdg_weak.dat"));
	parameter::set(parmap,string("B3D_RESONANCES_INFO_FILE"),string("progdata/madai/resinfo/resonances_pdg_weak.dat"));
	reslist=new CResList();
	reslist->ReadResInfo();
	DELRAPIDITY=2.0*parameter::getD(parmap,"B3D_ETAMAX",1.0);
	COLLISIONS=false;
	if(WRITEDETAILS){
		detailsdirname="analysis/"+run_name+"/"+qualifier+"/details";
		command="mkdir -p "+detailsdirname;
		system(command.c_str());
	}
	double meanpt_pion=0.0,meanpt_kaon=0.0,meanpt_proton=0.0,meanpt_omega=0.0;
	long long int npions=0,nkaons=0,nprotons=0,nomegas=0;
	nevents=0;
	do{
		KillAllParts();
		nparts=ReadOSCAR(nevents+1);
		if(nparts>0){
			nevents+=1;
			PerformAllActions(); // Decays unstable particles
			for(ppos=FinalPartMap.begin();ppos!=FinalPartMap.end();ppos++){
				part=ppos->second;
				pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
				et=sqrt(part->msquared+pt*pt);
				pz=et*sinh(part->y-part->eta);
				phi=acos(fabs(part->p[1]/pt));
				etot+=sqrt(et*et+part->p[3]*part->p[3]);
				ID=part->resinfo->code;
				if(ID!=111 && abs(ID)!=211 && abs(ID)!=2212 && abs(ID)!=2112 && ID!=22 && abs(ID)!=321 && abs(ID)!=311){
					printf("funny ID=%d\n",ID);
					exit(1);
				}
				weight=part->weight;
				reality=part->reality;
				ispecies=-1;
				if(ID==211 || ID==-211) ispecies=0;
				else if(ID==321 || ID==-321) ispecies=1;
				else if(ID==2112 || ID==2212 || ID==-2112 || ID==-2212) ispecies=2;
				else if(ID==3334 || ID==-3334) ispecies=3;
				//			else if(ID!=22 && ID!=-311){
				//				printf("why is particle here? ID=%d\n",ID);
				//			}
				//else if(abs(ID)!=311 &ma& abs(ID)!=3122 && abs(ID)!=22 && abs(ID)!=3222 
				// && abs(ID)!=3312 && abs(ID)!=3112 && abs(ID)!=3322 && abs(ID)!=3334) 
				//printf("ID=%d\n",ID);
				//else if(ID==3334 || ID==-3334) ispecies=3;
			
				if(ispecies==0){
					if(pt>ALICE_PTMIN[0] && pt<ALICE_PTMAX[0]){
						npions+=weight;
						meanpt_pion+=pt*weight;
					}
				}
				if(ispecies==1) {
					if(pt>ALICE_PTMIN[1] && pt<ALICE_PTMAX[1]){
						nkaons+=weight;
						meanpt_kaon+=pt*weight;
					}
				}
				if(ispecies==2){
					if(pt>ALICE_PTMIN[2] && pt<ALICE_PTMAX[2]){
						nprotons+=weight;
						meanpt_proton+=pt*weight;
					}
				}
				if(ispecies==3){
					nomegas+=weight;
					meanpt_omega+=pt*weight;
				}
			
				if(ispecies==0 || ispecies==1 || ispecies==2){
					ipt=lrint(floor(pt/ALICE_DELPT));
					iphi=lrint(floor(2.0*NPHI*phi/PI));
					if(iphi>=NPHI || iphi<0){
						printf("iphi=%d???\n",iphi);
						exit(1);
					}
					if(ipt<ALICE_NPTBINS){
						binnedALICE_spectra[ispecies][ipt]+=weight;
					}
					if(ispecies==0 && ipt<ALICE_NPTBINS){
						part->GetHBTPars(tau,rout,rside,rlong);
						//printf("tau=%g, rout=%g, rside=%g, rlong=%g\n",tau,rout,rside,rlong);
						//					printf("rout=%g, rside=%g, rlong=%g, tau=%g\n",rout,rside,rlong,tau);
						tau=part->tau0;
						wR_tot_ALICE[ipt][iphi]+=weight;
						w=exp(-pow(tau-20.0,2)/800.0);
						wR[ipt][iphi]+=w*weight;
						Rx2[ipt][iphi]+=w*rout*rout*weight;
						Ry2[ipt][iphi]+=w*rside*rside*weight;
						Rz2[ipt][iphi]+=w*rlong*rlong*weight;
						Rx[ipt][iphi]+=w*rout*weight;
					}
				}
			}
		}
	}while(!feof(oscarfile) && nevents<neventsmax);
	if(nevents!=neventsmax){
		printf("EVENT SHORTAGE?\n");
	}
	printf("nevents=%d, neventsmax=%d\n",nevents,neventsmax);
	printf("<pt>=%g for pions, %g for kaons and %g for protons\n",
	meanpt_pion/double(npions),meanpt_kaon/double(nkaons),meanpt_proton/double(nprotons));
	printf("npions=%lld, nkaons=%lld, nprotons=%lld\n",npions,nkaons,nprotons);
	printf("read in info for %d events\n",nevents);
	fclose(oscarfile);
	oscarfile=NULL;
	
	for(ispecies=0;ispecies<NSPECIES;ispecies++){
		if(ispecies==0)
			degen=2;
		else if(ispecies==1)
			degen=2;
		else if(ispecies==2)
			degen=4;
		else if (ispecies==3)
			degen=2;
		for(ipt=0;ipt<ALICE_NPTBINS;ipt++){
			pt=(0.5+ipt)*ALICE_DELPT;
			if(ispecies==0){
				binnedALICE_spectra[0][ipt]=0.0;
				pionspectra=0.0;
				for(iphi=0;iphi<NPHI;iphi++){
					Rx2[ipt][iphi]=Rx2[ipt][iphi]/wR[ipt][iphi];
					Ry2[ipt][iphi]=Ry2[ipt][iphi]/wR[ipt][iphi];
					Rz2[ipt][iphi]=Rz2[ipt][iphi]/wR[ipt][iphi];
					Rx[ipt][iphi]=Rx[ipt][iphi]/wR[ipt][iphi];
					//printf("ipt=%d, iphi=%d : Rx2/Ry2/Rz2=%g,%g,%g,wR=%g\n",ipt,iphi,Rx2[ipt][iphi],Ry2[ipt][iphi],Rz2[ipt][iphi],wR[ipt][iphi]);
					Rout_ALICE[ipt][iphi]=sqrt(Rx2[ipt][iphi]-Rx[ipt][iphi]*Rx[ipt][iphi]);
					Rside_ALICE[ipt][iphi]=sqrt(Ry2[ipt][iphi]);
					Rlong_ALICE[ipt][iphi]=sqrt(Rz2[ipt][iphi]);
					//printf("ipt=%d, iphi=%d : Rout/side/long=%g,%g,%g\n",ipt,iphi,Rout_ALICE[ipt][iphi],Rside_ALICE[ipt][iphi],Rlong_ALICE[ipt][iphi]);
					alpha=pow(2.0*PI,1.5)*HBARC*HBARC*HBARC/(Rout_ALICE[ipt][iphi]*Rlong_ALICE[ipt][iphi]*Rside_ALICE[ipt][iphi]);
					alpha*=wR[ipt][iphi]/(sqrt(pt*pt+139.57*139.57));
					alpha*=double(NPHI)/(2.0*PI*ALICE_DELPT*pt*DELRAPIDITY*degen*double(nevents*NSAMPLE));
					boseweight=0.0;
					for(n=1;n<10;n+=1){
						boseweight+=pow(double(n),-1.5)*pow(alpha,n-1);
					}
					boseweight=(1.0-wR[ipt][iphi]/wR_tot_ALICE[ipt][iphi])+wR[ipt][iphi]*boseweight/wR_tot_ALICE[ipt][iphi];
					//printf("ipt=%d, iphi=%d, boseweight=%g\n",ipt,iphi,boseweight);
					//boseweight=1.0;
					pionspectra=(wR[ipt][iphi]*boseweight+(wR_tot_ALICE[ipt][iphi]-wR[ipt][iphi]));
					binnedALICE_spectra[0][ipt]+=pionspectra;
				}
				//printf("pt=%g, alpha=%g, spectra0=%g, spectra=%g, boseweight=%g\n",pt,alpha,
				//binnedALICE_spectra[0][ipt],pionspectra,pionspectra/binnedALICE_spectra[0][ipt]);
			}
		}
	
		yield=meanpt=0.0;
		if(WRITEDETAILS){
			if(ispecies==0) detailsfilename=detailsdirname+"/ALICE_spectra_pion.dat";
			if(ispecies==1) detailsfilename=detailsdirname+"/ALICE_spectra_kaon.dat";
			if(ispecies==2) detailsfilename=detailsdirname+"/ALICE_spectra_proton.dat";
			printf("detailsfilename=%s\n",detailsfilename.c_str());
			detailsfile=fopen(detailsfilename.c_str(),"w");
			if(ispecies==0) fprintf(detailsfile,"# PION SPECTRA \n");
			if(ispecies==1) fprintf(detailsfile,"# KAON SPECTRA \n");
			if(ispecies==2) fprintf(detailsfile,"# PPBAR SPECTRA \n");
		}
		for(ipt=0;ipt<ALICE_NPTBINS;ipt++){
			pt=ALICE_DELPT*(0.5+ipt);
			dmult=binnedALICE_spectra[ispecies][ipt]/double(nevents);
			if(pt>ALICE_PTMIN[ispecies] && pt<ALICE_PTMAX[ispecies]){
				yield+=dmult;
				meanpt+=dmult*pt;
			}
			sigma=0.1*dmult*pt/1000.0;
			if(WRITE_SPECTRA){
				if(ispecies==0) fprintf(anal_output,"ALICE_SPECTRA_PION_pt%d %g %g\n",int(lrint(pt)),dmult,sigma);
				if(ispecies==1) fprintf(anal_output,"ALICE_SPECTRA_KAON_pt%d %g %g\n",int(lrint(pt)),dmult,sigma);
				if(ispecies==2) fprintf(anal_output,"ALICE_SPECTRA_PPBAR_pt%d %g %g\n",int(lrint(pt)),dmult,sigma);
			}
			if(WRITEDETAILS && (ispecies==0 || ispecies==1 || ispecies==2)){
				fprintf(detailsfile,"%g %g\n",pt/1000.0,
				binnedALICE_spectra[ispecies][ipt]/(degen*ALICE_DELPT*DELRAPIDITY*2*PI*pt/(1.E6)*double(nevents)));
			}	
		}
		meanpt=meanpt/yield;
		yield=yield/(DELRAPIDITY);
		printf("yield=%g, meanpt=%g\n",yield,meanpt);
		
		if(ispecies==0){
			fprintf(anal_output,"ALICE_SPECTRA_PION_YIELD %g %g\n",yield,0.06*yield);	
			fprintf(anal_output,"ALICE_SPECTRA_PION_MEANPT %g %g\n",meanpt,0.06*meanpt);
		}
		if(ispecies==1){
			fprintf(anal_output,"ALICE_SPECTRA_KAON_YIELD %g %g\n",yield,0.06*yield);	
			fprintf(anal_output,"ALICE_SPECTRA_KAON_MEANPT %g %g\n",meanpt,0.06*meanpt);
		}
		if(ispecies==2){
			fprintf(anal_output,"ALICE_SPECTRA_PPBAR_YIELD %g %g\n",yield,0.06*yield);	
			fprintf(anal_output,"ALICE_SPECTRA_PPBAR_MEANPT %g %g\n",meanpt,0.06*meanpt);
		}
		if(WRITEDETAILS){
			fflush(detailsfile);
			fclose(detailsfile);
		}
	}
	delete reslist;
	fclose(anal_output);
	return ALICEpionyield;
}

#endif
