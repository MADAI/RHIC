#ifndef __ANALYZE_V2_CC__
#define __ANALYZE_V2_CC__
#include "b3d.h"

void CB3D::CalcV2_STAR(){
	string detailsfilename,detailsdirname;
	string anal_output_filename;
	FILE *anal_output; 
	bool WRITEDETAILS=parameter::getB(parmap,"WRITEDETAILS",false);
	bool WRITE_V2=true;
	FILE *detailsfile;
	const int NSPECIES=3;
	double STAR_PTMAX[NSPECIES]={1040,720,1360},STAR_PTMIN=160.0;
	string species[NSPECIES]={"PION","KAON","PROTON"};
	double starv2_ptweight[NSPECIES]={0.0},starv2_ptweightnorm[NSPECIES]={0.0};
	const int NPTBINS=25;
	double DELPT=80.0;
	double error,weight;
	double v2[NSPECIES][NPTBINS]={0.0};
	double v2norm[NSPECIES][NPTBINS]={0.0};
	int nparts,nevents,ispecies,ID,ipt,n;
	double yield,pt,meanpt,spread,sigma,sigmapt,sigmaspread,etot=0.0,mass,pz,rapidity,et,phi;
	bool reality;
	
	double rout,rlong,rside,tau;
	double w,wtot;
	double alpha,boseweight;
	CPartMap::iterator ppos;
	CPart *part;
	string dirname,command;
	dirname="analysis/"+run_name+"/"+qualifier;
	command="mkdir -p "+dirname;
	system(command.c_str());
	anal_output_filename=dirname+"/results_star_v2.dat";
	anal_output=fopen(anal_output_filename.c_str(),"w");

	parameter::set(parmap,string("B3D_RESONANCES_DECAYS_FILE"),string("progdata/madai/resinfo/decays_pdg_weak.dat"));
	parameter::set(parmap,string("B3D_RESONANCES_INFO_FILE"),string("progdata/madai/resinfo/resonances_pdg_weak.dat"));
	reslist=new CResList();
	reslist->ReadResInfo();
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
			//printf("before, nparts=%d\n",nparts);
			PerformAllActions(); // Decays unstable particles
			//printf("after decays, nparts=%d\n",int(FinalPartMap.size()));
			for(ppos=FinalPartMap.begin();ppos!=FinalPartMap.end();ppos++){
				part=ppos->second;
				pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
				ipt=int(lrint(floor(pt/DELPT)));
				et=sqrt(part->msquared+pt*pt);
				pz=et*sinh(part->y-part->eta);
				phi=acos(fabs(part->p[1]/pt));
				etot+=sqrt(et*et+part->p[3]*part->p[3]);
				ID=part->resinfo->code;
				weight=part->weight;
				reality=part->reality;
				ispecies=-1;
				if(ID==211 || ID==-211) ispecies=0;
				else if(ID==321 || ID==-321) ispecies=1;
				else if(ID==2112 || ID==2212 || ID==-2112 || ID==-2212) ispecies=2;
				else if(ID==3334 || ID==-3334) ispecies=3;
				if(ispecies<NSPECIES && ispecies>=0 && ipt<NPTBINS){
					v2[ispecies][ipt]+=weight*cos(2.0*phi);
					v2norm[ispecies][ipt]+=weight;
				}
			}
		}
	} while(!feof(oscarfile));
	printf("read info for %d events\n",nevents);
	fclose(oscarfile);
	oscarfile=NULL;

	for(ispecies=0;ispecies<NSPECIES;ispecies++){
		for(ipt=0;ipt<NPTBINS;ipt++){
			v2[ispecies][ipt]=v2[ispecies][ipt]/v2norm[ispecies][ipt];
			pt=(0.5+ipt)*DELPT;
			//printf("%5.2f %g\n",pt,v2[ispecies][ipt]);
			if(pt>STAR_PTMIN && pt<STAR_PTMAX[ispecies]){
				starv2_ptweight[ispecies]+=pt*v2[ispecies][ipt];
				starv2_ptweightnorm[ispecies]+=pt;
			}
		}
		starv2_ptweight[ispecies]=starv2_ptweight[ispecies]/starv2_ptweightnorm[ispecies];
		sigma=sqrt(0.0001+pow(0.075*starv2_ptweight[ispecies],2));
		fprintf(anal_output,"STAR_V2_%s_PTWEIGHT %g %g\n",species[ispecies].c_str(),starv2_ptweight[ispecies],sigma);
	
		if(WRITE_V2){
			for(ipt=0;ipt<NPTBINS;ipt++){
				pt=(0.5+ipt)*DELPT;
				if(pt>STAR_PTMIN && pt<STAR_PTMAX[ispecies]){
					error=100.0/sqrt(2.0*v2norm[ispecies][ipt]);
					fprintf(anal_output,"STAR_%s_V2_pt%d %g %g\n",species[ispecies].c_str(),int(lrint(pt)),100.0*v2[0][ipt],error);
				}
			}
		}
		
		if(WRITEDETAILS){
			if(ispecies==0)
				detailsfilename=detailsdirname+"/star_v22_pion.dat";
			if(ispecies==1)
				detailsfilename=detailsdirname+"/star_v22_kaon.dat";
			if(ispecies==2)
				detailsfilename=detailsdirname+"/star_v22_ppbar.dat";
			detailsfile=fopen(detailsfilename.c_str(),"w");
			fprintf(detailsfile,"# STAR_%s_V2\n",species[ispecies].c_str());
			for(ipt=0;ipt<NPTBINS;ipt++){
				pt=(0.5+ipt)*DELPT;
				fprintf(detailsfile,"%6.4f %g\n",pt/1000.0,100.0*v2[ispecies][ipt]);
			}
			fflush(detailsfile);
			fclose(detailsfile);
		}
	}
	delete reslist;
	fclose(anal_output);
}

void CB3D::CalcV2_ALICE(){
	string detailsfilename,detailsdirname;
	string anal_output_filename;
	FILE *anal_output;
	bool WRITEDETAILS=parameter::getB(parmap,"WRITEDETAILS",false);
	bool WRITE_V2=true;
	FILE *detailsfile;
	const int NSPECIES=3;
	double ALICE_PTMAX[NSPECIES]={1200,1600,2000},ALICE_PTMIN=200.0;
	string species[NSPECIES]={"PION","KAON","PROTON"};
	double alicev2_ptweight[NSPECIES]={0.0},alicev2_ptweightnorm[NSPECIES]={0.0};
	const int NPTBINS=20;
	double DELPT=100.0;
	double error,weight;
	double v2[NSPECIES][NPTBINS]={0.0};
	double v2norm[NSPECIES][NPTBINS]={0.0};
	int nparts,nevents,ispecies,ID,ipt,n;
	double yield,pt,meanpt,spread,sigma,sigmapt,sigmaspread,etot=0.0,mass,pz,rapidity,et,phi;
	bool reality;
	
	double rout,rlong,rside,tau;
	double w,wtot;
	double alpha,boseweight;
	CPartMap::iterator ppos;
	CPart *part;
	string dirname,command;
	dirname="analysis/"+run_name+"/"+qualifier;
	command="mkdir -p "+dirname;
	system(command.c_str());
	anal_output_filename=dirname+"/results_alice_v2.dat";
	anal_output=fopen(anal_output_filename.c_str(),"w");

	parameter::set(parmap,string("B3D_RESONANCES_DECAYS_FILE"),string("progdata/madai/resinfo/decays_pdg_weak.dat"));
	parameter::set(parmap,string("B3D_RESONANCES_INFO_FILE"),string("progdata/madai/resinfo/resonances_pdg_weak.dat"));
	reslist=new CResList();
	reslist->ReadResInfo();
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
			//printf("before, nparts=%d\n",nparts);
			PerformAllActions(); // Decays unstable particles
			//printf("after decays, nparts=%d\n",int(FinalPartMap.size()));
			for(ppos=FinalPartMap.begin();ppos!=FinalPartMap.end();ppos++){
				part=ppos->second;
				pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
				ipt=int(lrint(floor(pt/DELPT)));
				et=sqrt(part->msquared+pt*pt);
				pz=et*sinh(part->y-part->eta);
				phi=acos(fabs(part->p[1]/pt));
				etot+=sqrt(et*et+part->p[3]*part->p[3]);
				ID=part->resinfo->code;
				weight=part->weight;
				reality=part->reality;
				ispecies=-1;
				if(ID==211 || ID==-211) ispecies=0;
				else if(ID==321 || ID==-321) ispecies=1;
				else if(ID==2112 || ID==2212 || ID==-2112 || ID==-2212) ispecies=2;
				else if(ID==3334 || ID==-3334) ispecies=3;
				if(ispecies<NSPECIES && ispecies>=0 && ipt<NPTBINS){
					v2[ispecies][ipt]+=weight*cos(2.0*phi);
					v2norm[ispecies][ipt]+=weight;
				}
			}
		}
	} while(!feof(oscarfile));
	printf("read info for %d events\n",nevents);
	fclose(oscarfile);
	oscarfile=NULL;

	for(ispecies=0;ispecies<NSPECIES;ispecies++){
		for(ipt=0;ipt<NPTBINS;ipt++){
			v2[ispecies][ipt]=v2[ispecies][ipt]/v2norm[ispecies][ipt];
			pt=(0.5+ipt)*DELPT;
			//printf("%5.2f %g\n",pt,v2[ispecies][ipt]);
			if(pt>ALICE_PTMIN && pt<ALICE_PTMAX[ispecies]){
				alicev2_ptweight[ispecies]+=pt*v2[ispecies][ipt];
				alicev2_ptweightnorm[ispecies]+=pt;
			}
		}
		alicev2_ptweight[ispecies]=alicev2_ptweight[ispecies]/alicev2_ptweightnorm[ispecies];
		sigma=sqrt(0.0001+pow(0.075*alicev2_ptweight[ispecies],2));
		fprintf(anal_output,"ALICE_V2_%s_PTWEIGHT %g %g\n",species[ispecies].c_str(),alicev2_ptweight[ispecies],sigma);
	
		if(WRITE_V2){
			for(ipt=0;ipt<NPTBINS;ipt++){
				pt=(0.5+ipt)*DELPT;
				if(pt>ALICE_PTMIN && pt<ALICE_PTMAX[ispecies]){
					error=100.0/sqrt(2.0*v2norm[ispecies][ipt]);
					fprintf(anal_output,"ALICE_%s_V2_pt%d %g %g\n",species[ispecies].c_str(),int(lrint(pt)),100.0*v2[0][ipt],error);
				}
			}
		}
		
		if(WRITEDETAILS){
			if(ispecies==0)
				detailsfilename=detailsdirname+"/alice_v22_pion.dat";
			if(ispecies==1)
				detailsfilename=detailsdirname+"/alice_v22_kaon.dat";
			if(ispecies==2)
				detailsfilename=detailsdirname+"/alice_v22_ppbar.dat";
			detailsfile=fopen(detailsfilename.c_str(),"w");
			fprintf(detailsfile,"# ALICE_%s_V2\n",species[ispecies].c_str());
			for(ipt=0;ipt<NPTBINS;ipt++){
				pt=(0.5+ipt)*DELPT;
				fprintf(detailsfile,"%6.4f %g\n",pt/1000.0,100.0*v2[ispecies][ipt]);
			}
			fflush(detailsfile);
			fclose(detailsfile);
		}
	}
	delete reslist;
	fclose(anal_output);
}
#endif
