#ifndef __ANALYZE_HBT_CC__
#define __ANALYZE_HBT_CC__
#include "b3d.h"
#include "coral.h"


void CB3D::CalcHBT_STAR(){
	string anal_output_filename;
	FILE *anal_output;
	double Rx=6.0,Ry=6.8,Rz=7.6,lambda=0.775;
	bool WRITEDETAILS=parameter::getB(parmap,"WRITEDETAILS",false);
	//string HBTPAIRS="KK";
	string HBTPAIRS="PIPI";
	string details_filename;
	FILE* detailsfile;
	
	string dirname,command;
	dirname="analysis/"+run_name+"/"+qualifier;
	command="mkdir -p "+dirname;
	system(command.c_str());
	anal_output_filename=dirname+"/results_star_hbt.dat";
	anal_output=fopen(anal_output_filename.c_str(),"w");
	if(WRITEDETAILS){
		dirname="analysis/"+run_name+"/"+qualifier+"/details";
		command="mkdir -p "+dirname;
		system(command.c_str());
		details_filename=dirname+"/hbt_star.dat";
		detailsfile=fopen(details_filename.c_str(),"w");
		fprintf(detailsfile,"# STAR HBT RADII\n");
	}
	if(qualifier=="rhic200_cent0to5"){
		Rx=10.0; Ry=6.8; Rz=7.6;
	}
	else if(qualifier=="rhic200_cent5to10"){
		Rx=9.7; Ry=6.6; Rz=7.4;
	}
	else if(qualifier=="rhic200_cent0to10"){
		Rx=10.2; Ry=6.5; Rz=7.5;
	}
	else if(qualifier=="rhic200_cent10to15"){
		Rx=8.6; Ry=5.9; Rz=6.5;
	}
	else if(qualifier=="rhic200_cent10to20"){
		Rx=7.4; Ry=5.3; Rz=5.6;
	}
	else if(qualifier=="rhic200_cent15to20"){
		Rx=8.2; Ry=5.7; Rz=6.6;
	}
	else if(qualifier=="rhic200_cent20to30"){
		Rx=7.8; Ry=5.3; Rz=6.2;
	}
	else if(qualifier=="rhic200_cent30to40"){
		Rx=6.9; Ry=4.6; Rz=5.3;
	}
	else if(qualifier=="rhic200_cent40to50"){
		Rx=5.8; Ry=4.0; Rz=4.8;
	}
	else{
		printf("Don't recognize qualifier %s, can't initialize radii for hbt radius search\n",qualifier.c_str());
		exit(1);
	}
	
		//printf("qualifier=%s, R=(%g,%g,%g)\n",qualifier.c_str(),Rx,Ry,Rz);
	int *idlista,*idlistb,nida,nidb;
	if(HBTPAIRS=="PIPI"){
		idlista=new int[3];
		idlistb=new int[3];
		idlista[0]=idlistb[0]=211;
		idlista[1]=idlistb[1]=-211;
		idlista[2]=idlistb[2]=111;
		nida=nidb=3;
	}
	if(HBTPAIRS=="KK"){
		idlista=new int [4];
		idlistb=new int [4];
		idlista[0]=idlistb[0]=311;
		idlista[1]=idlistb[1]=-311;
		idlista[2]=idlistb[2]=321;
		idlista[3]=idlistb[3]=-321;
		nida=nidb=4;
	}
	int nmc,ikt,iphi,istarbin;
	char sdirname[120],cfdirname[120];
	double kt,spectra,wtot=0.0,w,m=139.57,sigma,ktbar;
	CMCPRList ***list=NULL;
		// Initialize Source Calc Object
	//#ifdef __B3D_USE_HDF5__
	//CSourceCalc_HDF5_MultiBin *scalc=new CSourceCalc_HDF5_MultiBin();
	//#else
	CSourceCalc_OSCAR_MultiBin *scalc=new CSourceCalc_OSCAR_MultiBin();
	scalc->NPHIBINS=6;
	scalc->NPTBINS=18;
	scalc->DELPT=25.0;
	scalc->PTMIN=150.0;
	scalc->OSCARfilename=oscarfilename;
	if(BINARY_RW)
		scalc->B3D_BINARY_FORMAT=true;
	//#endif

	scalc->PTMAX=scalc->PTMIN+scalc->NPTBINS*scalc->DELPT;
	scalc->DELPHI=90.0/double(scalc->NPHIBINS);
	parameter::PrintPars(scalc->spars);
	
	scalc->SetIDs(idlista,nida,idlistb,nidb);
	int nptbins=scalc->NPTBINS;
	int nphibins=scalc->NPHIBINS;
	double **Rout=new double *[nptbins];
	double **Rside=new double *[nptbins];
	double **Rlong=new double *[nptbins];
	double *Routbar=new double[nptbins];
	double *Rsidebar=new double[nptbins];
	double *Rlongbar=new double[nptbins];
	for(ikt=0;ikt<nptbins;ikt++){
		Rout[ikt]=new double[nphibins];
		Rside[ikt]=new double[nphibins];
		Rlong[ikt]=new double[nphibins];
	}
	for(istarbin=0;istarbin<4;istarbin++){
		Routbar[istarbin]=0.0;
		Rsidebar[istarbin]=0.0;
		Rlongbar[istarbin]=0.0;
	}
	double xoff,yoff,zoff,gamma;
		
	scalc->CalcS(list,list);
	
	istarbin=0;
		//for(ikt=0;ikt<nptbins;ikt++){
	for(ikt=0;ikt<nptbins;ikt++){
		kt=scalc->PTMIN+scalc->DELPT*(0.5+ikt);
		printf("XXXXXXXXXXXX kt=%g XXXXXXXXXXXXXXX\n",kt);
		gamma=sqrt(kt*kt+m*m)/m;
		for(iphi=0;iphi<scalc->NPHIBINS;iphi++){
			nmc=list[ikt][iphi]->GetNMC();
			spectra=double(nmc)/(kt*sqrt(kt*kt+m*m));
			w=double(nmc)*spectra;
			wtot+=w;
				//printf("---- kt=%g, ikt=%d, istarbin=%d, iphi=%d, nmc=%d -----\n",kt,ikt,istarbin,iphi,nmc);
			scalc->CalcEffGaussParsPureBose(list[ikt][iphi],lambda,Rx,Ry,Rz);
			Rout[ikt][iphi]=Rx;
			Rside[ikt][iphi]=Ry;
			Rlong[ikt][iphi]=Rz;
			Routbar[istarbin]+=Rout[ikt][iphi]*w/gamma;
			Rsidebar[istarbin]+=Rside[ikt][iphi]*w;
			Rlongbar[istarbin]+=Rlong[ikt][iphi]*w;
			sigma=0.5;
		}
		if(ikt==3 || ikt==7 || ikt==11 || ikt==17){
			Routbar[istarbin]=Routbar[istarbin]/wtot;
			Rsidebar[istarbin]=Rsidebar[istarbin]/wtot;
			Rlongbar[istarbin]=Rlongbar[istarbin]/wtot;
			wtot=0.0;
			ktbar=scalc->PTMIN+(0.5+istarbin)*100.0;
			printf("----------- istarbin=%d, ktbar=%g ---------------\n",istarbin,ktbar);
			
			printf("ROUT_PION_STAR = %g +/- %g\n",Routbar[istarbin],sigma);
			printf("RSIDE_PION_STAR = %g +/- %g\n",Rsidebar[istarbin],sigma);
			printf("RLONG_PION_STAR = %g +/- %g\n",Rlongbar[istarbin],sigma);
			if(WRITEDETAILS) fprintf(detailsfile,"%6.2f %6.4f %6.4f %6.4f\n",ktbar,Routbar[istarbin],Rsidebar[istarbin],Rlongbar[istarbin]);
			istarbin+=1;
		}
	}
	sigma=0.5;
	Rx=(Routbar[0]+Routbar[1]+Routbar[2]+Routbar[3])/4.0;
	Ry=(Rsidebar[0]+Rsidebar[1]+Rsidebar[2]+Rsidebar[3])/4.0;
	Rz=(Rlongbar[0]+Rlongbar[1]+Rlongbar[2]+Rlongbar[3])/4.0;
	if(HBTPAIRS=="PIPI"){
		fprintf(anal_output,"STAR_ROUT_PION %g %g\n",Rx,sigma);
		fprintf(anal_output,"STAR_RSIDE_PION %g %g\n",Ry,sigma);
		fprintf(anal_output,"STAR_RLONG_PION %g %g\n",Rz,sigma);
	}
	else if(HBTPAIRS=="KK"){
		fprintf(anal_output,"STAR_ROUT_KAON %g %g\n",Rx,sigma);
		fprintf(anal_output,"STAR_RSIDE_KAON %g %g\n",Ry,sigma);
		fprintf(anal_output,"STAR_RLONG_KAON %g %g\n",Rz,sigma);
	}

	for(ikt=0;ikt<nptbins;ikt++){
		for(iphi=0;iphi<scalc->NPHIBINS;iphi++){
			delete list[ikt][iphi];
		}
		delete [] list[ikt];
	}
	delete [] list;
		//Asourcebarbar->ScaleArray(1.0/wtottot);
	if(WRITEDETAILS){
		fflush(detailsfile);
		fclose(detailsfile);
	}
	for(ikt=0;ikt<nptbins;ikt++){
		delete [] Rout[ikt];
		delete [] Rside[ikt];
		delete [] Rlong[ikt];
	}
	delete [] Rout;
	delete [] Rside;
	delete [] Rlong;
	delete [] Routbar;
	delete [] Rsidebar;
	delete [] Rlongbar;
	delete [] idlista;
	delete [] idlistb;
	fclose(anal_output);
}

void CB3D::CalcHBT_ALICE(){
	string anal_output_filename;
	FILE *anal_output;
	double Rx=6.0,Ry=6.8,Rz=7.6,lambda=0.775;
	bool WRITEDETAILS=parameter::getB(parmap,"WRITEDETAILS",false);
	//string HBTPAIRS="KK";
	string HBTPAIRS="PIPI";
	string details_filename;
	FILE* detailsfile;
	
	string dirname,command;
	dirname="analysis/"+run_name+"/"+qualifier;
	command="mkdir -p "+dirname;
	system(command.c_str());
	anal_output_filename=dirname+"/results_alice_hbt.dat";
	anal_output=fopen(anal_output_filename.c_str(),"w");

	if(WRITEDETAILS){
		dirname="analysis/"+run_name+"/"+qualifier+"/details";
		command="mkdir -p "+dirname;
		system(command.c_str());
		details_filename=dirname+"/hbt_alice.dat";
		detailsfile=fopen(details_filename.c_str(),"w");
		fprintf(detailsfile,"# ALICE HBT RADII\n");
	}
	if(qualifier=="lhc2760_cent0to5"){
		Rx=10.0; Ry=6.8; Rz=7.6;
	}
	else if(qualifier=="lhc2760_cent5to10"){
		Rx=9.7; Ry=6.6; Rz=7.4;
	}
	else if(qualifier=="lhc2760_cent0to10"){
		Rx=10.2; Ry=6.5; Rz=7.5;
	}
	else if(qualifier=="lhc2760_cent10to15"){
		Rx=8.6; Ry=5.9; Rz=6.5;
	}
	else if(qualifier=="lhc2760_cent10to20"){
		Rx=7.4; Ry=5.3; Rz=5.6;
	}
	else if(qualifier=="lhc2760_cent15to20"){
		Rx=8.2; Ry=5.7; Rz=6.6;
	}
	else if(qualifier=="lhc2760_cent20to30"){
		Rx=7.8; Ry=5.3; Rz=6.2;
	}
	else if(qualifier=="lhc2760_cent30to40"){
		Rx=6.9; Ry=4.6; Rz=5.3;
	}
	else if(qualifier=="lhc2760_cent40to50"){
		Rx=5.8; Ry=4.0; Rz=4.8;
	}
	else{
		printf("Don't recognize qualifier, can't initialize radii for hbt radius search\n");
		exit(1);
	}
	Rx*=1.2; Ry*=1.2; Rz*=1.2;
	
		//printf("qualifier=%s, R=(%g,%g,%g)\n",qualifier.c_str(),Rx,Ry,Rz);
	int *idlista,*idlistb,nida,nidb;
	if(HBTPAIRS=="PIPI"){
		idlista=new int[3];
		idlistb=new int[3];
		idlista[0]=idlistb[0]=211;
		idlista[1]=idlistb[1]=-211;
		idlista[2]=idlistb[2]=111;
		nida=nidb=3;
	}
	if(HBTPAIRS=="KK"){
		idlista=new int [4];
		idlistb=new int [4];
		idlista[0]=idlistb[0]=311;
		idlista[1]=idlistb[1]=-311;
		idlista[2]=idlistb[2]=321;
		idlista[3]=idlistb[3]=-321;
		nida=nidb=4;
	}
	int nmc,ikt,iphi,ialicebin;
	char sdirname[120],cfdirname[120];
	double kt,spectra,wtot=0.0,w,m=139.57,sigma,ktbar;
	CMCPRList ***list=NULL;
		// Initialize Source Calc Object
	//#ifdef __B3D_USE_HDF5__
	//CSourceCalc_HDF5_MultiBin *scalc=new CSourceCalc_HDF5_MultiBin();
	//#else
	CSourceCalc_OSCAR_MultiBin *scalc=new CSourceCalc_OSCAR_MultiBin();
	scalc->NPHIBINS=6;
	scalc->NPTBINS=24;
	scalc->DELPT=25.0;
	scalc->PTMIN=200.0;
	scalc->OSCARfilename=oscarfilename;
	if(BINARY_RW)
		scalc->B3D_BINARY_FORMAT=true;
	//#endif

	scalc->PTMAX=scalc->PTMIN+scalc->NPTBINS*scalc->DELPT;
	scalc->DELPHI=90.0/double(scalc->NPHIBINS);
	parameter::PrintPars(scalc->spars);
	
	scalc->SetIDs(idlista,nida,idlistb,nidb);
	int nptbins=scalc->NPTBINS;
	int nphibins=scalc->NPHIBINS;
	double **Rout=new double *[nptbins];
	double **Rside=new double *[nptbins];
	double **Rlong=new double *[nptbins];
	double *Routbar=new double[nptbins];
	double *Rsidebar=new double[nptbins];
	double *Rlongbar=new double[nptbins];
	for(ikt=0;ikt<nptbins;ikt++){
		Rout[ikt]=new double[nphibins];
		Rside[ikt]=new double[nphibins];
		Rlong[ikt]=new double[nphibins];
	}
	for(ialicebin=0;ialicebin<6;ialicebin++){
		Routbar[ialicebin]=0.0;
		Rsidebar[ialicebin]=0.0;
		Rlongbar[ialicebin]=0.0;
	}
	double xoff,yoff,zoff,gamma;
		
	scalc->CalcS(list,list);
	
	ialicebin=0;
		//for(ikt=0;ikt<nptbins;ikt++){
	for(ikt=0;ikt<nptbins;ikt++){
		kt=scalc->PTMIN+scalc->DELPT*(0.5+ikt);
		printf("XXXXXXXXXXXX kt=%g XXXXXXXXXXXXXXX\n",kt);
		gamma=sqrt(kt*kt+m*m)/m;
		for(iphi=0;iphi<scalc->NPHIBINS;iphi++){
			nmc=list[ikt][iphi]->GetNMC();
			spectra=double(nmc)/(kt*sqrt(kt*kt+m*m));
			w=double(nmc)*spectra;
			wtot+=w;
				//printf("---- kt=%g, ikt=%d, ialicebin=%d, iphi=%d, nmc=%d -----\n",kt,ikt,ialicebin,iphi,nmc);
			scalc->CalcEffGaussParsPureBose(list[ikt][iphi],lambda,Rx,Ry,Rz);
			Rout[ikt][iphi]=Rx;
			Rside[ikt][iphi]=Ry;
			Rlong[ikt][iphi]=Rz;
			Routbar[ialicebin]+=Rout[ikt][iphi]*w/gamma;
			Rsidebar[ialicebin]+=Rside[ikt][iphi]*w;
			Rlongbar[ialicebin]+=Rlong[ikt][iphi]*w;
			sigma=0.5;
		}
		if(ikt==3 || ikt==7 || ikt==11 || ikt==15 || ikt==19 || ikt==23){
			Routbar[ialicebin]=Routbar[ialicebin]/wtot;
			Rsidebar[ialicebin]=Rsidebar[ialicebin]/wtot;
			Rlongbar[ialicebin]=Rlongbar[ialicebin]/wtot;
			wtot=0.0;
			ktbar=scalc->PTMIN+(0.5+ialicebin)*100.0;
			printf("----------- ialicebin=%d, ktbar=%g ---------------\n",ialicebin,ktbar);
			
			printf("ROUT_PION_ALICE = %g +/- %g\n",Routbar[ialicebin],sigma);
			printf("RSIDE_PION_ALICE = %g +/- %g\n",Rsidebar[ialicebin],sigma);
			printf("RLONG_PION_ALICE = %g +/- %g\n",Rlongbar[ialicebin],sigma);
			if(WRITEDETAILS) fprintf(detailsfile,"%6.2f %6.4f %6.4f %6.4f\n",ktbar,Routbar[ialicebin],Rsidebar[ialicebin],Rlongbar[ialicebin]);
			ialicebin+=1;
		}
	}
	sigma=0.3;
	Rx=(Routbar[0]+Routbar[1]+Routbar[2]+Routbar[3]+Routbar[4]+Routbar[5])/6.0;
	Ry=(Rsidebar[0]+Rsidebar[1]+Rsidebar[2]+Rsidebar[3]+Rsidebar[4]+Rsidebar[5])/6.0;
	Rz=(Rlongbar[0]+Rlongbar[1]+Rlongbar[2]+Rlongbar[3]+Rlongbar[4]+Rlongbar[5])/6.0;
	if(HBTPAIRS=="PIPI"){
		fprintf(anal_output,"ALICE_ROUT_PION %g %g\n",Rx,sigma);
		fprintf(anal_output,"ALICE_RSIDE_PION %g %g\n",Ry,sigma);
		fprintf(anal_output,"ALICE_RLONG_PION %g %g\n",Rz,sigma);
	}
	else if(HBTPAIRS=="KK"){
		fprintf(anal_output,"ALICE_ROUT_KAON %g %g\n",Rx,sigma);
		fprintf(anal_output,"ALICE_RSIDE_KAON %g %g\n",Ry,sigma);
		fprintf(anal_output,"ALICE_RLONG_KAON %g %g\n",Rz,sigma);
	}

	for(ikt=0;ikt<nptbins;ikt++){
		for(iphi=0;iphi<scalc->NPHIBINS;iphi++){
			delete list[ikt][iphi];
		}
		delete [] list[ikt];
	}
	delete [] list;
		//Asourcebarbar->ScaleArray(1.0/wtottot);
	if(WRITEDETAILS){
		fflush(detailsfile);
		fclose(detailsfile);
	}
	for(ikt=0;ikt<nptbins;ikt++){
		delete [] Rout[ikt];
		delete [] Rside[ikt];
		delete [] Rlong[ikt];
	}
	delete [] Rout;
	delete [] Rside;
	delete [] Rlong;
	delete [] Routbar;
	delete [] Rsidebar;
	delete [] Rlongbar;
	delete [] idlista;
	delete [] idlistb;
	
	fclose(anal_output);
}

#endif
