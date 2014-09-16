#ifndef __INCLUDE_SPECIES_CC__
#define __INCLUDE_SPECIES_CC__
#include "species.h"
// CLASS CSpecies *********************************************

using namespace std;

void CSpecies::Print(){
	int iQ,i;
	printf("Nspecies=%d, Ncharges=%d\n",Nspecies,Ncharges);
	for(i=0;i<Nspecies;i++){
		printf("i=%2d: m=%g, degen=%g",i,m[i],degen[i]);
		for(iQ=0;iQ<Ncharges;iQ++) printf(", Q_%d=%4.1f",iQ,Q[iQ][i]);
		printf("\n");
	}
}

void CSpecies::InitStandardHadrons(){

	NID=new int[Nspecies];
	ID=new int*[Nspecies];
	m[0]=139.57; degen[0]=3; //pions
	m[1]=497.0; degen[1]=4;  //Kaons
	m[2]=549.0; degen[2]=1; // eta
	m[3]=770.0; degen[3]=9; // rho
	m[4]=783.0; degen[4]=3; // omega
	m[5]=892.0; degen[5]=12; //K*
	m[6]=958.0; degen[6]=3; //eta'
	m[7]=1020.0; degen[7]=3; //phi
	m[8]=939.0; degen[8]=8; // nucleons
	m[9]=1232.0; degen[9]=32; // Delta
	m[10]=1116.0; degen[10]=4; // Lambda
	m[11]=1193.0; degen[11]=12;  // Sigma
	m[12]=1383.0; degen[12]=24;  //Sigma *
	m[13]=1315.0; degen[13]=8; //Xi
	m[14]=1530.0; degen[14]=16;  //Xi*
	m[15]=1672.0; degen[15]=8;  //Omega

	NID[0]=3; ID[0]=new int[NID[0]];
	ID[0][0]=211; ID[0][1]=111; ID[0][2]=-211;

	NID[1]=4; ID[1]=new int[NID[1]];
	ID[1][0]=311; ID[1][1]=-311; ID[1][2]=321; ID[1][3]=-321;

	NID[2]=1; ID[2]=new int[NID[2]];
	ID[2][0]=221;

	NID[3]=3; ID[3]=new int[NID[3]];
	ID[3][0]=213; ID[3][1]=113; ID[3][2]=-213;

	NID[4]=1; ID[4]=new int[NID[4]];
	ID[4][0]=223;

	NID[5]=4; ID[5]=new int[NID[5]];
	ID[5][0]=313; ID[5][1]=-313; ID[5][2]=323; ID[5][3]=-323;

	NID[6]=1; ID[6]=new int[NID[6]];
	ID[6][0]=331;

	NID[7]=1; ID[7]=new int[NID[7]];
	ID[7][0]=333;

	NID[8]=4; ID[8]=new int[NID[8]];
	ID[8][0]=2212; ID[8][1]=-2212; ID[8][2]=2112; ID[8][3]=-2112;

	NID[9]=8; ID[9]=new int[NID[9]];
	ID[9][0]=2224; ID[9][1]=-2224; ID[9][2]=2214; ID[9][3]=-2214;
	ID[9][4]=2114; ID[9][5]=-2114; ID[9][6]=1114; ID[9][7]=-1114;

	NID[10]=2; ID[10]=new int[NID[10]];
	ID[10][0]=3122; ID[10][1]=-3122;

	NID[11]=6; ID[11]=new int[NID[11]];
	ID[11][0]=3222; ID[11][1]=-3222; ID[11][2]=3212; ID[11][3]=-3212;
	ID[11][4]=3112; ID[11][5]=-3112;

	NID[12]=6; ID[12]=new int[NID[12]];
	ID[12][0]=3224; ID[12][1]=-3224; ID[12][2]=3214; ID[12][3]=-3214;
	ID[12][4]=3114; ID[12][5]=-3114;

	NID[13]=4; ID[13]=new int[NID[13]];
	ID[13][0]=3322; ID[13][1]=-3322; ID[13][2]=3312; ID[13][3]=-3312;

	NID[14]=4; ID[14]=new int[NID[14]];
	ID[14][0]=3324; ID[14][1]=-3324; ID[14][2]=3314; ID[14][3]=-3314;

	NID[15]=2; ID[15]=new int[NID[15]];
	ID[15][0]=3334; ID[15][1]=-3334;
}

CSpecies_QGP::CSpecies_QGP(){
	Nspecies=Ncharges=0;
}

CSpecies_RelGas::CSpecies_RelGas(){
	int iQ,i;
	Nspecies=1;
	Ncharges=1;
	m=new double[Nspecies];
	degen=new double[Nspecies];
	NID=new int[Nspecies];
	Q=new double *[Ncharges];
	for(iQ=0;iQ<Ncharges;iQ++){
		Q[iQ]=new double[Nspecies];
		for(i=0;i<Nspecies;i++) Q[iQ][i]=0.0;
	}
	m[0]=1.0E-10; degen[0]=1;
	for(iQ=0;iQ<Ncharges;iQ++) Q[iQ][iQ]=1.0;
}

CSpecies_StandardHadrons_Equil::CSpecies_StandardHadrons_Equil(){
	int iQ,i;
	Nspecies=16;
	Ncharges=0;
	m=new double[Nspecies];
	degen=new double[Nspecies];
	InitStandardHadrons();

}

CSpecies_StandardHadrons::CSpecies_StandardHadrons(){
	int iQ,i;
	Nspecies=16;
	Ncharges=16;
	m=new double[Nspecies];
	degen=new double[Nspecies];
	InitStandardHadrons();
	for(iQ=0;iQ<Ncharges;iQ++) Q[iQ][iQ]=1.0;
}

CSpecies_StandardHadrons5Q::CSpecies_StandardHadrons5Q(){
	using namespace std;
	int iQ,i;
	Nspecies=16;
	Ncharges=5;
	m=new double[Nspecies];
	degen=new double[Nspecies];
	InitStandardHadrons();
	Q=new double *[Ncharges];
	for(iQ=0;iQ<Ncharges;iQ++){
		Q[iQ]=new double[Nspecies];
		for(i=0;i<Nspecies;i++) Q[iQ][i]=0.0;
	}
	// Baryon+Antibaryon # Conserved
	Q[0][8]=Q[0][9]=Q[0][10]=Q[0][11]=Q[0][12]=Q[0][13]=Q[0][14]=1.0;
	// Strange+Antistrange # Conserved
	Q[1][1]=Q[1][5]=Q[1][10]=Q[1][11]=Q[1][12]=1.0;
	Q[1][7]=Q[1][13]=Q[1][14]=2.0;
	Q[1][15]=3.0;
	// Pion # Conserved
	Q[2][0]=1.0;
	Q[2][3]=2.0;
	Q[2][5]=1.0;
	Q[2][9]=1.0;
	Q[2][12]=1.0;
	Q[2][14]=1.0;
	// Eta # Conserved
	Q[3][2]=Q[3][6]=1.0;
	// Omega # Conserved
	Q[4][4]=1.0;
}

CSpecies_PionsOnly::CSpecies_PionsOnly(){
	int iQ,i;
	Nspecies=1;
	Ncharges=0;
	m=new double[Nspecies];
	degen=new double[Nspecies];
	NID=new int[Nspecies];
	ID=new int*[Nspecies];
	Q=new double *[Ncharges];
	for(iQ=0;iQ<Ncharges;iQ++){
		Q[iQ]=new double[Nspecies];
		for(i=0;i<Nspecies;i++) Q[iQ][i]=0.0;
	}
	m[0]=139.57; degen[0]=3;
	// Conserve Pion Number
	
	NID[0]=3; ID[0]=new int[NID[0]];
	ID[0][0]=211; ID[0][1]=111; ID[0][2]=-211;

}

#endif
