#ifndef __COLLIDE_CC__
#define __COLLIDE_CC__
#include "b3d.h"
using namespace std;
//#define __SCATTERING_ON__

// If two particles pass one another, Collide will determine whether to scatter and how

int CB3D::Collide(CPart *part1,CPart *part2){
	CPartMap::iterator ppos;
	const double g[4]={1,-1,-1,-1};
	double sigma=0.0,sigma_annihilation,sigma_inel,Gamma,G,G2,MR,M,m1,m2,b,q2,q3,q4,qR2,tan2delta,scompare;
	double mt,P[4],q[4],r[4],qdotr,P2,Pdotq,Pdotr,rsquared;
	const int NWMAX=5000;
	double weight[NWMAX]={0.0};
	double inel_d=0.0,q_prime,temp;
	int ir1,ir2,irflip,alpha,G_Value,L_merge, pmq=0, pmb=0, pms=0;
	CMerge *merge;
	list<CInelasticInfo>::iterator inel;
	list<CInelasticInfo> inel_list;
	bool G_Parity = false;
	int netb=0,netq=0,nets=0,ndaughters;
	int qpions,iK,ipair,npaircheck;
	double Plab,p1dotp2;
	int itau;
	double jR,j1,j2,j1_i,j2_i,rstrange;

	ir1=part1->resinfo->ires; ir2=part2->resinfo->ires;
	if(ir1>ir2){
		irflip=ir1; ir1=ir2; ir2=irflip;
	}
	bool bjtranslate=false;
	if(BJORKEN && ((part1->cell->ieta==0 && part2->cell->ieta==2*NETA-1) || (part1->cell->ieta==2*NETA-1 && part2->cell->ieta==0))){
		bjtranslate=true;
		part1->BjorkenTranslate();
	}

	//Calc PI*b^2
	P2=Pdotq=Pdotr=0.0;
	for(alpha=0;alpha<4;alpha++){
		P[alpha]=part1->p[alpha]+part2->p[alpha];
		q[alpha]=part1->p[alpha]-part2->p[alpha];
		r[alpha]=part1->r[alpha]-part2->r[alpha];
		P2+=g[alpha]*P[alpha]*P[alpha];
		Pdotq+=g[alpha]*P[alpha]*q[alpha];
		Pdotr+=g[alpha]*P[alpha]*r[alpha];
	}
	q2=qdotr=rsquared=0.0;
	for(alpha=0;alpha<4;alpha++){
		q[alpha]-=Pdotq*P[alpha]/P2;
		r[alpha]-=Pdotr*P[alpha]/P2;
		q2+=g[alpha]*q[alpha]*q[alpha];
		qdotr+=g[alpha]*q[alpha]*r[alpha];
	}
	for(alpha=0;alpha<4;alpha++){
		r[alpha]-=qdotr*q[alpha]/q2;
		rsquared-=g[alpha]*r[alpha]*r[alpha];
	}
	scompare=PI*rsquared;

	// Use Fixed Xsection for s-wave scattering
	sigma+=SIGMADEFAULT/double(NSAMPLE);
	if(scompare<sigma){
		Scatter(part1,part2);
		if(bjtranslate)
			part1->BjorkenUnTranslate();
		part1->tau_lastint=part1->tau0;
		part1->actionmother=nactions;
		part1->KillActions();
		part1->FindActions();
		part2->tau_lastint=part2->tau0;
		part2->actionmother=nactions;
		part2->KillActions();
		part2->FindActions();
		return 2;
	}

	m2=part2->GetMass();
	m1=part1->GetMass();
	if(ANNIHILATION_CHECK && (part1->resinfo->baryon*part2->resinfo->baryon)<0){
		p1dotp2=0.25*(P2-q2);
		Plab=sqrt((p1dotp2*p1dotp2/(m2*m2))-m1*m1);
		//printf("Plab=%g, sqrt(p1dotp2)=%g, m1=%g, m2=%g\n",Plab,sqrt(p1dotp2),m1,m2);
		sigma_annihilation=6.7*pow(Plab/1000.0,-0.7)/double(NSAMPLE);
		rstrange=0.5*sqrt(sigma_annihilation);
		rstrange*=pow(ANNIHILATION_SREDUCTION,abs(part1->resinfo->strange))+pow(ANNIHILATION_SREDUCTION,abs(part2->resinfo->strange));
		sigma_annihilation=rstrange*rstrange;
		if(scompare<sigma && tau<TAUCOLLMAX){
			itau=lrint(floor(tau));
			annihilation_array[itau]+=1.0;
			ndaughters=Annihilate(part1,part2);
			return 4;
		}
	}

	//Calculate quantities used for both inel scattering and merging
	merge=reslist->MergeArray[ir1][ir2];
	if(merge!=NULL||inel!=inel_list.end()){
		j1=part1->resinfo->spin;
		j2=part2->resinfo->spin;
		M=0.0;
		for(alpha=0;alpha<4;alpha++) M+=g[alpha]*pow(part1->p[alpha]+part2->p[alpha],2);
			M=sqrt(M);
		if(M<m1+m2){
			printf("SCREWY MASSES: m1=%g, m2=%g, M=%g\n",m1,m2,M);
			//part1->Print();
			//part2->Print();
			q2=0.0;
		}
		else q2=Misc::triangle(M,m1,m2);
	}

	//Check for merging
	while(merge!=NULL){
		Gamma=merge->resinfo->width;
		b=merge->branching;
		jR=merge->resinfo->spin;
		MR=merge->resinfo->mass;
		L_merge = merge->L;


		if(m1+m2<MR){
			qR2=Misc::triangle(MR,m1,m2);
			q3=pow(q2/qR2,(2*L_merge + 1)/2);
			q4=pow(q2/qR2,(2*L_merge)/2);
			G=Gamma*(MR/M)*q3*1.2/(1.0+0.2*q4);
			tan2delta=0.25*G*G/((M-MR)*(M-MR));

			sigma+=b*((4.0*PI*HBARC*HBARC/q2)*(tan2delta/(1.0+tan2delta))
				*((2.0*jR+1.0)/((2.0*j1+1.0)*(2.0*j2+1.0))))/double(NSAMPLE);
			if(sigma>scompare){
				if (Merge(part1,part2,merge->resinfo)){
					if(bjtranslate)
						part1->BjorkenUnTranslate();
					part1->tau_lastint=part1->tau0;
					part1->actionmother=nactions;
					part1->KillActions();
					part1->FindActions();
					return 1;
				}
			}
		}
		merge=merge->next;
	}

	//Check for Inelastic Scatering
	//inel_d = (2.0*j1+1.0)*(2.0*j2+1.0)*q2;
	if(INELASTIC){
		if(sigma+(SIGMAINELASTIC/double(NSAMPLE))>scompare){
			if(part1->resinfo->G_Parity && part2->resinfo->G_Parity){
				G_Parity = true;
				G_Value = part1->resinfo->G_Parity * part2->resinfo->G_Parity;
			}
			// First calculate denominator
			int iw;

			if(inelasticlist->UseInelasticArray){
				inel_list = inelasticlist->InelasticArray[ir1][ir2];
			}else{
				netq = part1->resinfo->charge+part2->resinfo->charge;
				netb = part1->resinfo->baryon+part2->resinfo->baryon;
				nets = part1->resinfo->strange+part2->resinfo->strange;
				inel_list = inelasticlist->ThermalArray[abs(netb)][abs(netq)][abs(nets)][Misc::Sign(netb)][Misc::Sign(netq)][Misc::Sign(nets)];
			}
			inel_d = Q0*Q0;
			inel = inel_list.begin();
			iw=0;
			while(inel!=inel_list.end() && (M>inel->min_mass)){
				if((G_Parity && (inel->resinfo_1->G_Parity * inel->resinfo_2->G_Parity == G_Value)) || (!G_Parity)){
					j1_i=inel->resinfo_1->spin;
					j2_i=inel->resinfo_2->spin;
					if(inel->resinfo_1->mass+inel->resinfo_2->mass<M){
						//if(inel->resinfo_1->baryon==0 || inel->resinfo_2->baryon==0){
						q_prime = Misc::triangle(M,inel->resinfo_1->mass, inel->resinfo_2->mass);
						weight[iw] = (2.0*j1_i+1.0)*(2.0*j2_i+1.0)*q_prime;
						//}
					}
					inel_d += weight[iw];
				}
				iw+=1;
				if(iw==NWMAX){
					printf("MUST INCREASE NWMAX in int CB3D::Collide\n");
					exit(1);
				}
				inel++;
			}

			// now thumb through
			inel = inel_list.begin();
			iw=0;
			while(inel!=inel_list.end() && (M>inel->min_mass)){
				if((G_Parity && (inel->resinfo_1->G_Parity * inel->resinfo_2->G_Parity == G_Value)) || !G_Parity){
					sigma+=SIGMAINELASTIC*weight[iw]/(inel_d*double(NSAMPLE));
					if(sigma > scompare){
						netb=part1->resinfo->baryon+part2->resinfo->baryon;
						netq=part1->resinfo->charge+part2->resinfo->charge;
						nets=part1->resinfo->strange+part2->resinfo->strange;
						InelasticScatter(part1,part2, *inel);
						if(bjtranslate)
							part1->BjorkenUnTranslate();
						part1->tau_lastint=part1->tau0;
						part1->actionmother=nactions;
						part1->FindActions();
						part2->tau_lastint=part2->tau0;
						part2->actionmother=nactions;
						part2->FindActions();
						netb-=part1->resinfo->baryon+part2->resinfo->baryon;
						netq-=part1->resinfo->charge+part2->resinfo->charge;
						nets-=part1->resinfo->strange+part2->resinfo->strange;
						if(netb!=0 || netq!=0 || nets!=0){
							printf("charge not conserved in inel collision, netb=%d, netq=%d, nets=%d\n",netb,netq,nets);
						}
						return 3;
					}
				}
				iw+=1;
				inel++;
			}
		}
	}

	// Nothing happens
	part1->tau_lastint=part1->tau0;
 	part1->actionmother=nactions;
 	part2->tau_lastint=part2->tau0;
 	part2->actionmother=nactions;
 	if(bjtranslate)
 		part1->BjorkenUnTranslate();
 	
	return 0;
}

#endif