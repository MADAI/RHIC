#ifndef __COLLIDE_CC__
#define __COLLIDE_CC__
#include "b3d.h"
//#define __SCATTERING_ON__

// If two particles pass one another, Collide will determine whether to scatter and how

int CB3D::Collide(CPart *part1,CPart *part2,int &nproducts,array<CPart*,5> &product,double pibsquared){
	CPartMap::iterator ppos;
	CPart *part3,*part4;
	const double g[4]={1,-1,-1,-1};
	double sigma=0.0,sigma_annihilation,sigma_inel,Gamma,G,G2,MR,M,b,q2,q3,q4,qR2;
	double tan2delta;
	double mt,P2,Pdotq,vrel;
	const int NWMAX=5000;
	double inel_weight[NWMAX]={0.0};
	double inel_d=0.0,q_prime;
	int ir1,ir2,irflip,alpha,G_Value,L_merge, pmq=0, pmb=0, pms=0;
	CMerge *merge;
	list<CInelasticInfo>::iterator inel;
	list<CInelasticInfo> inel_list;
	bool G_Parity = false;
	int netb=0,netq=0,nets=0,ndaughters;
	int qpions,iK,ipair,npaircheck;
	double Plab,p1dotp2;
	int itau;
	double jR,j1,j2,j1_i,j2_i,rstrange,Pcheck;

	ir1=part1->resinfo->ires; ir2=part2->resinfo->ires;
	if(ir1>ir2){
		irflip=ir1; ir1=ir2; ir2=irflip;
	}
	bool bjtranslate=false;
	if(BJORKEN && ((part1->cell->ieta==0 && part2->cell->ieta==2*NETA-1) || (part1->cell->ieta==2*NETA-1 && part2->cell->ieta==0))){
		bjtranslate=true;
		part1->BjorkenTranslate();
	}
	merge=reslist->MergeArray[ir1][ir2];
	if(merge!=NULL || (BARYON_ANNIHILATION && (part1->resinfo->baryon*part2->resinfo->baryon)<0) || INELASTIC){
		p1dotp2=0.0;
		for(alpha=0;alpha<4;alpha++){
			p1dotp2+=part1->p[alpha]*part2->p[alpha]*g[alpha];
		}
		P2=part1->msquared+part2->msquared+2.0*p1dotp2;
		q2=Misc::triangle2(P2,part1->msquared,part2->msquared);
		M=sqrt(P2);
	}
	
	// Use Fixed Xsection for s-wave scattering
	sigma+=SIGMADEFAULT/double(NSAMPLE);
	if(pibsquared<sigma){
		part3=product[0];
		part4=product[1];
		part3->CopyPositionInfo(part1);
		part4->CopyPositionInfo(part2);
		nproducts=2;

		Scatter(part1,part2,part3,part4);
		
		if(bjtranslate){
			part1->BjorkenUnTranslate();
			part3->BjorkenUnTranslate();
		}
		return 2;
	}
	//
	// Annihilation
	if(BARYON_ANNIHILATION && (part1->resinfo->baryon*part2->resinfo->baryon)<0){
		sigma_annihilation=GetAnnihilationSigma(part1,part2,vrel);
		sigma+=sigma_annihilation;
		if(pibsquared<sigma && tau<TAUCOLLMAX){
			itau=lrint(floor(tau));
			//annihilation_array[itau]+=1.0;
			Annihilate(part1,part2,nproducts,product);
			return 4;
		}
	}
	
	//Calculate quantities used for both inel scattering and merging
	if(merge!=NULL || inel!=inel_list.end()){
		j1=part1->resinfo->spin;
		j2=part2->resinfo->spin;
	}

	//Check for merging
	
	while(merge!=NULL){
		Gamma=merge->resinfo->width;
		b=merge->branching;
		jR=merge->resinfo->spin;
		MR=merge->resinfo->mass;
		L_merge = merge->L;
		qR2=Misc::triangle2(MR*MR,part1->msquared,part2->msquared);
		q3=pow(q2/qR2,(2*L_merge + 1)/2);
		q4=pow(q2/qR2,(2*L_merge)/2);
		//q3=q4=1.0;
		G=Gamma*(MR/M)*q3*1.2/(1.0+0.2*q4);
		tan2delta=pow(0.5*G/(M-MR),2);

		sigma+=b*((4.0*PI*HBARC*HBARC/q2)*(tan2delta/(1.0+tan2delta))
			*((2.0*jR+1.0)/((2.0*j1+1.0)*(2.0*j2+1.0))))/double(NSAMPLE);
		if(sigma>pibsquared){
			part3=product[0];
			part3->CopyPositionInfo(part1);
			if(Merge(part1,part2,part3,merge->resinfo)){
				if(bjtranslate){
					part1->BjorkenUnTranslate();
					part3->BjorkenUnTranslate();
				}
				nproducts=1;
				return 1;
			}
			else{
				if(bjtranslate)
					part1->BjorkenUnTranslate();
				nproducts=0;
				return 0;
			}
		}
		merge=merge->next;
	}
	
	//Check for Inelastic Scatering
	//inel_d = (2.0*j1+1.0)*(2.0*j2+1.0)*q2;
	if(INELASTIC){
		if(sigma+(SIGMAINELASTIC/double(NSAMPLE))>pibsquared){
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
						inel_weight[iw] = (2.0*j1_i+1.0)*(2.0*j2_i+1.0)*q_prime;
						//}
					}
					inel_d += inel_weight[iw];
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
					sigma+=SIGMAINELASTIC*inel_weight[iw]/(inel_d*double(NSAMPLE));
					if(sigma>pibsquared){
						netb=part1->resinfo->baryon+part2->resinfo->baryon;
						netq=part1->resinfo->charge+part2->resinfo->charge;
						nets=part1->resinfo->strange+part2->resinfo->strange;
						InelasticScatter(part1,part2,part3,part4,*inel);
						if(bjtranslate){
							part1->BjorkenUnTranslate();
							part3->BjorkenUnTranslate();
						}
						if(netb!=0 || netq!=0 || nets!=0){
							printf("charge not conserved in inel collision, netb=%d, netq=%d, nets=%d\n",netb,netq,nets);
						}
						nproducts=2;
						return 3;
					}
				}
				iw+=1;
				inel++;
			}
		}
	}
	// Nothing happens
	if(bjtranslate)
		part1->BjorkenUnTranslate();
	nproducts=0;
	return 0;
}

#endif