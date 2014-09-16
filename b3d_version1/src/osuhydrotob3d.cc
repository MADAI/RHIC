#ifndef __hydrotob3d_cc__
#define __hydrotob3d_cc__
#include "b3d.h"
using namespace std;

CB3D *CFreezeoutCell::b3d=NULL;
COSUHydrotoB3D *CFreezeoutCell::osuhydrotob3d=NULL;

void COSUHydrotoB3D::Init(){
	if(initialization==true){
		printf("COSUHydrotoB3D already initialized\n");
		exit(1);
	}
	int ires,iring;
	string inputfilename;
	double pi,ei,sigma2i,dedti,degen;
	T=1000.0*parameter::getD(b3d->parmap,"HYDRO_FOTEMP",160.0);
	netdens=0.0;
	ETAMAX=b3d->ETAMAX;
	randy=b3d->randy;
	MC_Ntarget=0.0;
	MC_Nbar=0.0;
	MC_NWrite=0;
	reslist=b3d->reslist;
	NRES=reslist->NResonances;
	//printf("epsilon_H=%g\n",intrinsic->epsilon);
	density=new double[NRES];
	epsilon=P=0.0;
	CResInfo *resinfo=reslist->GfirstResInfoptr;
	for(ires=0;ires<NRES;ires++){
		if(resinfo->code!=22){
			degen=2.0*resinfo->spin+1.0;
			//printf("ires=%d, ID=%d, degen=%g, mass=%g, T=%g\n",ires,resinfo->code,degen,resinfo->mass,T);
			b3d->freegascalc_onespecies(resinfo->mass,T,pi,ei,density[ires],sigma2i,dedti);
			//printf("mass=%g, T=%g, pi=%g, ei=%g, di=%g\n",resinfo->mass,T,pi,ei,density[ires]);
			density[ires]*=degen;
			netdens+=density[ires];
			P+=degen*pi;
			epsilon+=degen*ei;
		}
		resinfo=resinfo->nextResInfoptr;
	}
	GetLambdaFact();
	//TestLambdaFact();
	initialization=true;
	NCELLS=0;
	cell.resize(10000);
	cell[0].b3d=b3d;
	cell[0].osuhydrotob3d=this;
}

void COSUHydrotoB3D::ReadInput(){
	if(!initialization) Init();
	string inputfilename="output/"+b3d->run_name+"/"+b3d->qualifier+"/osuFOS.dat";
	FILE *input=fopen(inputfilename.c_str(),"r");
	char dummy[200];
	int icell=0;
	while(!feof(input)){
		if(icell==cell.size()) cell.resize(cell.size()+10000);
		if(b3d->HYDRO_PURE_BJORKEN){
			cell[icell].Read(input);	
		}
		else{
			cell[icell].Read3D(input);	
		}
		icell+=1;
	}
	NCELLS=icell;
	//printf("NCELLS=%lld\n",NCELLS);
	fclose(input);
}

bool CFreezeoutCell::Read(FILE *fptr){
	double pi33,pi00,pi01,pi02,pi11,pi12,pi22,**pi;
	double V,V2,Vmag,gamma,udotV;
	int alpha,beta;
	double vx,vy,eps_dec,T_dec,P_dec,muB_dec,muS_dec,bdens_dec,b;
	bool success=false;
	pi=new double *[4];
	char dummy[200];
	for(alpha=0;alpha<4;alpha++){
		pi[alpha]=new double[4];
	}
	for(alpha=0;alpha<4;alpha++){
		u[alpha]=0.0;
		r[alpha]=0.0;
		for(beta=0;beta<4;beta++){
			pi[alpha][beta]=0.0;
		}
	}
	fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&r[0],&r[1],&r[2],&Omega[0],&Omega[1],&Omega[2],&vx,&vy,&eps_dec,&bdens_dec,&T_dec,&muB_dec,&muS_dec,&P_dec,&pi[3][3],&pi[0][0],&pi[0][1],&pi[0][2],&pi[1][1],&pi[1][2],&pi[2][2]);
	if(!feof(fptr)){
		success=true;
		fscanf(fptr,"%lf",&b);
	//fgets(dummy,150,fptr);
		T=1000.0*T_dec;
		u[0]=1.0/sqrt(1.0-vx*vx-vy*vy);
		u[1]=u[0]*vx;
		u[2]=u[0]*vy;
		u[3]=0.0;
	//printf("pixx=%g, piyy=%g, pixy=%g, pizz=%g\n",pi[1][1],pi[2][2],pi[1][2],pi[3][3]);
	//printf("u=(%g,%g,%g,%g), tau^2=%g\n",u[0],u[1],u[2],u[3],r[0]*r[0]-r[3]*r[3]);
	//Omega[1]=-Omega[1]*2.0*r[0]*b3d->ETAMAX;
	//Omega[2]=-Omega[2]*2.0*r[0]*b3d->ETAMAX;
	//Omega[0]*=2.0*r[0]*b3d->ETAMAX;
		Omega[1]=-Omega[1]*2.0*b3d->ETAMAX*r[0];
		Omega[2]=-Omega[2]*2.0*b3d->ETAMAX*r[0];
		Omega[0]*=2.0*b3d->ETAMAX*r[0];
		Omega[3]=0.0;
		pi[2][1]=pi[1][2]; pi[1][0]=pi[0][1]; pi[2][0]=pi[0][2];
		Misc::BoostToCM(u,pi,dToverH);

		for(alpha=0;alpha<4;alpha++){
			for(beta=0;beta<4;beta++){
				dToverH[alpha][beta]=dToverH[alpha][beta]/(eps_dec+P_dec);
			//printf("%g ",dToverH[alpha][beta]);
			}
		//printf("\n");
		}

		for(alpha=0;alpha<4;alpha++) delete [] pi[alpha];
			delete [] pi;

		for(alpha=0;alpha<4;alpha++){
			for(beta=0;beta<4;beta++){
				lambda[alpha][beta]=0.0;
			}
		}
		for(alpha=1;alpha<4;alpha++){
			for(beta=1;beta<4;beta++){
				lambda[alpha][beta]=dToverH[alpha][beta]/osuhydrotob3d->lambdafact;
			}
			lambda[alpha][alpha]+=1.0;
		}

		if(b3d->HYDRO_OCTANT_SYMMETRY==2) for(alpha=0;alpha<4;alpha++) Omega[alpha]*=4.0;

		V2=Omega[0]*Omega[0]-Omega[1]*Omega[1]-Omega[2]*Omega[2]-Omega[3]*Omega[3];
		udotV=u[0]*Omega[0]-u[1]*Omega[1]-u[2]*Omega[2];
		V=sqrt(fabs(V2));
		if(V2>0){
			gamma=fabs(udotV/V);
		//printf("timelike Omega=(%g,%g,%g,%g), tau=%g\n",Omega[0],Omega[1],Omega[2],Omega[3],r[0]);
		}
		else{
			Vmag=sqrt(Omega[1]*Omega[1]+Omega[2]*Omega[2]+Omega[3]*Omega[3]);
			gamma=(u[0]*Vmag-(Omega[0]/Vmag)*(Omega[1]*u[1]+Omega[2]*u[2]+Omega[3]*u[3]))/V;
		//printf("spacelike Omega=(%g,%g,%g,%g), tau=%g, rvec.Omegavec=%g\n",Omega[0],Omega[1],Omega[2],Omega[3],r[0],Omega[1]*r[1]+Omega[2]*r[2]+Omega[3]*r[3]);
		}
		if(gamma<1.0){
			printf("Disaster: gamma =%g, but should be >1\n",gamma);
			printf("V2=%g, V=%g, gamma=%g, udotV=%g wmax=%g\n",V2,V,gamma,udotV,wmax);
			printf("fabs(MC_Nbar)=%g\n",fabs(osuhydrotob3d->MC_Nbar));
			exit(1);
		}

		wmax=V*(gamma+sqrt(gamma*gamma-1.0));
		//printf("V=%g, wmax=%g\n",V,wmax);
		return success;
	}
	else{
		return false;
	}
}

bool CFreezeoutCell::Read3D(FILE *fptr){
	double pi33,pi00,pi01,pi02,pi11,pi12,pi22,**pi;
	double V,V2,Vmag,gamma,udotV,b;
	bool success=false;
	int alpha,beta;
	double eps_dec,T_dec,P_dec,muB_dec,muS_dec,bdens_dec;
	double tau,eta,y,u_eta;
	pi=new double *[4];
	char dummy[200];
	for(alpha=0;alpha<4;alpha++){
		pi[alpha]=new double[4];
	}
	for(alpha=0;alpha<4;alpha++){
		u[alpha]=0.0;
		r[alpha]=0.0;
		for(beta=0;beta<4;beta++){
			pi[alpha][beta]=0.0;
		}
	}
	fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&tau,&r[1],&r[2],&eta,&Omega[0],&Omega[1],&Omega[2],&Omega[3],
		&u[1],&u[2],&u_eta,&eps_dec,&bdens_dec,&T_dec,&muB_dec,&muS_dec,&P_dec,
		&pi[0][0],&pi[0][1],&pi[0][2],&pi[0][3],&pi[1][1],&pi[1][2],&pi[1][3],&pi[2][2],&pi[2][3],&pi[3][3],&b);
	double udotpi[4],g[4]={1.0,-1.0,-1.0,-1.0};
	if(!feof(fptr)){
		T=1000.0*T_dec;
		r[0]=tau*cosh(eta);
		r[3]=tau*sinh(eta);
		u[3]=u_eta*tau;
		u[0]=sqrt(1.0+u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
			//printf("------ u=(%g,%g,%g,%g), tau/x/eta=(%g,%g,%g,%g), Omega=(%g,%g,%g,%g)\n",u[0],u[1],u[2],u[3],tau,r[1],r[2],eta,Omega[0],Omega[1],Omega[2],Omega[3]);
		pi[3][3]*=tau*tau;
		pi[2][3]*=tau;
		pi[1][3]*=tau;
		pi[0][3]*=tau;
		Omega[3]*=tau;
		pi[2][1]=pi[1][2]; pi[1][0]=pi[0][1]; pi[2][0]=pi[0][2];
		pi[3][1]=pi[1][3]; pi[3][2]=pi[2][3]; pi[3][0]=pi[0][3];
		for(alpha=0;alpha<4;alpha++){
			udotpi[alpha]=0.0;
			for(beta=0;beta<4;beta++){
				udotpi[alpha]+=g[beta]*u[beta]*pi[alpha][beta];
			}
		}
		Misc::BoostToCM(u,pi,dToverH);
		for(alpha=0;alpha<4;alpha++) delete [] pi[alpha];
			delete [] pi;

		for(alpha=0;alpha<4;alpha++){
			for(beta=0;beta<4;beta++){
				lambda[alpha][beta]=0.0;
			}
		}
		for(alpha=1;alpha<4;alpha++){
			for(beta=1;beta<4;beta++){
				lambda[alpha][beta]=dToverH[alpha][beta]/osuhydrotob3d->lambdafact;
			}
			lambda[alpha][alpha]+=1.0;
		}
		if(b3d->HYDRO_OCTANT_SYMMETRY==3) for(alpha=0;alpha<4;alpha++) Omega[alpha]*=8.0;
		if(b3d->HYDRO_OCTANT_SYMMETRY==2) for(alpha=0;alpha<4;alpha++) Omega[alpha]*=4.0;
		if(b3d->HYDRO_OCTANT_SYMMETRY==1) for(alpha=0;alpha<4;alpha++) Omega[alpha]*=2.0;
		//printf("HYDRO_OCTANT_SYMMETRY=%d\n",b3d->HYDRO_OCTANT_SYMMETRY);

		V2=Omega[0]*Omega[0]-Omega[1]*Omega[1]-Omega[2]*Omega[2]-Omega[3]*Omega[3];
		udotV=u[0]*Omega[0]-u[1]*Omega[1]-u[2]*Omega[2]-u[3]*Omega[3];
		V=sqrt(fabs(V2));
		if(V2>0){
			gamma=fabs(udotV/V);
			//printf("timelike Omega=(%g,%g,%g,%g), tau=%g\n",Omega[0],Omega[1],Omega[2],Omega[3],r[0]);
		}
		else{
			Vmag=sqrt(Omega[1]*Omega[1]+Omega[2]*Omega[2]+Omega[3]*Omega[3]);
			gamma=(u[0]*Vmag-(Omega[0]/Vmag)*(Omega[1]*u[1]+Omega[2]*u[2]+Omega[3]*u[3]))/V;
		//printf("spacelike Omega=(%g,%g,%g,%g), tau=%g, rvec.Omegavec=%g\n",Omega[0],Omega[1],Omega[2],Omega[3],r[0],Omega[1]*r[1]+Omega[2]*r[2]+Omega[3]*r[3]);
		}
		if(gamma<1.0){
			printf("Disaster: gamma =%g, but should be >1\n",gamma);
			printf("V2=%g, V=%g, gamma=%g, udotV=%g wmax=%g\n",V2,V,gamma,udotV,wmax);
			printf("Omega=(%g,%g,%g,%g)\n",Omega[0],Omega[1],Omega[2],Omega[3]);
			printf("fabs(MC_Nbar)=%g\n",fabs(osuhydrotob3d->MC_Nbar));
			exit(1);
		}

		wmax=V*(gamma+sqrt(gamma*gamma-1.0));
		//printf("V=%g, wmax=%g, Omega=(%g,%g,%g,%g), u=(%g,%g,%g,%g)\n",V,wmax,Omega[0],Omega[1],Omega[2],Omega[3],u[0],u[1],u[2],u[3]);
		// Note all of Josh's quantities were in mesh frame, we now boost to lab frame.
		return success;
	}
	else{
		return false;
	}
}

int COSUHydrotoB3D::MakeEvent(){
	meanpt=0.0;
	meanu=0.0;
	normpt=0;
	b3d->Reset();
	int nparts=0,iring;
	MC_NWrite=0;
	MC_Nbar=0.0;
	nparts=0;
	MC_Ntarget=-log(randy->ran());
	for(int icell=0;icell<NCELLS;icell++){
		if(b3d->HYDRO_PURE_BJORKEN)	nparts+=cell[icell].GenerateParticles();
		else nparts+=cell[icell].GenerateParticles3D();
	}
	printf("nparts=%d\n",nparts);
	return nparts;
}

int CFreezeoutCell::GenerateParticles(){
	CRandom *randy=osuhydrotob3d->randy;
	CResInfo *resinfo;
	double dNbarmax,rg1,rg2;
	int ncalls;
	int ires,iquad;
	int alpha,beta,nparts=0;
	double zsize,pdotV,pdotu,udotV;
	double ptilde[4],p[4],pt,y,et,weight,minv,mass,eta,tau,x[4];	 
	int ID,ipart;
	nparts=0;
		// For MC procedure, weights cannot exceed unity. These factors are added into densities and weights
		// to ensure MC weights do not exceed unity for any p. wmax is calculated so that the maximum weight
		// will be unity for a given u and V as p is varied. If wmax were larger, the answer would not change
		// but the efficiency would suffer. We check to make sure weights never exceed unity
	resinfo=b3d->reslist->GfirstResInfoptr;
	for(ires=0;ires<osuhydrotob3d->NRES;ires++){
		ID=resinfo->code;
		if(ID!=22){
			dNbarmax=osuhydrotob3d->density[ires]*wmax*double(b3d->NSAMPLE);
			//printf("ires=%d, dens=%g, dNbarmax=%g, \n",ires,osuhydrotob3d->density[ires],dNbarmax);
				//if(ID==211 || ID==-211 || ID==111) printf("dNBarmax=%g\n",dNbarmax);
			osuhydrotob3d->MC_Nbar+=dNbarmax;
			while(osuhydrotob3d->MC_Nbar>osuhydrotob3d->MC_Ntarget){
				osuhydrotob3d->MC_Ntarget-=log(randy->ran());
				mass=resinfo->mass;
				randy->generate_boltzmann(mass,T,ptilde);
				
				for(alpha=1;alpha<4;alpha++){
					p[alpha]=0.0;
					for(beta=1;beta<4;beta++) p[alpha]+=ptilde[beta]*lambda[alpha][beta];
						p[0]=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
				}

				for(alpha=0;alpha<4;alpha++) ptilde[alpha]=p[alpha];
					Misc::Boost(u,ptilde,p);
				pt=sqrt(p[1]*p[1]+p[2]*p[2]);

				pdotu=p[0]*u[0]-p[1]*u[1]-p[2]*u[2]-p[3]*u[3];
				pdotV=p[0]*Omega[0]-p[1]*Omega[1]-p[2]*Omega[2]-p[3]*Omega[3];
				if(pdotV>0.0){
					weight=(pdotV/pdotu)/wmax;
					//weight=1.0;
					//
					/**
					if(weight>0.99){
						printf("weight=%g\n",weight);
						Misc::Pause();
					}
					*/
					if(weight>1.00001){
						printf("Weight too large, weight=%g\n",weight);
						printf("Omega=(%g,%g,%g,%g)\n",Omega[0],Omega[1],Omega[2],Omega[3]);
						printf("p=(%g,%g,%g,%g), pdotV=%g\n",p[0],p[1],p[2],p[3],pdotV);
							//exit(1);
					}
					if(randy->ran()<weight){
						osuhydrotob3d->MC_NWrite+=1;
						nparts+=1;
						
						randy->gauss2(&rg1,&rg2);
						x[1]=r[1]+0.1*rg1;   // fix this, add spread below to x[0] and x[3] so that tau stays same
						x[2]=r[2]+0.1*rg2;

						x[0]=r[0];
						x[3]=r[3];

						eta=-b3d->ETAMAX+2.0*randy->ran()*b3d->ETAMAX;  // not necessary in 3D code
						tau=x[0];
						x[0]=tau*cosh(eta);
						x[3]=tau*sinh(eta);
						et=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]);
						y=0.5*log((p[0]+p[3])/(p[0]-p[3]));
						y+=eta;
						p[0]=et*cosh(y);
						p[3]=et*sinh(y);
							//
						if(fabs(sqrt(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3])-mass)>p[0]*1.0E-4){
							printf("invariant mass screwed up, =%g !=%g\n",sqrt(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3]),mass);
							printf("p=(%g,%g,%g,%g)\n",p[0],p[1],p[2],p[3]);
							exit(1);
						}

						if(ID==111 || ID==211 || ID==-211){
							osuhydrotob3d->meanpt+=sqrt(p[1]*p[1]+p[2]*p[2]);
							osuhydrotob3d->normpt+=1;
						}
						if(b3d->HYDRO_OCTANT_SYMMETRY==2){   // need to add possibility to reflect along z axis
							
							iquad=lrint(floor(4*randy->ran()));
							if(iquad==1 || iquad==2){
								x[1]=-x[1];
								p[1]=-p[1];
							}
							if(iquad==2 || iquad==3){
								x[2]=-x[2];
								p[2]=-p[2];
							}
							ipart=int(b3d->PartMap.size());
						}
						if(b3d->DeadPartMap.size()==0){
							printf("MUST INCREASE NPARTSMAX!!!!\n");
							exit(1);
						}
						b3d->partarray[ipart]->Init(ID,x[1],x[2],tau,eta,p[1],p[2],mass,y);
					}
				}
				//else{
					//printf("pdotV < 0, =%g\n",pdotV);
				//}
			}
		}
		resinfo=resinfo->nextResInfoptr;
	}
	//printf("created %d parts, wmax*netdens=%g, MC_Nbar=%g, MC_Ntarget=%g\n",nparts,wmax*osuhydrotob3d->netdens,osuhydrotob3d->MC_Nbar,osuhydrotob3d->MC_Ntarget);
	return nparts;
}

int CFreezeoutCell::GenerateParticles3D(){
	CRandom *randy=osuhydrotob3d->randy;
	CResInfo *resinfo;
	double dNbarmax,rg1,rg2;
	int ncalls;
	int ires,iquad;
	int alpha,beta,nparts=0;
	double zsize,pdotV,pdotu,udotV;
	double ptilde[4],p[4],pt,y,et,weight,minv,mass,eta,tau,x[4];	 
	int ID,ipart;
	nparts=0;
		// For MC procedure, weights cannot exceed unity. These factors are added into densities and weights
		// to ensure MC weights do not exceed unity for any p. wmax is calculated so that the maximum weight
		// will be unity for a given u and V as p is varied. If wmax were larger, the answer would not change
		// but the efficiency would suffer. We check to make sure weights never exceed unity
	resinfo=b3d->reslist->GfirstResInfoptr;
	for(ires=0;ires<osuhydrotob3d->NRES;ires++){
		ID=resinfo->code;
		if(ID!=22){
			dNbarmax=osuhydrotob3d->density[ires]*wmax*double(b3d->NSAMPLE);
				//printf("ires=%d, dens=%g, dNbarmax=%g, \n",ires,osuhydrotob3d->density[ires],dNbarmax);
				//if(ID==211 || ID==-211 || ID==111) printf("dNBarmax=%g\n",dNbarmax);
			osuhydrotob3d->MC_Nbar+=dNbarmax;
			while(osuhydrotob3d->MC_Nbar>osuhydrotob3d->MC_Ntarget){
				osuhydrotob3d->MC_Ntarget-=log(randy->ran());
				mass=resinfo->mass;
				randy->generate_boltzmann(mass,T,ptilde);
				
				for(alpha=1;alpha<4;alpha++){
					p[alpha]=0.0;
					for(beta=1;beta<4;beta++){
						p[alpha]+=ptilde[beta]*lambda[alpha][beta];
					}
				}
				p[0]=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);

				for(alpha=0;alpha<4;alpha++){
					ptilde[alpha]=p[alpha];
				}
				Misc::Boost(u,ptilde,p);
				pt=sqrt(p[1]*p[1]+p[2]*p[2]);

				pdotu=p[0]*u[0]-p[1]*u[1]-p[2]*u[2]-p[3]*u[3];
				pdotV=p[0]*Omega[0]-p[1]*Omega[1]-p[2]*Omega[2]-p[3]*Omega[3];
				//printf("weight=%g, Omega=(%g,%g,%g,%g), Omega^2=%g\n",pdotV/(pdotu*wmax),
				//Omega[0],Omega[1],Omega[2],Omega[3],Omega[0]*Omega[0]-Omega[1]*Omega[1]-Omega[2]*Omega[2]-Omega[3]*Omega[3]);
				if(pdotV>0.0){
					weight=(pdotV/pdotu)/wmax;
						//weight=1.0;
						//
					/**
					 if(weight>0.99){
					 printf("weight=%g\n",weight);
					 Misc::Pause();
					 }
					 */
					 if(weight>1.01){
					 	printf("Weight too large, weight=%g\n",weight);
					 	printf("Omega=(%g,%g,%g,%g)\n",Omega[0],Omega[1],Omega[2],Omega[3]);
					 	printf("p=(%g,%g,%g,%g), pdotV=%g\n",p[0],p[1],p[2],p[3],pdotV);
							//exit(1);
					 }
					 if(randy->ran()<weight){
					 	osuhydrotob3d->MC_NWrite+=1;
					 	nparts+=1;

					 	randy->gauss2(&rg1,&rg2);
					 	x[1]=r[1]+0.1*rg1;   
					 	x[2]=r[2]+0.1*rg2;

					 	tau=sqrt(r[0]*r[0]-r[3]*r[3]);				
					 	eta=asinh(r[3]/tau);
					 	randy->gauss2(&rg1,&rg2);
					 	eta+=0.1*rg1/tau;
					 	x[0]=tau*cosh(eta);
					 	x[3]=tau*sinh(eta);
					 	et=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]);
					 	y=0.5*log((p[0]+p[3])/(p[0]-p[3]));
					 	p[0]=et*cosh(y);
					 	p[3]=et*sinh(y);
							//
					 	if(fabs(sqrt(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3])-mass)>p[0]*1.0E-4){
					 		printf("invariant mass screwed up, =%g !=%g\n",sqrt(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3]),mass);
					 		printf("p=(%g,%g,%g,%g)\n",p[0],p[1],p[2],p[3]);
					 		exit(1);
					 	}

					 	if(ID==111 || ID==211 || ID==-211){
					 		osuhydrotob3d->meanpt+=sqrt(p[1]*p[1]+p[2]*p[2]);
					 		osuhydrotob3d->normpt+=1;
					 	}
					 	if(b3d->HYDRO_OCTANT_SYMMETRY==3){
					 		iquad=lrint(floor(8*randy->ran()));
					 		if(iquad==0 || iquad==1 || iquad==4 || iquad==5){
					 			x[1]=-x[1];
					 			p[1]=-p[1];
					 		}
					 		if(iquad==1 || iquad==2 || iquad==5 || iquad==6){
					 			x[2]=-x[2];
					 			p[2]=-p[2];
					 		}
					 		if(iquad==4 || iquad==5 || iquad==6 || iquad==7){
					 			eta=-eta;
					 			y=-y;
					 		}
					 	}
					 	if(b3d->HYDRO_OCTANT_SYMMETRY==2){
					 		iquad=lrint(floor(4*randy->ran()));
					 		if(iquad==0 || iquad==1){
					 			x[1]=-x[1];
					 			p[1]=-p[1];
					 		}
					 		if(iquad==1 || iquad==2){
					 			x[2]=-x[2];
					 			p[2]=-p[2];
					 			eta=-eta;
					 			y=-y;
					 		}
					 	}
					 	if(b3d->HYDRO_OCTANT_SYMMETRY==1){
					 		iquad=lrint(floor(2*randy->ran()));
					 		if(iquad==0){
					 			x[1]=-x[1];
					 			p[1]=-p[1];
					 		}
					 	}

					 	ipart=int(b3d->PartMap.size());
					 	if(b3d->DeadPartMap.size()==0){
					 		printf("MUST INCREASE NPARTSMAX!!!!\n");
					 		exit(1);
					 	}
							//printf("ID=%d, x=(%g,%g,%g,%g), pt=(%g,%g), y=%g\n",ID,tau,x[1],x[2],tau*eta,p[1],p[2],y);
					 	// Since Josh's calcs were in mesh frame, we need to boost by eta
					 	b3d->partarray[ipart]->Init(ID,x[1],x[2],tau,eta,p[1],p[2],mass,y+eta);
					}
				}
			}
		}
		resinfo=resinfo->nextResInfoptr;
	}
		//printf("created %d parts, netdens=%g, wmax*netdens=%g, MC_Nbar=%g, MC_Ntarget=%g\n",nparts,osuhydrotob3d->netdens,wmax*osuhydrotob3d->netdens,osuhydrotob3d->MC_Nbar,osuhydrotob3d->MC_Ntarget);
	return nparts;
}

void COSUHydrotoB3D::GetLambdaFact(){
	int ires,iQ,n,i;
	const int nmax=70;
	double G[nmax+5];
	double m,degen,z,Ipp=0.0,Ipptest=0.0,dIpp,Ptest=0.0,J,nfact,sign,alpha;
	double dIpptest=0.0,dp=4.0,p,e;
	CResInfo *resinfo=reslist->GfirstResInfoptr;
	for(ires=0;ires<NRES;ires++){
		if(resinfo->code!=22){
			m=resinfo->mass;
			degen=(2.0*resinfo->spin+1);
			z=m/T;
			alpha=0.0;

			G[0]=gsl_sf_gamma_inc(5,z)*pow(z,-5);
			for(int i=1;i<nmax+5;i++){
				n=5-2*i;
				if(n!=-1)	G[i]=(-exp(-z)/n)+(G[i-1]*z*z-z*exp(-z))/((n+1.0)*n);
				else G[i]=gsl_sf_gamma_inc(-1,z)*z;
			}
			J=0.0;
			nfact=1.0;
			sign=1.0;
			for(n=0;n<nmax;n+=1){
				if(n>0) sign=-1.0;
				J+=sign*nfact*(G[n]-2.0*G[n+1]+G[n+2]);
				nfact=nfact*0.5/(n+1.0);
				if(n>0) nfact*=(2.0*n-1.0);
			}
			dIpp=degen*exp(alpha)*pow(m,4)*(-z*J+15.0*gsl_sf_bessel_Kn(2,z)/(z*z));
			dIpp=dIpp/(60.0*PI*PI*HBARC*HBARC*HBARC);
			/*
			dIpptest=0.0;
			for(p=0.5*dp;p<3000;p+=dp){
				e=sqrt(m*m+p*p);
				dIpptest+=degen*(4.0*PI/(8.0*PI*PI*PI*HBARC*HBARC*HBARC))*p*p*dp*exp(-e/T)*( (2.0/3.0)*(p*p/e) - (2.0/15.0)*pow(p,4)/pow(e,3) );
				//dIpptest+=degen*(4.0*PI/(8.0*PI*PI*PI*HBARC*HBARC*HBARC))*p*p*dp*exp(-e/T)*((2.0/15.0)*pow(p,4)/pow(e,3) );
			}
			Ipptest+=dIpptest;
			*/
			
			Ipp+=dIpp;
			//Ptest+=Ipptest;
			//printf("dIpptest=%g =? %g, ratio=%g\n",dIpptest,dIpp,dIpptest/dIpp);
		}
		resinfo=resinfo->nextResInfoptr;
	}
	lambdafact=(2.0*P-4.0*Ipp)/(P+epsilon);
	//printf("Ipp=%g =? %g, lambdafact=%g\n",Ipptest,2.0*P-4.0*Ipp,lambdafact);
}

void COSUHydrotoB3D::TestLambdaFact(){
	int ires,i,j;
	double dp0=0.01,dp,pmag,e,bfact,phi,ctheta,stheta;
	double degen,m,a=0.1;
	double stress[4][4]={0.0},lamb[4][4]={0.0},pprime[4],p[4];
	lamb[1][1]=lamb[2][2]=lamb[3][3]=1.0;
	lamb[1][1]+=a;
	lamb[2][2]-=a;
	CResInfo *resinfo=reslist->GfirstResInfoptr;
	for(ires=0;ires<NRES;ires++){
		if(resinfo->code!=22){
			m=resinfo->mass;
			degen=(2.0*resinfo->spin+1);
			dp=sqrt(m*T)/10000.0;
			for(pmag=0.5*dp;pmag<20*sqrt(m*T+T*T);pmag+=dp){
				e=sqrt(m*m+pmag*pmag);
				phi=2.0*PI*randy->ran();
				ctheta=1.0-2.0*randy->ran();
				p[3]=pmag*ctheta;
				stheta=sqrt(1.0-ctheta*ctheta);
				p[1]=pmag*stheta*cos(phi);
				p[2]=pmag*stheta*sin(phi);
				bfact=(degen/(2.0*PI*PI*HBARC*HBARC*HBARC))*pmag*pmag*dp*exp(-e/T);
				pprime[1]=p[1]*(1.0+a);
				pprime[2]=p[2]*(1.0-a);
				pprime[3]=p[3];
				pprime[0]=sqrt(m*m+pprime[1]*pprime[1]+pprime[2]*pprime[2]+pprime[3]*pprime[3]);
				for(i=1;i<4;i++){
					for(j=1;j<4;j++){
						stress[i][j]+=bfact*((pprime[i]*pprime[j]/pprime[0])-(p[i]*p[j]/e));
					}
				}
			}
		}
		resinfo=resinfo->nextResInfoptr;
	}
	printf("a*lambdafact=%g, lambdafact=%g\n",a*lambdafact,lambdafact);
	for(i=1;i<4;i++){
		for(j=1;j<4;j++){
			printf("%10.3e ",stress[i][j]/(P+epsilon));
		}
		printf("\n");
	}
}

CFreezeoutCell::CFreezeoutCell(){
	dToverH=new double *[4];
	for(int alpha=0;alpha<4;alpha++) dToverH[alpha]=new double[4];
}

#endif
