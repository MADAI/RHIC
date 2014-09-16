#ifndef __hydrotob3d_cc__
#define __hydrotob3d_cc__
#define __JOSH_FORMAT__

#include "hydrotob3d.h"
using namespace std;
//initialize static variables
double CPRCell::T=0.0;
double CPRCell::P=0.0;
double CPRCell::epsilon=0.0;
double CPRCell::omega_xyspacing=0.0;
double CPRCell::omega_tauspacing=0.0;
double CPRCell::etamax=0.0;
CPRCell::CPRCell(){
	ux=uy=x=y=0.0;
}

void CHYDROtoB3D::InitPR(){
	if(initialization==true){
		printf("CHYDROtoB3D already initialized\n");
		exit(1);
	}
	int ires,iring;
	string inputfilename;
	double pi,ei,sigma2i,dedti,degen;
	T=parameter::getD(b3d->parmap,"HYDRO_FOTEMP",150.0);
	ETAMAX=b3d->ETAMAX;
	randy=b3d->randy;
	MC_Ntarget=0.0;
	MC_Nbar=0.0;
	MC_NWrite=0;
	nsample=b3d->NSAMPLE;
	reslist=b3d->reslist;
	nres=reslist->NResonances;
	//printf("epsilon_H=%g\n",intrinsic->epsilon);
	density=new double[nres];
	prcell=new CPRCell[b3d->NPRCELLSMAX];
	epsilon=P=0.0;
	CResInfo *resinfo=reslist->GfirstResInfoptr;
	for(ires=0;ires<nres;ires++){
		if(resinfo->code!=22){
			degen=2.0*resinfo->spin+1.0;
			printf("ires=%d, ID=%d, degen=%g, mass=%g, T=%g\n",ires,resinfo->code,degen,resinfo->mass,T);
			b3d->freegascalc_onespecies(resinfo->mass,T,pi,ei,density[ires],sigma2i,dedti);
			density[ires]*=degen;
			//printf("mass=%g, T=%g, pi=%g, ei=%g, di=%g\n",resinfo->mass,T,pi,ei,density[ires]);
			P+=degen*pi;
			epsilon+=degen*ei;
		}
		resinfo=resinfo->nextResInfoptr;
	}
	GetLambdaFact();
	//TestLambdaFact();
	CPRCell::T=T;
	CPRCell::epsilon=epsilon;
	CPRCell::P=P;
	CPRCell::etamax=b3d->ETAMAX;
	CPRCell::omega_xyspacing=parameter::getD(b3d->parmap,"B3D_PR_OMEGA_XYSPACING",0.0);
	CPRCell::omega_tauspacing=parameter::getD(b3d->parmap,"B3D_PR_OMEGA_TAUSPACING",0.0);
	initialization=true;
}

void CHYDROtoB3D::ReadInputPR(){
	if(!initialization) InitPR();
	string inputfilename="output/"+b3d->run_name+"/"+b3d->qualifier+"/freezeout.dat";
	printf("freezeout info from %s\n",inputfilename.c_str());
	input=fopen(inputfilename.c_str(),"r");
	if(input != NULL){
		nprcells=0;
		do{
			prcell[nprcells].Read(input);
			nprcells+=1;
		}while(!feof(input));
		fclose(input);
	}
}

int CHYDROtoB3D::MakeEventPR(){
	meanpt=0.0;
	meanu=0.0;
	normpt=0;
	etot=vtot=0.0;
	b3d->Reset();
	int nparts=0,iprcell;
	MC_NWrite=0;
	MC_Nbar=0.0;
	Ncheck=0.0;
	nparts=0;
	MC_Ntarget=-log(randy->ran());
	for(iprcell=0;iprcell<nprcells;iprcell++){
		nparts+=GenerateParticlesPR(iprcell);
	}
	printf("hydrotob3dPR made %d parts\n",nparts);
	return nparts;
}

int CHYDROtoB3D::GenerateParticlesPR(int iprcell){
	CResInfo *resinfo;
	CPRCell *cell=&prcell[iprcell];
	double dNbarmax,gspread=0.2;
	int ncalls,ID,ipart;
	int ires,iquad,alpha,beta;
	int nparts=0;
	double V,V2,V0,Vx,Vy,Vmag,g1,g2;
	double lambda[4][4],u[4]={1.0,0.0,0.0,0.0},pdotV,pdotu,udotV,wmax,gamma;
	double ptilde[4],p[4],pt,y,et,weight,minv,mass,eta,tau,wx,wy,x[4];
	//printf("deadpartmap.size=%d, partmap.size=%d, finalpartmap.size=%d,NPARTSMAX=%d\n",
	//	int(b3d->DeadPartMap.size()),int(b3d->PartMap.size()),int(b3d->FinalPartMap.size()),int(b3d->NPARTSMAX));
	randy->gauss2(&g1,&g2);
	x[1]=cell->x+gspread*g1;
	x[2]=cell->y+gspread*g2;
	x[0]=cell->tau;
	tau=cell->tau;
	u[1]=cell->ux;
	u[2]=cell->uy;
	u[0]=sqrt(1.0+u[1]*u[1]+u[2]*u[2]);
	eta=ETAMAX-2.0*randy->ran()*ETAMAX;

	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++)
			lambda[alpha][beta]=(cell->dToverH[alpha][beta])/lambdafact;
		lambda[alpha][alpha]+=1.0;
	}
	V0=cell->Omega_0;
	Vx=-cell->Omega_x;
	Vy=-cell->Omega_y;
	V2=V0*V0-Vx*Vx-Vy*Vy;
	udotV=u[0]*V0-u[1]*Vx-u[2]*Vy;
	V=sqrt(fabs(V2));
	if(V2>0){
		gamma=fabs(udotV/V);
	}
	else{
		Vmag=sqrt(Vx*Vx+Vy*Vy);
		gamma=(u[0]*Vmag-(V0/Vmag)*(Vx*u[1]+Vy*u[2]))/V;
	}
	wmax=V*(gamma+sqrt(gamma*gamma-1.0));
	vtot+=wmax;
	//printf("Vmag=%g, V=%g, wmax=%g, Omega=(%g,%g,%g), vtot=%g\n",Vmag,V, wmax,V0,Vx,Vy,vtot);
	if(gamma<1.0){
		printf("Disaster: gamma =%g, but should be >1\n",gamma);
		printf("V0=%g, Vx=%g, Vy=%g\n",V0,Vx,Vy);
		printf("V2=%g, V=%g, gamma=%g, udotV=%g wmax=%g\n",V2,V,gamma,udotV,wmax);
		printf("fabs(MC_Nbar)=%g\n",fabs(MC_Nbar));
		exit(1);
	}
	resinfo=reslist->GfirstResInfoptr;
	for(ires=0;ires<nres;ires++){
		ID=resinfo->code;
		if(ID!=22){
			dNbarmax=density[ires]*wmax*double(nsample);
			//printf("dNbarmax=%g, density[%d]=%g\n",dNbarmax,ires,density[ires]);
			MC_Nbar+=dNbarmax;
			while(MC_Nbar>MC_Ntarget){
				MC_Ntarget-=log(randy->ran());
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

				pdotu=p[0]*u[0]-p[1]*u[1]-p[2]*u[2];
				pdotV=p[0]*V0-p[1]*Vx-p[2]*Vy;
				if(pdotV>0.0){
					pdotu=p[0]*u[0]-p[1]*u[1]-p[2]*u[2];
					weight=(pdotV/pdotu)/wmax;
					if(weight>1.00000001){
						printf("Weight too large, weight=%g\n",weight);
						printf("V0=%g, Vx=%g, Vy=%g\n",V0,Vx,Vy);
						printf("p=(%g,%g,%g,%g), V2=%g, pdotV=%g, V/udotV=%g\n",p[0],p[1],p[2],p[3],V2,pdotV,(V/udotV));
					}
					if(randy->ran()<weight){
						MC_NWrite+=1;
						etot+=p[0];
						nparts+=1;


						x[0]=tau*cosh(eta);
						x[3]=tau*sinh(eta);
						et=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]);
						y=0.5*log((p[0]+p[3])/(p[0]-p[3]));
						y+=eta;
						p[0]=et*cosh(y);
						p[3]=et*sinh(y);
						if(fabs(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3]-mass*mass)>1.0E-2){
							printf("invariant mass screwed up, =%g !=%g\n",sqrt(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3]),mass);
							printf("p=(%g,%g,%g,%g)\n",p[0],p[1],p[2],p[3]);
							exit(1);
						}

						if(ID==111 || ID==211 || ID==-211){
							meanpt+=sqrt(p[1]*p[1]+p[2]*p[2]);
							meanu+=(p[1]*u[1]+p[2]*u[2])/pt;
							normpt+=1;
						}

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
						if(b3d->DeadPartMap.size()==0){
							printf("MUST INCREASE NPARTSMAX!!!! NPARTSMAX=%d\n",b3d->NPARTSMAX);
							exit(1);
						}
						b3d->partarray[ipart]->Init(ID,x[1],x[2],tau,eta,p[1],p[2],mass,y);
					}
				}
			}
		}
		resinfo=resinfo->nextResInfoptr;
	}
	return nparts;
}

void CPRCell::Read(FILE *fptr){
	double pixx_overh,pixy_overh,piyy_overh,zsize;
	int idirection;
	fscanf(fptr,"%lf %lf %lf %d %lf %lf %lf %lf %lf %lf",
		&x,&y,&tau,&idirection,&ux,&uy,&pixx_overh,&pixy_overh,&piyy_overh,&T);
	Omega_x=Omega_y=Omega_0=0.0;
	zsize=2.0*etamax*tau;
	//printf("zsize=%g, etamax=%g, tau=%g\n",zsize,etamax,tau);
	//printf("omega_xyspacing+%g, omega_tauspacing=%g\n",omega_xyspacing,omega_tauspacing);
	//printf("idirection=%d\n",idirection);
	if(idirection==1){
		Omega_x=omega_xyspacing*omega_tauspacing*zsize;
	}
	else if(idirection==-1){
		Omega_x=-omega_xyspacing*omega_tauspacing*zsize;
	}
	else if(idirection==2){
		Omega_y=omega_xyspacing*omega_tauspacing*zsize;
	}
	else if(idirection==-2){
		Omega_y=-omega_xyspacing*omega_tauspacing*zsize;
	}
	else if(idirection==3){
		Omega_0=omega_xyspacing*omega_xyspacing*zsize;
	}
	else if(idirection==-3){
		Omega_0=-omega_xyspacing*omega_xyspacing*zsize;
	}
	//printf("Omega=(%g,%g,%g)\n",Omega_0,Omega_x,Omega_y);
	GetPiTildeOverH(pixx_overh,pixy_overh,piyy_overh);
}

void CPRCell::GetPiTildeOverH(double pixx,double pixy,double piyy){ // This gets Pi/h in fluid frame from pixx/h,pixy/h,piyy/h, ux,uy
	// then boosts to get PiTilde
	double Pi[4][4],u[4];
	double u0=sqrt(1.0+ux*ux+uy*uy);
	u[0]=u0;
	u[1]=ux;
	u[2]=uy;
	u[3]=0.0;
	Pi[0][0]=(ux*ux*pixx+2*ux*uy*pixy+uy*uy*piyy)/(u0*u0);
	Pi[0][1]=(ux*pixx+uy*pixy)/u0;
	Pi[0][2]=(ux*pixy+uy*piyy)/u0;
	Pi[0][3]=0.0;
	Pi[1][0]=Pi[0][1];
	Pi[1][1]=pixx;
	Pi[1][2]=pixy;
	Pi[1][3]=0.0;
	Pi[2][0]=Pi[0][2];
	Pi[2][1]=Pi[1][2];
	Pi[2][2]=piyy;
	Pi[2][3]=0.0;
	Pi[3][0]=0.0;
	Pi[3][1]=0.0;
	Pi[3][2]=0.0;
	Pi[3][3]=Pi[0][0]-Pi[1][1]-Pi[2][2];
	//printf("Trace Pi=%g =? 0\n",Pi[0][0]-Pi[1][1]-Pi[2][2]-Pi[3][3]);
	// Now boost to com frame
	int alpha,beta,gamma,delta;
	double Linv[4][4],L[4][4];  // Lorentz Boost Matrices
	double ucontra[4]={u[0],-u[1],-u[2],-u[3]},n[4]={1.0,0,0,0};
	for(alpha=0;alpha<4;alpha++){
		for(beta=0;beta<4;beta++){
			Linv[alpha][beta]=2.0*n[alpha]*ucontra[beta]-(u[alpha]+n[alpha])*(ucontra[beta]+n[beta])/(1.0+u[0]);
			if(alpha==beta) Linv[alpha][beta]+=1.0;
			L[beta][alpha]=Linv[alpha][beta];
		}
	}
	for(alpha=0;alpha<4;alpha++){
		for(delta=0;delta<4;delta++){
			dToverH[alpha][delta]=0.0;
			for(beta=0;beta<4;beta++){
				for(gamma=0;gamma<4;gamma++){
					dToverH[alpha][delta]+=Linv[alpha][beta]*Pi[beta][gamma]*L[gamma][delta];
				}
			}
		}
	}
}


#endif
