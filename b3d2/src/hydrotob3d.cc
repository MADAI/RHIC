#ifndef __hydrotob3d_cc__
#define __hydrotob3d_cc__
#define __JOSH_FORMAT__

#include "hydrotob3d.h"

void CHYDROtoB3D::Init(){
	int ires,iring;
	string inputfilename;
	double pi,ei,sigma2i,dedti,degen;
	CResInfoMap::iterator rpos;
	CResInfo *resinfo;
	T=1000.0*parameter::getD(b3d->parmap,"HYDRO_FOTEMP",160.0);
	ETAMAX=b3d->ETAMAX;
	randy=b3d->randy;
	MC_Ntarget=0.0;
	MC_Nbar=0.0;
	MC_NWrite=0;
	nsample=b3d->NSAMPLE;
	reslist=b3d->reslist;
	nres=reslist->resmap.size();
#ifdef __HYDROTOB3D_XY_WRITE__
	xyfptr=fopen("xy.dat","w");
#endif
	//printf("epsilon_H=%g\n",intrinsic->epsilon);
	density.resize(nres);
	ring.resize(b3d->NRINGSMAX+1);
	epsilon=P=0.0;
	rpos=reslist->resmap.begin();
	for(ires=0;ires<nres;ires++){
		resinfo=rpos->second;
		if(resinfo->code!=22){
			degen=2.0*resinfo->spin+1.0;
			//printf("ires=%d, ID=%d, degen=%g, mass=%g, T=%g\n",ires,resinfo->code,degen,resinfo->mass,T);
			reslist->freegascalc_onespecies(resinfo->mass,T,ei,pi,density[ires],sigma2i,dedti);
			//printf("mass=%g, T=%g, pi=%g, ei=%g, di=%g\n",resinfo->mass,T,pi,ei,density[ires]);
			density[ires]*=degen;
			P+=degen*pi;
			epsilon+=degen*ei;
		}
		rpos++;
	}
	GetLambdaFact();
	//TestLambdaFact();
	initialization=true;
}

void CHYDROtoB3D::Init3D(){
	if(initialization==true){
		printf("CHYDROtoB3D already initialized\n");
		exit(1);
	}
	int ires,iring,ieta;
	string inputfilename;
	double pi,ei,sigma2i,dedti,degen;
	CResInfo *resinfo;
	CResInfoMap::iterator rpos;
	T=1000.0*parameter::getD(b3d->parmap,"HYDRO_FOTEMP",160.0);
	//ETAMAX=b3d->ETAMAX;
	NETA=parameter::getI(b3d->parmap,"HYDRO_NSIZE",1);
	DETA=parameter::getD(b3d->parmap,"HYDRO_DN",0.1);
	ETAMAX=NETA*DETA;
	randy=b3d->randy;
	MC_Ntarget=0.0;
	MC_Nbar=0.0;
	MC_NWrite=0;
	nsample=b3d->NSAMPLE;
	reslist=b3d->reslist;
	nres=reslist->resmap.size();
	//printf("epsilon_H=%g\n",intrinsic->epsilon);
	density.resize(nres);
	ring3d.resize(NETA);
	nrings3d.resize(NETA);
	for(ieta=0;ieta<NETA;ieta++){
		ring3d[ieta].resize(b3d->NRINGSMAX+1);
		nrings3d[ieta]=0;
	}
	epsilon=P=0.0;
	rpos=reslist->resmap.begin();
	for(ires=0;ires<nres;ires++){
		resinfo=rpos->second;
		if(resinfo->code!=22){
			degen=2.0*resinfo->spin+1.0;
			//printf("ires=%d, ID=%d, degen=%g, mass=%g, T=%g\n",ires,resinfo->code,degen,resinfo->mass,T);
			reslist->freegascalc_onespecies(resinfo->mass,T,ei,pi,density[ires],sigma2i,dedti);
			//printf("mass=%g, T=%g, pi=%g, ei=%g, di=%g\n",resinfo->mass,T,pi,ei,density[ires]);
			density[ires]*=degen;
			P+=degen*pi;
			epsilon+=degen*ei;
		}
		rpos++;
	}
	GetLambdaFact();
	//TestLambdaFact();
	initialization=true;
}

// ---------------------------------------------------------------

void CHYDROtoB3D::ReadInput(){
	if(!initialization){
		Init();
	}
	string inputfilename="output/"+b3d->run_name+"/"+b3d->qualifier+"/freezeout.dat";
	printf("freezeout info from %s\n",inputfilename.c_str());
	input=fopen(inputfilename.c_str(),"r");
	if(input != NULL){
		ReadHeader(input);
		nrings=0;
		do{
			ring[nrings].Read(input);
			
			if(ring[nrings].nread==0){
				ring[nrings].tau=ring[nrings-1].tau;
			}
			
			nrings++;
		}while(ring[nrings-1].nread!=0);
		fclose(input);
	}
}

void CHYDROtoB3D::ReadInput3D(){
	double tau,etaread;
	char dummy[120];
	string dummystring;
	int ieta,nread,iring;
	if(!initialization) Init3D();
	string inputfilename="output/"+b3d->run_name+"/"+b3d->qualifier+"/freezeout.dat";
	//	printf("freezeout info from %s\n",inputfilename.c_str());
	input=fopen(inputfilename.c_str(),"r");
	ReadHeader3D(input);
	fscanf(input,"%s",dummy);
	dummystring.assign(dummy);
	if(dummystring!="time:"){
		printf("YIKES, dummystring should have been time:, =%s\n",dummystring.c_str());
		exit(1);
	}
	fscanf(input,"%lf",&tau);
	for(ieta=0;ieta<NETA;ieta++){
		nrings3d[ieta]=0;
		//printf("nrings3d[%d]=%d\n",ieta,nrings3d[ieta]);
	}
	//printf("ReadInput3D, tau=%g\n",tau);
	while(!feof(input)){
		fscanf(input,"%s",dummy);
		if(!feof(input)){
			dummystring.assign(dummy);
			if(dummystring=="time:"){
				fscanf(input,"%lf",&tau);
				fscanf(input,"%s",dummy);
				dummystring.assign(dummy);
				fscanf(input,"%lf",&etaread);
			}
			else{
				if(dummystring!="eta:"){
					printf("dummystring should have been eta:, =%s, dummy=%s\n",dummystring.c_str(),dummy);
					exit(1);
				}
				fscanf(input,"%lf",&etaread);
			}
			fscanf(input,"%s",dummy);
			dummystring.assign(dummy);
			if(dummystring!="size:"){
				printf("YIKES, dummystring should have been size:\n");
				exit(1);
			}
			fscanf(input,"%d",&nread);		
			ieta=lrint(etaread/DETA);
			iring=nrings3d[ieta];
			ring3d[ieta][iring].tau=tau;
			ring3d[ieta][iring].etamin=etaread-0.5*DETA;
			if(ieta==0) ring3d[ieta][iring].etamin=0.0;
			ring3d[ieta][iring].etamax=etaread+0.5*DETA;
			ring3d[ieta][iring].nread=nread;
			printf("Read3D, tau=%g, etaread=%g, ieta=%d, iring=%d, nread=%d\n",tau,etaread,ieta,iring,nread);
			printf("---------------------------------------\n");
			
			ring3d[ieta][iring].Read3D(input);
			nrings3d[ieta]+=1;
			//fgets(dummy,100,input);
		}
	}
	fclose(input);
}

// -------------------------------------------------------------------

int CHYDROtoB3D::MakeEvent(){
	meanpt=0.0;
	normpt=0;
	etot=vtot=0.0;
	b3d->Reset();
	int nparts=0,iring;
	MC_NWrite=0;
	MC_Nbar=0.0;
	Ncheck=0.0;
	nparts=0;
	MC_Ntarget=-log(randy->ran());
	for(iring=1;iring<nrings;iring++){
		nparts+=GenerateParticles(iring);
	}
	
	//	printf("before b3d, nparts=%d, meanpt for pions=%g\n",nparts,meanpt/double(normpt));
	/** printf("nparts=%d, MC_Nbar=%g\n",nparts,MC_Nbar);
	printf("total energy should be %g, vtot=%g\n",0.5*etot/b3d->ETAMAX,0.5*vtot/b3d->ETAMAX); */
	printf("hydrotob3d made %d parts\n",nparts);
	return nparts;
}

int CHYDROtoB3D::MakeEvent3D(){
	meanpt=0.0;
	normpt=0;
	etot=vtot=0.0;
	b3d->Reset();
	int nparts=0,iring,ieta;
	MC_NWrite=0;
	MC_Nbar=0.0;
	Ncheck=0.0;
	nparts=0;
	MC_Ntarget=-log(randy->ran());
	for(ieta=0;ieta<NETA;ieta++){
		for(iring=1;iring<nrings3d[ieta];iring++){
			nparts+=GenerateParticles3D(ieta,iring);
		}
	}
	return nparts;
}

// -----------------------------------------------

int CHYDROtoB3D::GenerateParticles(int iring){
	CHB3DRing *earlier=&ring[iring-1];
	CHB3DRing *later=&ring[iring];
	CResInfo *resinfo;
	CResInfoMap::iterator rpos;
	double dNbarmax;
	int ncalls;
	int ires,iquad;
	int iphi1,iphi2,alpha,beta,nphi=later->nphi,nparts=0;
	double V0,Vx,Vy,V,V2,Vmag,x1[3],x2[3],y1[3],y2[3],a[3],b[3];
	double smallestr,biggestr,cphi1,sphi1,cphi2,sphi2,phi1,phi2,rx,ry,w,r;
	double zsize,lambda[4][4],u[4]={1.0,0.0,0.0,0.0},pdotV,pdotu,udotV,wmax,gamma;
	double ptilde[4],p[4],pt,y,et,weight,minv,mass,eta,tau,wx,wy,x[4];
	int ID,intweight;
	bool reality;
	CPart *newpart;
	for(iphi1=0;iphi1<nphi;iphi1++){
		iphi2=iphi1+1;
		phi1=0.5*PI*iphi1/double(nphi);
		phi2=0.5*PI*iphi2/double(nphi);
		cphi1=cos(phi1);
		sphi1=sin(phi1);
		cphi2=cos(phi2);
		sphi2=sin(phi2);
		x1[0]=later->tau; x1[1]=later->r[iphi1]*cphi1; x1[2]=later->r[iphi1]*sphi1;
		x2[0]=later->tau; x2[1]=later->r[iphi2]*cphi2; x2[2]=later->r[iphi2]*sphi2;
		y1[0]=earlier->tau; y1[1]=earlier->r[iphi1]*cphi1; y1[2]=earlier->r[iphi1]*sphi1;
		y2[0]=earlier->tau; y2[1]=earlier->r[iphi2]*cphi2; y2[2]=earlier->r[iphi2]*sphi2;
		for(alpha=1;alpha<4;alpha++){
			for(beta=1;beta<4;beta++)
				lambda[alpha][beta]=(0.25/lambdafact)*(later->dToverH[iphi1][alpha][beta]+later->dToverH[iphi2][alpha][beta]
					+earlier->dToverH[iphi1][alpha][beta]+earlier->dToverH[iphi2][alpha][beta]);
			lambda[alpha][alpha]+=1.0;
		}
		for(alpha=0;alpha<3;alpha++){
			a[alpha]=0.5*(x2[alpha]-x1[alpha]+y2[alpha]-y1[alpha]);
			b[alpha]=0.5*(x2[alpha]-y2[alpha]+x1[alpha]-y1[alpha]);
		}
		//a[0]=b[0]=0.0;
		zsize=ETAMAX*(earlier->tau+later->tau);
		u[1]=0.25*(earlier->ux[iphi1]+earlier->ux[iphi2]+later->ux[iphi1]+later->ux[iphi2]);
		u[2]=0.25*(earlier->uy[iphi1]+earlier->uy[iphi2]+later->uy[iphi1]+later->uy[iphi2]);
		u[0]=sqrt(1.0+u[1]*u[1]+u[2]*u[2]);
		V0=(a[1]*b[2]-b[1]*a[2])*zsize;
		Vx=-(a[2]*b[0]-b[2]*a[0])*zsize;
		Vy=-(a[0]*b[1]-b[0]*a[1])*zsize;
		V0*=4.0; Vx*=4.0; Vy*=4.0;
		V2=V0*V0-Vx*Vx-Vy*Vy;
		udotV=u[0]*V0-u[1]*Vx-u[2]*Vy;
		wmax=udotV+sqrt(-V2+udotV*udotV);
		
		if(gamma<1.0){
			printf("Disaster: gamma =%g, but should be >1\n",gamma);
			printf("x1=(%g,%g,%g), x2=(%g,%g,%g), y1=(%g,%g,%g), y2=(%g,%g,%g)\n",
			x1[0],x1[1],x1[2],x2[0],x2[1],x2[2],y1[0],y1[2],y1[2],y2[0],y2[1],y2[2]);
			printf("V0=%g, Vx=%g, Vy=%g\n",V0,Vx,Vy);
			printf("V2=%g, V=%g, gamma=%g, udotV=%g wmax=%g\n",V2,V,gamma,udotV,wmax);
			printf("fabs(MC_Nbar)=%g\n",fabs(MC_Nbar));
			exit(1);
		}
		rpos=reslist->resmap.begin();
		for(ires=0;ires<nres;ires++){
			resinfo=rpos->second;
			ID=resinfo->code;
			if(ID!=22){
				dNbarmax=density[ires]*wmax*double(nsample);
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
					pdotu=p[0]*u[0]-p[1]*u[1]-p[2]*u[2];
					weight=(pdotV/pdotu)/wmax;
					if(weight<0){
						if(b3d->COOPERFRYE_CREATENEGPARTS){
							intweight=-1;
							reality=false;
						}
						else{
							weight=0.0;
						}
					}
					else{
						intweight=1;
						reality=true;
					}
					if(fabs(weight)>1.00000001){
						printf("Weight too large, weight=%g\n",weight);
						printf("V0=%g, Vx=%g, Vy=%g\n",V0,Vx,Vy);
						printf("p=(%g,%g,%g,%g), V2=%g, pdotV=%g, V/udotV=%g\n",p[0],p[1],p[2],p[3],V2,pdotV,(V/udotV));
						//exit(1);
					}
					if(randy->ran()<fabs(weight)){
						MC_NWrite+=1;
						etot+=p[0];
						nparts+=1;
							
						rx=sqrt( pow(0.5*(x1[1]+x2[1]),2) + pow(0.5*(x1[2]+x2[2]),2) );
						ry=sqrt( pow(0.5*(y1[1]+y2[1]),2) + pow(0.5*(y1[2]+y2[2]),2) );

						w=randy->ran();
						r=sqrt(w*rx*rx+(1.0-w)*ry*ry);
						if((r>rx && r>ry) || (r<rx && r<ry)){
							printf("r is out of range in CHYDROtoB3D::GenerateParticles, = %g, rx=%g,ry=%g\n",r,rx,ry);
							exit(1);
						}
						wx=fabs(ry-r)/fabs(ry-rx);
						wy=1.0-wx;

						for(alpha=1;alpha<3;alpha++){
							x[alpha]=0.5*(wx*x1[alpha]+wx*x2[alpha]+wy*y1[alpha]+wy*y2[alpha]);
							x[alpha]+=(0.5-randy->ran())*a[alpha];
						}
						eta=ETAMAX-2.0*randy->ran()*ETAMAX;
						tau=wx*later->tau+wy*earlier->tau;
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
							
						iquad=lrint(floor(4*randy->ran()));
						if(iquad==1 || iquad==2){
							x[1]=-x[1];
							p[1]=-p[1];
						}
						if(iquad==2 || iquad==3){
							x[2]=-x[2];
							p[2]=-p[2];
						}
						newpart=b3d->GetDeadPart();
#ifdef __HYDROTOB3D_XY_WRITE__
						fprintf(xyfptr,"%g %g %g %g\n",tau,x[1],x[2],eta);
#endif
						newpart->Init(ID,x[1],x[2],tau,eta,p[1],p[2],mass,y,intweight,reality);
					}
				}
			}
			rpos++;
		}
	}
	return nparts;
}

int CHYDROtoB3D::GenerateParticles3D(int ieta,int iring){
	CHB3DRing *earlier=&ring3d[ieta][iring-1];
	CHB3DRing *later=&ring3d[ieta][iring];
	//printf("iring1=%d, tau1=%g, iring2=%d, tau2=%g\n",iring-1,earlier->tau,iring,later->tau);
	CResInfo *resinfo;
	CResInfoMap::iterator rpos;
	double dNbarmax,etamin,etamax;
	CPart *newpart;
	int ncalls;
	//MC_Ntarget-=log(randy->ran())/nsample;
	int ires,iquad;
	int iphi1,iphi2,alpha,beta,nphi=later->nphi,nparts=0;
	double V0,Vx,Vy,V,V2,Vmag,x1[3],x2[3],y1[3],y2[3],a[3],b[3];
	double smallestr,biggestr,cphi1,sphi1,cphi2,sphi2,phi1,phi2;
	double zsize,lambda[4][4],u[4]={1.0,0.0,0.0,0.0},pdotV,pdotu,udotV,wmax,gamma;
	double ptilde[4],p[4],pt,y,et,weight,minv,mass,eta,tau,wx,wy,x[4];	 
	int ID,intweight;
	bool reality;
	etamin=ring3d[ieta][iring].etamin;
	etamax=ring3d[ieta][iring].etamax;
	//printf("ieta=%d, iring=%d, etamin=%g, etamax=%g\n",ieta,iring,etamin,etamax);
	for(iphi1=0;iphi1<nphi;iphi1++){
		//printf("iphi1=%d, nphi=%d, ",iphi1,nphi);
		//printf("r1=%g r2=%g\n",earlier->r[iphi1],later->r[iphi1]);
		iphi2=iphi1+1;
		//printf("-----------------------grep : iphi1=%d, iphi2=%d, nphi=%d\n",iphi1,iphi2,nphi);
		phi1=0.5*PI*iphi1/double(nphi);
		phi2=0.5*PI*iphi2/double(nphi);
		cphi1=cos(phi1);
		sphi1=sin(phi1);
		cphi2=cos(phi2);
		sphi2=sin(phi2);
		x1[0]=later->tau; x1[1]=later->r[iphi1]*cphi1; x1[2]=later->r[iphi1]*sphi1;
		x2[0]=later->tau; x2[1]=later->r[iphi2]*cphi2; x2[2]=later->r[iphi2]*sphi2;
		y1[0]=earlier->tau; y1[1]=earlier->r[iphi1]*cphi1; y1[2]=earlier->r[iphi1]*sphi1;
		y2[0]=earlier->tau; y2[1]=earlier->r[iphi2]*cphi2; y2[2]=earlier->r[iphi2]*sphi2;
		for(alpha=1;alpha<4;alpha++){
			for(beta=1;beta<4;beta++)
				lambda[alpha][beta]=(0.25/lambdafact)*(later->dToverH[iphi1][alpha][beta]+later->dToverH[iphi2][alpha][beta]
					+earlier->dToverH[iphi1][alpha][beta]+earlier->dToverH[iphi2][alpha][beta]);
			//lambda[alpha][beta]=0.0;
			lambda[alpha][alpha]+=1.0;
		}
		for(alpha=0;alpha<3;alpha++){
			a[alpha]=0.5*(x2[alpha]-x1[alpha]+y2[alpha]-y1[alpha]);
			b[alpha]=0.5*(x2[alpha]-y2[alpha]+x1[alpha]-y1[alpha]);
		}
		//a[0]=b[0]=0.0;
		zsize=fabs(etamax-etamin)*(earlier->tau+later->tau);
		u[1]=0.25*(earlier->ux[iphi1]+earlier->ux[iphi2]+later->ux[iphi1]+later->ux[iphi2]);
		u[2]=0.25*(earlier->uy[iphi1]+earlier->uy[iphi2]+later->uy[iphi1]+later->uy[iphi2]);
		u[0]=sqrt(1.0+u[1]*u[1]+u[2]*u[2]);
		V0=(a[1]*b[2]-b[1]*a[2])*zsize;
		Vx=-(a[2]*b[0]-b[2]*a[0])*zsize;
		Vy=-(a[0]*b[1]-b[0]*a[1])*zsize;
		V0*=4.0; Vx*=4.0; Vy*=4.0;
		V2=V0*V0-Vx*Vx-Vy*Vy;
		udotV=u[0]*V0-u[1]*Vx-u[2]*Vy;
		V=sqrt(fabs(V2));
		//etot+=V*(304.0*u[0]*u[0]+50.3*(u[0]*u[0]-1.0));
		vtot+=V;
		//printf("phi1=%g, phi2=%g, V=%g\n",phi1*180/PI,phi2*180/PI,V);
		//printf("iphi1=%d, V=%g=(%g,%g,%g), u=(%g,%g,%g)\n",iphi1,V,V0,Vx,Vy,u[0],u[1],u[2]);
		//printf("a=(%g,%g,%g), b=(%g,%g,%g)\n",a[0],a[1],a[2],b[0],b[1],b[2]);
		//printf("x1=(%g,%g,%g)\n",x1[0],x1[1],x1[2]);
		//printf("y1=(%g,%g,%g)\n",y1[0],y1[1],y1[2]);
		//printf("_________________________________\n");
		if(V2>0){
			gamma=fabs(udotV/V);
		}
		else{
			Vmag=sqrt(Vx*Vx+Vy*Vy);
			gamma=(u[0]*Vmag-(V0/Vmag)*(Vx*u[1]+Vy*u[2]))/V;
		}
		wmax=V*(gamma+sqrt(gamma*gamma-1.0));
		//printf("iphi1=%d, wmax=%g\n",iphi1,wmax);
		// For MC procedure, weights cannot exceed unity. These factors are added into densities and weights
		// to ensure MC weights do not exceed unity for any p. wmax is calculated so that the maximum weight
		// will be unity for a given u and V as p is varied. If wmax were larger, the answer would not change
		// but the efficiency would suffer. We check to make sure weights never exceed unity
		if(gamma<1.0){
			printf("Disaster: gamma =%g, but should be >1\n",gamma);
			printf("x1=(%g,%g,%g), x2=(%g,%g,%g), y1=(%g,%g,%g), y2=(%g,%g,%g)\n",
			x1[0],x1[1],x1[2],x2[0],x2[1],x2[2],y1[0],y1[2],y1[2],y2[0],y2[1],y2[2]);
			printf("V0=%g, Vx=%g, Vy=%g\n",V0,Vx,Vy);
			printf("V2=%g, V=%g, gamma=%g, udotV=%g wmax=%g\n",V2,V,gamma,udotV,wmax);
			printf("fabs(MC_Nbar)=%g\n",fabs(MC_Nbar));
			exit(1);
		}
		rpos=reslist->resmap.begin();
		for(ires=0;ires<nres;ires++){
			resinfo=rpos->second;
			ID=resinfo->code;
			if(ID!=22){
				dNbarmax=density[ires]*wmax*double(nsample);
				//if(ID==211 || ID==-211 || ID==111) printf("dNBarmax=%g\n",dNbarmax);
				MC_Nbar+=dNbarmax;
				while(MC_Nbar>MC_Ntarget){
					MC_Ntarget-=log(randy->ran());
					//MC_Ntarget+=1.0;
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
					//if(pdotV>0.0){
					pdotu=p[0]*u[0]-p[1]*u[1]-p[2]*u[2];
					weight=(pdotV/pdotu)/wmax;
					if(weight<0){
						if(b3d->COOPERFRYE_CREATENEGPARTS){
							intweight=-1;
							reality=false;
						}
						else{
							weight=0.0;
						}
					}
					else{
						intweight=1;
						reality=true;
					}
					if(fabs(weight)>1.00000001){
						printf("Weight too large, weight=%g\n",weight);
						printf("V0=%g, Vx=%g, Vy=%g\n",V0,Vx,Vy);
						printf("p=(%g,%g,%g,%g), V2=%g, pdotV=%g, V/udotV=%g\n",p[0],p[1],p[2],p[3],V2,pdotV,(V/udotV));
						//exit(1);
					}
					if(randy->ran()<fabs(weight)){
						MC_NWrite+=1;
						etot+=p[0];
						nparts+=1;
						wx=randy->ran();
						wy=1.0-wx;
						for(alpha=1;alpha<3;alpha++){
							x[alpha]=0.5*(wx*x1[alpha]+wx*x2[alpha]+wy*y1[alpha]+wy*y2[alpha]);
							x[alpha]+=(0.5-randy->ran())*a[alpha];
						}
						eta=etamin+randy->ran()*DETA;
						if(randy->ran()<0.5) eta=-eta;
						tau=wx*later->tau+wy*earlier->tau;
						x[0]=tau*cosh(eta);
						x[3]=tau*sinh(eta);
						//printf("success: p=(%g,%g,%g,%g), x=(%g,%g,%g,%g), tau=%g, eta=%g\n",p[0],p[1],p[2],p[3],x[0],x[1],x[2],x[3],tau,eta);
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
							meanpt+=intweight*sqrt(p[1]*p[1]+p[2]*p[2]);
							normpt+=intweight;
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
						newpart=b3d->GetDeadPart();
						newpart->Init(ID,x[1],x[2],tau,eta,p[1],p[2],mass,y,intweight,reality);
					}
				}
			}
			rpos++;
		}
	}
	return nparts;
}

// -----------------------------------------------

void CHYDROtoB3D::GetLambdaFact(){
	int iQ,n,i;
	const int nmax=70;
	double G[nmax+5];
	double m,degen,z,Ipp=0.0,Ipptest=0.0,dIpp,Ptest=0.0,J,nfact,sign,alpha;
	double dIpptest=0.0,dp=4.0,p,e;
	CResInfo *resinfo;
	CResInfoMap::iterator rpos;
	for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();rpos++){
		resinfo=rpos->second;
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
	}
	lambdafact=(2.0*P-4.0*Ipp)/(P+epsilon);
	//printf("P=%g, epsilon=%g\n",P,epsilon);
	//printf("Ipp=%g =? %g, lambdafact=%g\n",Ipptest,2.0*P-4.0*Ipp,lambdafact);
}

void CHYDROtoB3D::TestLambdaFact(){
	int i,j;
	double dp0=0.01,dp,pmag,e,bfact,phi,ctheta,stheta;
	double degen,m,a=0.1;
	double stress[4][4]={0.0},lamb[4][4]={0.0},pprime[4],p[4];
	lamb[1][1]=lamb[2][2]=lamb[3][3]=1.0;
	lamb[1][1]+=a;
	lamb[2][2]-=a;
	CResInfo *resinfo;
	CResInfoMap::iterator rpos;
	for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->code!=22){
			m=resinfo->mass;
			printf("m=%g\n",m);
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
	}
	printf("a*lambdafact=%g, lambdafact=%g\n",a*lambdafact,lambdafact);
	for(i=1;i<4;i++){
		for(j=1;j<4;j++){
			printf("%10.3e ",stress[i][j]/(P+epsilon));
		}
		printf("\n");
	}
}

void CHYDROtoB3D::ReadHeader(FILE *fptr){
	char dummy[180];
	string dummystring;
	do{
		fgets(dummy,120,fptr);
		dummystring.assign(dummy);
		//printf("dummy header %s",dummy);
	}while(dummystring!="END_OF_HEADER\n");
	fscanf(fptr,"%s",dummy);
	dummystring.assign(dummy);
	if(dummystring!="time:"){
		printf("bad header\n");
		exit(1);
	}
}

void CHYDROtoB3D::ReadHeader3D(FILE *fptr){
	char dummy[180];
	string dummystring;
	do{
		fgets(dummy,120,fptr);
		dummystring.assign(dummy);
	}while(dummystring!="END_OF_HEADER\n");
}

CHB3DRing::CHB3DRing(){
	int iphi,alpha,beta;
	for(iphi=0;iphi<=nphi;iphi++){
		for(alpha=0;alpha<4;alpha++){
			for(beta=0;beta<4;beta++){
				dToverH[iphi][alpha][beta]=0.0;
			}
		}
	}
}

void CHB3DRing::Read(FILE *fptr){
	int i;
	string dummystring;
	char dummy[120];
	const int NPTS=300;
	double x[NPTS],y[NPTS],uxi[NPTS],uyi[NPTS];
	//double **PiOverh,**PiOverhtilde,Pixxoverh,Pixyoverh,Piyyoverh;
	double dToverH_xx[NPTS],dToverH_xy[NPTS],dToverH_yy[NPTS];
	int i1,i2,iphi,alpha,beta;
	double nt,nx,ny,a1,a2,a3,a4,a5,b=0.0,epsilon,P,Rqgp,T;
	const double root3=sqrt(3.0);
	i=-1;	
	fscanf(fptr,"%lf",&tau);
	fscanf(fptr,"%lf",&x[0]);
	do{
		i+=1;
		if(i==300){
			printf("CHB3DRing::Read, arrays too small\n");
			exit(1);
		}
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &y[i],&epsilon,&P,&T,&Rqgp,&uxi[i],&uyi[i],&nt,&nx,&ny,&a1,&a2,&a3,&a4,&a5);
		//printf("&&&& x[%d]=%g, y[%d]=%g, tau=%g\n",i,x[i],i,y[i],tau);
		//a1=a2=a3=a4=a4=0.0;
		if(!feof(fptr)){
			//if(i==0) printf("reading along y axis, tau=%g, y=%g, u=(%g,%g)\n",tau,y[0],uxi[0],uyi[0]);
			fgets(dummy,100,fptr);
			//FillOutPi(PiOverh,dToverH_xx[i],dToverH_xy[i],dToverH_yy[i],uxi[i],uyi[i]);
			dToverH_xx[i]=(b+a1+a2/root3);
			dToverH_xy[i]=a3;
			dToverH_yy[i]=(b-a1+a2/root3);
			fscanf(fptr,"%s",dummy);
			dummystring.assign(dummy);
			//printf("x,y=(%g,%g)\n",x[i],y[i]);
			if(dummystring!="time:") x[i+1]=atof(dummy);
		}
	}while(dummystring!="time:" && feof(fptr)==false);
	nread=i;
	FillOutArrays(x,y,uxi,uyi,dToverH_xx,dToverH_xy,dToverH_yy);
	//if(nread>0) printf("TIME=%g, %d points for Cooper-Frye ring\n",tau,nread);
}

void CHB3DRing::Read3D(FILE *fptr){
	string dummystring;
	char dummy[120];
	const int NPTS=300;
	double x[NPTS],y[NPTS],uxi[NPTS],uyi[NPTS];
	//double **PiOverh,**PiOverhtilde,Pixxoverh,Pixyoverh,Piyyoverh;
	double dToverH_xx[NPTS],dToverH_xy[NPTS],dToverH_yy[NPTS];
	int i,i1,i2,iphi,alpha,beta,ieta;
	double nt,nx,ny,a1,a2,a3,a4,a5,b=0.0,epsilon,P,Rqgp,T;
	const double root3=sqrt(3.0);
	if(nread>300){
		printf("CHB3DRing::Read, arrays too small\n");
		exit(1);
	}
	for(i=0;i<nread;i++){
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &x[i],&y[i],&epsilon,&P,&T,&Rqgp,&uxi[i],&uyi[i],&nt,&nx,&ny,&a1,&a2,&a3,&a4,&a5);
		fgets(dummy,100,fptr);
		dToverH_xx[i]=(b+a1+a2/root3);
		dToverH_xy[i]=a3;
		dToverH_yy[i]=(b-a1+a2/root3);
	}
	FillOutArrays(x,y,uxi,uyi,dToverH_xx,dToverH_xy,dToverH_yy);
}

void CHB3DRing::FillOutArrays(double *x,double *y,double *uxi,double *uyi,double *dToverH_xx,double *dToverH_xy,double *dToverH_yy){
	//char filename[120];
	//FILE *ringinfo;
	double phi,cphi,sphi,dphi,dphi1,dphi2,dphi1_min,dphi2_min,xx,yy,r1,r2,w1,w2;
	int i,i1,i2,iphi,alpha,beta;
	if(nread>0){
		//sprintf(filename,"plotdata/ring_tau%g.dat",tau);
		//ringinfo=fopen(filename,"w");
		for(iphi=0;iphi<=nphi;iphi++){
			//printf("x[%d]=%g, y[%d]=%g\n",iphi,x[iphi],iphi,y[iphi]);
			phi=iphi*0.5*PI/nphi;
			cphi=cos(phi);
			sphi=sin(phi);
			i1=i2=0;
			dphi1=10000*PI;
			dphi2=-10000*PI;
			for(i=0;i<nread;i++){
				xx=x[i]*cphi+y[i]*sphi;
				yy=x[i]*sphi-y[i]*cphi;
				dphi=atan2(yy,xx);
				if(dphi<0.0) dphi+=2.0*PI;
				if(dphi<dphi1){
					dphi1=dphi;
					i1=i;
				}
				dphi-=2.0*PI;
				if(dphi>dphi2){
					dphi2=dphi;
					i2=i;
				}
			}
			w1=fabs(dphi2)/(dphi1-dphi2);
			w2=1.0-w1;
			r1=sqrt(x[i1]*x[i1]+y[i1]*y[i1]);
			r2=sqrt(x[i2]*x[i2]+y[i2]*y[i2]);
			r[iphi]=w1*r1+w2*r2;
			ux[iphi]=w1*uxi[i1]+w2*uxi[i2];
			uy[iphi]=w1*uyi[i1]+w2*uyi[i2];
			dToverH[iphi][1][1]=w1*dToverH_xx[i1]+w2*dToverH_xx[i2];
			dToverH[iphi][1][2]=w1*dToverH_xy[i1]+w2*dToverH_xy[i2];
			dToverH[iphi][2][1]=dToverH[iphi][1][2];
			dToverH[iphi][2][2]=w1*dToverH_yy[i1]+w2*dToverH_yy[i2];
			dToverH[iphi][3][3]=-dToverH[iphi][2][2]-dToverH[iphi][1][1];
		}
	}
	else{
		for(iphi=0;iphi<=nphi;iphi++){
			ux[iphi]=uy[iphi]=0.0;
			r[iphi]=0.0;
			dToverH[iphi][1][1]=0.0;
			dToverH[iphi][1][2]=dToverH[iphi][2][1]=0.0;
			dToverH[iphi][2][2]=0.0;
			dToverH[iphi][3][3]=0.0;
		}
	}
}

#endif
