#ifndef __JOSHCONVERT_CC__
#define __JOSHCONVERT_CC__
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include "b3d.h"
#define __EVALUATE_SURFACE_VELOCITIES__

using namespace std;
CjoshConvert *CRing::jc=NULL;

void CjoshConvert::Init(string run_name_set){
	CRing::jc=this;
	int ires,iring;
	run_name=run_name_set;
	string parsfilename="model_output/"+run_name+"/fixed_parameters.dat";
	parameter::ReadParsFromFile(parmap,parsfilename);
	reslist=new CResList(&parmap);
	densityf.resize(reslist->resmap.size());
	string inputfilename;
	double pi,ei,sigma2i,dedti,degen,cs2;
	CResInfoMap::iterator rpos;
	CResInfo *resinfo;
	Tf=1000.0*parameter::getD(parmap,"HYDRO_FOTEMP",160.0);
	ETAMAX=parameter::getD(parmap,"B3D_ETAMAX",1.0);
	NRINGSMAX=parameter::getD(parmap,"B3D_NRINGSMAX",100);
	NBOSE=parameter::getI(parmap,"B3D_NBOSE",1);
	boseweight.resize(NBOSE);
	randy=new CRandom(-1234);
	ring.resize(NRINGSMAX+1);
	reslist->CalcEoS(Tf,epsilonf,Pf,nhadronsf,cs2,densityf);
	reslist->CalcEoS(Tf,epsilonf,Pf,nhadronsf,densityf,boseweight);
}

void CjoshConvert::ReadInput(string qualifier_set){
	int alpha,beta,iphi;
	qualifier=qualifier_set;
	string inputfilename="model_output/"+run_name+"/"+qualifier+"/freezeout.dat";
	printf("freezeout info from %s\n",inputfilename.c_str());
	input=fopen(inputfilename.c_str(),"r");
	ReadHeader(input);
	nrings=-1;
	do{
		nrings+=1;
		if(nrings==NRINGSMAX){
			printf("Increase NRINGSMAX\n");
			exit(1);
		}
		ring[nrings].Read(input);
		if(ring[nrings].nread==0){
			ring[nrings].tau=ring[nrings-1].tau;
			for(iphi=0;iphi<ring[nrings].nphi;iphi++){
				ring[nrings].T[iphi]=ring[nrings-1].T[iphi];
				ring[nrings].ux[iphi]=0.0;
				ring[nrings].uy[iphi]=0.0;
				for(alpha=1;alpha<4;alpha++){
					for(beta=alpha;beta<4;beta++)
						ring[nrings].dToverH[iphi][alpha][beta]=ring[nrings-1].dToverH[iphi][alpha][beta];
				}
			}
		}
	}while(ring[nrings].nread!=0);
	fclose(input);
}

void CjoshConvert::WriteNewFormat(){
	int iring;
	string filename="model_output/"+run_name+"/"+qualifier+"/triangles2D.dat";
	triangle_output=fopen(filename.c_str(),"w");
	filename="model_output/"+run_name+"/"+qualifier+"/vertices2D.dat";
	vertex_output=fopen(filename.c_str(),"w");
	double OmegaSum=0;
	
	for(iring=1;iring<=nrings;iring++){
		OmegaSum+=WriteRing(iring);
	}
	fclose(triangle_output);
	fclose(vertex_output);
	//printf("------ OmegaSum=%g -------\n",OmegaSum);
}

double CjoshConvert::WriteRing(int iring){
	double Xmax,Ymax;
	CRing *earlier=&ring[iring-1];
	CRing *later=&ring[iring];
	CResInfo *resinfo;
	CResInfoMap::iterator rpos;
	double dNbarmax,T,Xscale;
	double tau,deltau=0.0,x,delx,y,dely,eta,deleta;
	int ncalls;
	int ires,iquad,itau;
	long long int iv_b1,iv_b2,iv_a1,iv_a2;
	int iphi1,iphi2,alpha,beta,nphi=later->nphi,nparts=0;
	double Omega[3],a[3],b[3];
	double cphi1,sphi1,cphi2,sphi2,phi1,phi2;
	double zsize,dToverH;
	double ux,uy;
	double Omegasum=0.0;
#ifdef __EVALUATE_SURFACE_VELOCITIES__
	double OmegaSquared,vOmega,r;
#endif
	Xmax=Ymax=0.0;
	
	for(iphi1=0;iphi1<nphi;iphi1++){
		iphi2=iphi1+1;
		phi1=0.5*PI*iphi1/double(nphi);
		phi2=0.5*PI*iphi2/double(nphi);
		cphi1=cos(phi1);
		sphi1=sin(phi1);
		cphi2=cos(phi2);
		sphi2=sin(phi2);
		
		iv_b1=(iring-1)*(nphi+1)+iphi1;
		iv_b2=iv_b1+1;
		iv_a1=(iring)*(nphi+1)+iphi1;
		iv_a2=iv_a1+1;
		
		fprintf(vertex_output,"%lld %g %g %g %g %g\n",
		iv_b1,earlier->tau,earlier->r[iphi1]*cphi1,earlier->r[iphi1]*sphi1,earlier->ux[iphi1],earlier->uy[iphi1]);
		if(iphi2==nphi){
			fprintf(vertex_output,"%lld %g %g %g %g %g\n",
			iv_b2,earlier->tau,earlier->r[iphi2]*cphi2,earlier->r[iphi2]*sphi2,earlier->ux[iphi2],earlier->uy[iphi2]);
		}
		a[0]=later->tau-earlier->tau;
		a[1]=earlier->r[iphi1]*cphi1-later->r[iphi1]*cphi1;
		a[2]=earlier->r[iphi1]*sphi1-later->r[iphi1]*sphi1;
		b[0]=a[0];
		b[1]=earlier->r[iphi2]*cphi2-later->r[iphi1]*cphi1;
		b[2]=earlier->r[iphi2]*sphi2-later->r[iphi1]*sphi1;	
		zsize=2.0*ETAMAX*(2.0*earlier->tau+later->tau)/3.0;
		Omega[0]=0.5*(a[1]*b[2]-a[2]*b[1])*zsize;
		Omega[1]=0.5*(a[2]*b[0]-a[0]*b[2])*zsize;
		Omega[2]=0.5*(a[0]*b[1]-a[1]*b[0])*zsize;
#ifdef __EVALUATE_SURFACE_VELOCITIES__
		OmegaSquared=Omega[0]*Omega[0]-Omega[1]*Omega[1]-Omega[2]*Omega[2];
		if(OmegaSquared<0.0){
			vOmega=Omega[0]/sqrt(Omega[1]*Omega[1]+Omega[2]*Omega[2]);
		}
		else
			vOmega=sqrt(Omega[1]*Omega[1]+Omega[2]*Omega[2])/Omega[0];
		x=0.25*(earlier->r[iphi1]*cphi1+earlier->r[iphi2]*cphi2+later->r[iphi1]*cphi1+later->r[iphi2]*cphi2);
		y=0.25*(earlier->r[iphi1]*sphi1+earlier->r[iphi2]*sphi2+later->r[iphi1]*sphi1+later->r[iphi2]*sphi2);
		r=sqrt(x*x+y*y);
		if(fabs(x)>Xmax)
			Xmax=x;
		if(fabs(y)>Ymax)
			Ymax=y;
		
#endif
		Omega[0]*=4; // reflection symmetry
		Omega[1]*=4;
		Omega[2]*=4;
		Omegasum+=Omega[0];
		T=(earlier->T[iphi1]+earlier->T[iphi2]+later->T[iphi1])/3.0;
		Xscale=(earlier->Xscale[iphi1]+earlier->Xscale[iphi2]+later->Xscale[iphi1])/3.0;
		fprintf(triangle_output,"%lld %lld %lld %g %g %g %g %g ",
		iv_b1,iv_b2,iv_a1,Omega[0],Omega[1],Omega[2],T,Xscale);
		for(alpha=1;alpha<4;alpha++){
			for(beta=alpha;beta<4;beta++){
				dToverH=(earlier->dToverH[iphi1][alpha][beta]
					+earlier->dToverH[iphi2][alpha][beta]+later->dToverH[iphi1][alpha][beta])/3.0;
				fprintf(triangle_output,"%g ",dToverH);
			}
		}
		fprintf(triangle_output,"\n");
		
		if(iring!=nrings){
			a[0]=later->tau-earlier->tau;
			a[1]=earlier->r[iphi2]*cphi2-later->r[iphi1]*cphi1;
			a[2]=earlier->r[iphi2]*sphi2-later->r[iphi1]*sphi1;
			b[0]=a[0];
			b[1]=earlier->r[iphi2]*cphi2-later->r[iphi2]*cphi2;
			b[2]=earlier->r[iphi2]*sphi2-later->r[iphi2]*sphi2;
			zsize=2.0*ETAMAX*(earlier->tau+2.0*later->tau)/3.0;
			Omega[0]=-0.5*(a[1]*b[2]-a[2]*b[1])*zsize;
			Omega[1]=-0.5*(a[2]*b[0]-a[0]*b[2])*zsize;
			Omega[2]=-0.5*(a[0]*b[1]-a[1]*b[0])*zsize;
			Omega[0]*=4; // reflection symmetry
			Omega[1]*=4;
			Omega[2]*=4;
			Omegasum+=Omega[0];
			T=(earlier->T[iphi2]+later->T[iphi1]+later->T[iphi2])/3.0;
			Xscale=(earlier->Xscale[iphi2]+later->Xscale[iphi1]+later->Xscale[iphi2])/3.0;
			fprintf(triangle_output,"%lld %lld %lld %g %g %g %g %g ",
			iv_b2,iv_a1,iv_a2,Omega[0],Omega[1],Omega[2],T,Xscale);
			for(alpha=1;alpha<4;alpha++){
				for(beta=alpha;beta<4;beta++){
					dToverH=(earlier->dToverH[iphi2][alpha][beta]
						+later->dToverH[iphi1][alpha][beta]+later->dToverH[iphi2][alpha][beta])/3.0;
					fprintf(triangle_output,"%g ",dToverH);
				}
			}
			fprintf(triangle_output,"\n");
		}
	}
	if(iring==nrings){
		for(iphi1=0;iphi1<=nphi;iphi1++){
			iv_b1=nrings*(nphi+1)+iphi1;
			fprintf(vertex_output,"%lld %g %g %g %g %g\n",iv_b1,later->tau,0.0,0.0,0.0,0.0);
		}
	}
	//printf("iring=%d, tau=%7.4f, DelOmegaSum=%g, (Xmax,Ymax)=(%g,%g)\n",
	//iring,0.5*(earlier->tau+later->tau),Omegasum,Xmax,Ymax);
	return Omegasum;
}

// -----------------------------------------------

void CjoshConvert::ReadHeader(FILE *fptr){
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

CRing::CRing(){
	int iphi,alpha,beta;
	for(iphi=0;iphi<=nphi;iphi++){
		Xscale[iphi]=1.0;
		T[iphi]=jc->Tf;
		for(alpha=0;alpha<4;alpha++){
			for(beta=0;beta<4;beta++){
				dToverH[iphi][alpha][beta]=0.0;
			}
		}
	}
}

void CRing::Read(FILE *fptr){
	int i;
	string dummystring;
	char dummy[120];
	const int NPTS=300;
	double x[NPTS],y[NPTS],uxi[NPTS],uyi[NPTS];
	//double **PiOverh,**PiOverhtilde,Pixxoverh,Pixyoverh,Piyyoverh;
	double dToverH_xx[NPTS],dToverH_xy[NPTS],dToverH_yy[NPTS],XX[NPTS];
	int i1,i2,iphi,alpha,beta;
	double nt,nx,ny,a1,a2,a3,a4,a5,b=0.0,epsilon,P,Rqgp,T;
	const double root3=sqrt(3.0);
	i=-1;	
	fscanf(fptr,"%lf",&tau);
	fscanf(fptr,"%lf",&x[0]);
	do{
		i+=1;
		if(i==300){
			printf("CRing::Read, arrays too small\n");
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
			XX[i]=1000.0*epsilon/jc->epsilonf;
			//printf("XX[%d]=%g\n",i,XX[i]);
			fscanf(fptr,"%s",dummy);
			dummystring.assign(dummy);
			//printf("x,y=(%g,%g)\n",x[i],y[i]);
			if(dummystring!="time:") x[i+1]=atof(dummy);
		}
	}while(dummystring!="time:" && feof(fptr)==false);
	nread=i;
	FillOutArrays(x,y,uxi,uyi,dToverH_xx,dToverH_xy,dToverH_yy,XX);
	//if(nread>0) printf("TIME=%g, %d points for Cooper-Frye ring\n",tau,nread);
}

void CRing::FillOutArrays(double *x,double *y,double *uxi,double *uyi,double *dToverH_xx,double *dToverH_xy,double *dToverH_yy,double *XX){
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
			Xscale[iphi]=w1*XX[i1]+w2*XX[i2];
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

void CRing::Print(){
	int iphi,alpha,beta;
	printf("----- tau=%g ------\n",tau);
	for(iphi=0;iphi<nphi;iphi++){
		printf("%3d: %g %g %g %g %g\n",iphi,T[iphi],Xscale[iphi],r[iphi],ux[iphi],uy[iphi]);
		for(alpha=1;alpha<4;alpha++){
			for(beta=alpha;beta<4;beta++)
				printf("%g ",dToverH[iphi][alpha][beta]);
		}
		printf("\n");
	}
}

#endif
