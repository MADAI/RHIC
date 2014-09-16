#ifndef INCLUDE_SOURCECALC_CC
#define INCLUDE_SOURCECALC_CC

#include "constants.h"
#include "coralutils.h"
#include "sourcecalc.h"

using namespace std;

void CSourceCalc::CalcEffGaussPars(CMCPRList *list,double &Rx,double &Ry,double &Rz,double &Xoff,double &Yoff,double &Zoff){
	int i,norm=0,n=list->GetNMC();
	double *r,*p,t,rr;
	Xoff=Yoff=Zoff=Rx=Ry=Rz=0.0;
	for(i=0;i<n;i++){
		r=list->GetR(i);
		rr=sqrt(r[1]*r[1]+r[2]*r[2]+r[3]*r[3]);
		if(rr<60){
			Xoff+=r[1];
			Yoff+=r[2];
			Zoff+=r[3];
			Rx+=r[1]*r[1];
			Ry+=r[2]*r[2];
			Rz+=r[3]*r[3];
			norm+=1;
		}
	}
	Xoff=Xoff/double(norm);
	Yoff=Yoff/double(norm);
	Zoff=Zoff/double(norm);
	Rx=sqrt((Rx/norm)-Xoff*Xoff);
	Ry=sqrt((Ry/norm)-Yoff*Yoff);
	Rz=sqrt((Rz/norm)-Zoff*Zoff);
	//printf("Rx,y,z=(%g,%g,%g), offset=(%g,%g,%g)\n",Rx,Ry,Rz,Xoff,Yoff,Zoff);
}

void CSourceCalc::CalcEffGaussParsQ2(CMCPRList *list,double &Rx,double &Ry,double &Rz){
	int i,j,alpha,n=list->GetNMC();
	double *ri,*rj,r[4],expfact,rr;
	double numerx=0.0,numery=0.0,numerz=0.0,denom=0.0;
	double afact=1.0;
	double ax2=afact*Rx*Rx,ay2=afact*Ry*Ry,az2=afact*Rz*Rz;
	//double ax2=81.0,ay2=36.0,az2=36.0;
	for(int iterate=0;iterate<4;iterate++){
		for(i=1;i<n;i++){
			ri=list->GetR(i);
			for(j=0;j<i;j++){
				rj=list->GetR(j);
				for(alpha=0;alpha<4;alpha++) r[alpha]=ri[alpha]-rj[alpha];
				rr=r[1]*r[1]+r[2]*r[2]+r[3]*r[3];
				expfact=exp(-0.5*((r[1]*r[1]/ax2)+(r[2]*r[2]/ay2)+(r[3]*r[3]/az2)));
				denom+=expfact;
				numerx+=r[1]*r[1]*expfact;
				numery+=r[2]*r[2]*expfact;
				numerz+=r[3]*r[3]*expfact;
			}
		}
		numerx=numerx/denom;
		numery=numery/denom;
		numerz=numerz/denom;
		Rx=sqrt(0.5*numerx/(1.0-numerx/ax2));
		Ry=sqrt(0.5*numery/(1.0-numery/ay2));
		Rz=sqrt(0.5*numerz/(1.0-numerz/az2));
		printf("Rx,y,z=(%g,%g,%g)\n",Rx,Ry,Rz);
		ax2=afact*Rx*Rx;
		ay2=afact*Ry*Ry;
		az2=afact*Rz*Rz;
		numerx=numery=numerz=denom=0.0;
	}
	printf("---------------------------\n");
}

void CSourceCalc::CalcEffGaussPars(CCHArray *A){
	double Rx,Ry,Rz,Xoff,Yoff,Zoff;
	CalcEffGaussPars(A,Rx,Ry,Rz,Xoff,Yoff,Zoff);
}

void CSourceCalc::CalcEffGaussPars(CCHArray *A,double &Rx,double &Ry,double &Rz,double &Xoff,double &Yoff,double &Zoff){
	double xbar,ybar,zbar,x2bar,y2bar,z2bar,r2bar,r,r2,r3,r4,DELR;
	int ir,NRMAX;
	bool XSYM,YSYM,ZSYM;
	NRMAX=A->GetNRADIAL();
	DELR=A->GetRADSTEP();
	XSYM=A->GetXSYM();
	YSYM=A->GetYSYM();
	ZSYM=A->GetZSYM();
	const double PI=4.0*atan(1.0);
	xbar=ybar=zbar=x2bar=y2bar=z2bar=r2bar=0.0;
	for(ir=0;ir<NRMAX;ir++){
		r=(0.5+ir)*DELR;
		r2=r*r; 
		r3=r2*r;
		r4=r2*r2;
		if(!XSYM) xbar+=r3*A->GetElement(1,0,0,ir);
		if(!YSYM) ybar+=r3*A->GetElement(0,1,0,ir);
		if(!ZSYM) zbar+=r3*A->GetElement(0,0,1,ir);
		x2bar+=r4*A->GetElement(2,0,0,ir);
		y2bar+=r4*A->GetElement(0,2,0,ir);
		z2bar+=r4*A->GetElement(0,0,2,ir);
		r2bar+=r4*A->GetElement(0,0,0,ir);
	}
	xbar*=4.0*PI*DELR/3.0;
	ybar*=4.0*PI*DELR/3.0;
	zbar*=4.0*PI*DELR/3.0;
	double norm=GetNorm(A);
	xbar=xbar/norm;
	ybar=ybar/norm;
	zbar=zbar/norm;
	x2bar=x2bar/norm;
	y2bar=y2bar/norm;
	z2bar=z2bar/norm;
	printf("__________  EFFECTIVE GAUSSIAN PARAMETERS ____________\n");
	printf("Rinv=%g\n",sqrt(2.0*PI*DELR*r2bar/3.0));
	x2bar=4.0*PI*DELR*(2.0*x2bar/15.0+(r2bar/3.0))-xbar*xbar;
	y2bar=4.0*PI*DELR*(2.0*y2bar/15.0+(r2bar/3.0))-ybar*ybar;
	z2bar=4.0*PI*DELR*(2.0*z2bar/15.0+(r2bar/3.0))-zbar*zbar;
	printf("Gaussian distribution with same offsets and 1-part. radii\n");
	printf("offset_xyz=(%g,%g,%g), R_xyz=(%g,%g,%g)\n",
	xbar,ybar,zbar,
	sqrt(fabs(0.5*x2bar)),sqrt(fabs(0.5*y2bar)),sqrt(fabs(0.5*z2bar)));
	printf("______________________________________________________\n");
	
	Xoff=xbar;
	Yoff=ybar;
	Zoff=zbar;
	Rx=sqrt(fabs(0.5*x2bar));
	Ry=sqrt(fabs(0.5*y2bar));
	Rz=sqrt(fabs(0.5*z2bar));
	
}

void CSourceCalc::CalcEffGaussParsPureBose(CMCPRList *list,double &lambda,double &Rx,double &Ry,double &Rz){
	int alpha,beta;
	bool success;
	CGSLMatrix_Real *gslmatrix=new CGSLMatrix_Real(4);
	double **M=new double *[4]; // = dy_i/dx_j, y_i=dchi2/dx_i
	for(alpha=0;alpha<4;alpha++){
		M[alpha]=new double[4]; 
	}
	double *y=new double[4],*x=new double[4],*delx=new double[4],bestx[4];
	double chi2,bestchi2=1.0E50,R2out,R2side,R2long,qmag,progress,qrarg;
	x[0]=lambda; x[1]=Rx*Rx/(HBARC*HBARC); x[2]=Ry*Ry/(HBARC*HBARC); x[3]=Rz*Rz/(HBARC*HBARC);
	
	//* Calculate Model Points */
	int ix,iy,iz,nxyz=25;
	double delq=2.0,qx,qy,qz,*r,q[4],expfact,delxmag,ymag;
	complex<double> ci(0.0,1.0);
	complex<double> ***Ci=new complex<double> **[nxyz];
	double ***C=new double **[nxyz];
	for(ix=0;ix<nxyz;ix++){
		Ci[ix]=new complex<double> *[nxyz];
		C[ix]=new double *[nxyz];
		for(iy=0;iy<nxyz;iy++){
			Ci[ix][iy]=new complex<double>[nxyz];
			C[ix][iy]=new double[nxyz];
			for(iz=0;iz<nxyz;iz++){
				Ci[ix][iy][iz]=0.0;
				C[ix][iy][iz]=0.0;
			}
		}
	}
	int i,n=list->GetNMC();
	printf("nmc=%d\n",n);
	for(ix=0;ix<nxyz;ix++){
		qx=(0.5+ix)*delq;
		for(iy=0;iy<nxyz;iy++){
			qy=(0.5+iy)*delq;
			for(iz=0;iz<nxyz;iz++){
				qz=(0.5+iz)*delq;
				for(i=0;i<n;i++){
					r=list->GetR(i);
					qrarg=2.0*(qx*r[1]+qy*r[2]+qz*r[3])/HBARC;
					if(fabs(qrarg)<200)	Ci[ix][iy][iz]+=exp(ci*qrarg);
				}
				C[ix][iy][iz]=real(Ci[ix][iy][iz]*conj(Ci[ix][iy][iz]));
				C[ix][iy][iz]=(C[ix][iy][iz]-double(n))/double(n*(n-1));
				//printf("%2d %2d %2d: %g\n",ix,iy,iz,C[ix][iy][iz]);
			}
		}
	}
	/*
	for(ix=0;ix<nxyz;ix++){
	qx=(0.5+ix)*delq;
	printf("%6.1f %7.4f %7.4f %7.4f\n",qx,C[ix][0][0],C[0][ix][0],C[0][0][ix]);
	}
	*/
	
	int niterations=0;
	do{
		niterations+=1;
		CalcEffGaussParsPureBose_GetYChi2(C,x,nxyz,y,chi2);
		if(niterations==1) printf("^^^^^ First try, lambda=%g, R=(%g,%g,%g), chi2=%g,  ^^^^^^\n",x[0],HBARC*sqrt(x[1]),HBARC*sqrt(x[2]),HBARC*sqrt(x[3]),chi2);
		CalcEffGaussParsPureBose_GetM(C,x,nxyz,M);
		
		gslmatrix->SolveLinearEqs(y,M,delx);
		delxmag=0.0;
		for(alpha=0;alpha<4;alpha++)
			delxmag+=fabs(delx[alpha]/x[alpha]);
		if(delxmag>0.1){
			for(alpha=0;alpha<4;alpha++)
				delx[alpha]=(0.1/delxmag)*delx[alpha];
		}
		
		progress=0.0;
		for(alpha=0;alpha<4;alpha++){
			progress+=y[alpha]*delx[alpha];
		}
		if(progress<0.0){
			printf("!!!!!!!! negative progress !!!!!!!!!\n");
			for(alpha=0;alpha<4;alpha++) delx[alpha]=-delx[alpha];
		}
		
		for(alpha=0;alpha<4;alpha++) x[alpha]-=delx[alpha];
		
		CalcEffGaussParsPureBose_GetYChi2(C,x,nxyz,y,chi2);
		success=false;
		if(chi2<bestchi2){
			success=true;
			bestchi2=chi2;
			for(alpha=0;alpha<4;alpha++){
				bestx[alpha]=x[alpha];
			}
		}
		
		ymag=0.0;
		for(alpha=0;alpha<4;alpha++){
			ymag+=fabs(y[alpha]);
		}
		
		//printf("delxmag=%g, ymag=%g, chi2=%g, lambda=%g, R=(%g,%g,%g)\n",delxmag,ymag,chi2,x[0],sqrt(x[1])*HBARC,sqrt(x[2])*HBARC,sqrt(x[3])*HBARC);
	}while(niterations<100 && delxmag>1.0E-6);
	
	lambda=bestx[0];
	Rx=HBARC*sqrt(bestx[1]);
	Ry=HBARC*sqrt(bestx[2]);
	Rz=HBARC*sqrt(bestx[3]);
	printf("______ lambda=%g, R=(%g,%g,%g), chi2=%g _______\n",lambda,Rx,Ry,Rz,bestchi2);
		
	//* Clean Up */
	
	delete [] x;
	delete [] y;
	delete [] delx;
	for(alpha=0;alpha<4;alpha++) delete [] M[alpha];
	delete [] M;
	for(ix=0;ix<nxyz;ix++){
		for(iy=0;iy<nxyz;iy++){
			delete [] C[ix][iy];
			delete [] Ci[ix][iy];
		}
		delete [] C[ix];
		delete [] Ci[ix];
	}
	delete [] C;
	delete [] Ci;
	
	if(niterations==100){
		printf("FAILURE in CSourceCalc::CalcEffGaussParsPureBose: niterations=%d\n",niterations);
		exit(1);
	}

	
}

void CSourceCalc::CalcEffGaussParsPureBose_GetM(double ***C,double *x,int nxyz,double **M){
	int alpha,beta,ix,iy,iz;
	double q[4],expfact,qmag,delq=2.0;
	for(alpha=0;alpha<4;alpha++){
		for(beta=0;beta<4;beta++) M[alpha][beta]=0.0;
	}
	for(ix=0;ix<nxyz;ix++){
		q[1]=(0.5+ix)*delq;
		for(iy=0;iy<nxyz;iy++){
			q[2]=(0.5+iy)*delq;
			for(iz=0;iz<nxyz;iz++){
				q[3]=(0.5+iz)*delq;
				qmag=sqrt(q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
				/*
				Only relative momenta Qinv>8 MeV/c are used in fit
				*/
				if(qmag>4.0){
					expfact=0.0;
					for(alpha=1;alpha<4;alpha++)
						expfact+=4.0*q[alpha]*q[alpha]*x[alpha];
					expfact=exp(-expfact);
					M[0][0]+=expfact*expfact;
					for(alpha=1;alpha<4;alpha++){
						M[0][alpha]+=4.0*expfact*q[alpha]*q[alpha]
							*(-2.0*x[0]*expfact+C[ix][iy][iz]);
						for(beta=alpha;beta<4;beta++){
							M[alpha][beta]+=16.0*x[0]*expfact*q[alpha]*q[alpha]*q[beta]*q[beta]
								*(2.0*x[0]*expfact-C[ix][iy][iz]);
						}
					}
				}
			}
		}
	}
	for(alpha=1;alpha<4;alpha++){
		for(beta=0;beta<alpha;beta++) M[alpha][beta]=M[beta][alpha];
	}
}

void CSourceCalc::CalcEffGaussParsPureBose_GetYChi2(double ***C,double *x,int nxyz,double *y,double &chi2){
	int alpha,ix,iy,iz;
	double q[4],expfact,qmag,delq=2.0;
	for(alpha=0;alpha<4;alpha++){
		y[alpha]=0.0;
	}
	chi2=0.0;
	for(ix=0;ix<nxyz;ix++){
		q[1]=(0.5+ix)*delq;
		for(iy=0;iy<nxyz;iy++){
			q[2]=(0.5+iy)*delq;
			for(iz=0;iz<nxyz;iz++){
				q[3]=(0.5+iz)*delq;
				qmag=sqrt(q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
				if(qmag>4.0){
					expfact=0.0;
					for(alpha=1;alpha<4;alpha++)
						expfact+=4.0*q[alpha]*q[alpha]*x[alpha];
					expfact=exp(-expfact);
					chi2+=pow(C[ix][iy][iz]-x[0]*expfact,2.0);
					y[0]+=(x[0]*expfact-C[ix][iy][iz])*expfact;
					for(alpha=1;alpha<4;alpha++){
						y[alpha]-=4.0*x[0]*(x[0]*expfact-C[ix][iy][iz])*q[alpha]*q[alpha]*expfact;
					}
				}
			}
		}
	}
}

#endif
