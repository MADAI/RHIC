#ifndef __RHICSTATMCMC_OPT_CC__
#define __RHICSTATMCMC_OPT_CC__
#include "rhicstat.h"
using namespace std;

void CRHICStat::MetropolisOpt(long long int nburn,long long int nmcmc,double stepsize){
	NMCMC=nmcmc;
	NBURN=nburn;
	MCMC_STEPSIZE=stepsize;
	MetropolisOpt();
}

void CRHICStat::MetropolisOpt(){
	CRunInfo *bestrun=new CRunInfo(NX,NY);
	unsigned long long int isample,iburn,nsuccess=0,norm=0;
	int NBINS=50,nwrite,ixx,ix,ibin,iz;
	double ll,oldll,bestll=-1.0E99,xhatnorm;
	double *x,*oldx,*realx,*oldrealx,*bestx,*xmcmcbar,*rgauss;
	double root12=sqrt(12),*bestz;
	double xpmax=1.0; // 1.0 corresponds to xmin and xmax at boundaries
	bool success,boundcheck;
	double **spread,**xdist;
	if(xmcmc==NULL){
		xmcmc=new double[NX];
		for(ix=0;ix<NX;ix++){
			xmcmc[ix]=xbar[ix];
		}
	}
	if(xhatmcmc==NULL){
		xhatmcmc=new double*[NX];
		for(iz=0;iz<NX;iz++){
			xhatmcmc[iz]=new double[NX];
			for(ix=0;ix<NX;ix++)
				xhatmcmc[iz][ix]=dxdz_inv[iz][ix];
		}
		for(iz=0;iz<NX;iz++){
			xhatnorm=0.0;
			for(ix=0;ix<NX;ix++)
				xhatnorm+=pow(xhatmcmc[iz][ix],2);
			xhatnorm=sqrt(xhatnorm);
			for(ix=0;ix<NX;ix++)
				xhatmcmc[iz][ix]=xhatmcmc[iz][ix]/sqrt(xhatnorm);
		}
	}
	x=new double[NX];
	rgauss=new double [NX+1];
	bestx=new double[NX];
	realx=new double[NX];
	oldx=new double[NX];
	oldrealx=new double[NX];
	xmcmcbar=new double[NX];
	spread=new double *[NX];
	xdist=new double *[NX];
	for(ix=0;ix<NX;ix++){
		spread[ix]=new double[NX]();
		xdist[ix]=new double[NBINS]();
	}
	FILE *mcmc=fopen("mcmctrace.dat","w");
	ll=1.0E-99;
	for(ix=0;ix<NX;ix++){
		oldx[ix]=(xmcmc[ix]-xbar[ix])*root12/(xmax[ix]-xmin[ix]);
		oldrealx[ix]=xbar[ix]+oldx[ix]*(xmax[ix]-xmin[ix])/root12;
	}
	oldll=GetLL(oldx);
	iburn=0;
	while(iburn<NBURN){
		iburn+=1;
		success=false;
		for(iz=0;iz<NX;iz+=2)
			randy->gauss2(&rgauss[iz],&rgauss[iz+1]);
		for(ix=0;ix<NX;ix++){
			x[ix]=oldx[ix];
			for(iz=0;iz<NX;iz++){
				x[ix]+=MCMC_STEPSIZE*rgauss[iz]*xhatmcmc[iz][ix]/sqrt(eigenvalyy[iz]+1.0);
			}
			realx[ix]=xbar[ix]+x[ix]*(xmax[ix]-xmin[ix])/root12;
		}
		boundcheck=true;
		ix=0;
		while(ix<NX && boundcheck==true){
			if(realx[ix]>xmax[ix] || realx[ix]<xmin[ix])
				boundcheck=false;
			ix+=1;
		}
		if(boundcheck){
			ll=GetLL(x);
			if((ll>oldll || randy->ran()<exp(ll-oldll)) && (ll-oldll>-40)){
				success=true;
				nsuccess+=1;
				oldll=ll;
				for(ix=0;ix<NX;ix++){
					oldx[ix]=x[ix];
					oldrealx[ix]=realx[ix];
				}
				if(ll>bestll){
					bestll=ll;
					for(ix=0;ix<NX;ix++) bestx[ix]=x[ix];
						printf("bestll=%g\n",bestll);
				}
			}
		}
	}
	printf("BURN IN FINISHED\n");

	norm=0;
	nwrite=0;
	isample=0;
	while(isample<NMCMC){
		isample+=1;
		if((10*isample%NMCMC)==0){
			printf("----- finished %g percent of MCMC search -----\n",double(100*isample/NMCMC));
			printf("----- nsuccess=%lld, success rate=%g, bestll=%g\n",nsuccess,double(nsuccess)/double(isample),bestll);
		}
		nwrite+=1;
		success=false;
		for(iz=0;iz<NX;iz++)
			randy->gauss2(&rgauss[iz],&rgauss[iz+1]);
		for(ix=0;ix<NX;ix++){
			x[ix]=oldx[ix];
			for(iz=0;iz<NX;iz++){
				x[ix]+=MCMC_STEPSIZE*rgauss[iz]*xhatmcmc[iz][ix]/sqrt(eigenvalyy[iz]+1.0);
				realx[ix]=xbar[ix]+x[ix]*(xmax[ix]-xmin[ix])/root12;
			}
		}
		boundcheck=true;
		ix=0;
		while(ix<NX && boundcheck==true){
			if(realx[ix]>xmax[ix] || realx[ix]<xmin[ix])
				boundcheck=false;
			ix+=1;
		}
		if(boundcheck){
			ll=GetLL(x);
			if((ll>oldll || randy->ran()<exp(ll-oldll)) && (ll-oldll>-40)){
				success=true;
				oldll=ll;
				for(ix=0;ix<NX;ix++){
					oldx[ix]=x[ix];
					oldrealx[ix]=realx[ix];
				}
				if(ll>bestll){
					bestll=ll;
					for(ix=0;ix<NX;ix++) bestx[ix]=x[ix];
						printf("bestll=%g\n",bestll);
				}
			}
		}
		if(nwrite==3){
			for(ix=0;ix<NX;ix++){
				fprintf(mcmc,"%7.4f ",(oldrealx[ix]-xmin[ix])/(xmax[ix]-xmin[ix]));
			}
			fprintf(mcmc,"\n");
			nwrite=0;
			norm+=1;
			for(ix=0;ix<NX;ix++){
				xmcmcbar[ix]+=oldrealx[ix];
				for(ixx=0;ixx<NX;ixx++){
					spread[ix][ixx]+=oldx[ix]*oldx[ixx];
				}
			}
			for(ix=0;ix<NX;ix++){
				ibin=lrint(floor(NBINS*(oldrealx[ix]-xmin[ix])/(xmax[ix]-xmin[ix])));
				if(ibin<0 || ibin>=NBINS){
					printf("ibin=%d is out of range\n",ibin);
					exit(1);
				}
				xdist[ix][ibin]+=1.0;
			}
		}
	}


	//printf("----- nsuccess=%lld, success rate=%g, bestll=%g\n",nsuccess,double(nsuccess)/double(NMCMC),bestll);
	bestz=new double[NZ];
	for(ix=0;ix<NX;ix++){
		bestrun->x[ix]=bestx[ix];
			//printf("%25s= %7.4f\n",xname[ix].c_str(),xbar[ix]+(xmax[ix]-xmin[ix])*bestx[ix]/root12);
	}
	zgetter->GetZ(bestx,bestz);
	for(iz=0;iz<NZ;iz++){
		bestrun->z[iz]=bestz[iz];
	}
	GetYFromZ(bestrun);
	PrintX(bestrun);
	PrintY(bestrun);
	delete [] bestz;

	printf("  xbar=(");
		for(ix=0;ix<NX;ix++){
			xmcmcbar[ix]=xmcmcbar[ix]/double(norm);
			if(ix!=0) printf(",");
			printf("%7.4f",xmcmcbar[ix]);
		}
		printf(")\n");
		printf("SPREAD = \n");
		for(ix=0;ix<NX;ix++){
			for(ixx=0;ixx<NX;ixx++){
				spread[ix][ixx]=spread[ix][ixx]/(12.0*norm);
				spread[ix][ixx]-=xmcmcbar[ix]*xmcmcbar[ixx];
				if(ixx==0) printf(" ");
				printf("%7.4f",12.0*spread[ix][ixx]);
			}
			printf("\n");
		}
		fclose(mcmc);

		FILE *xdfile;
		char filename[80];
		for(ix=0;ix<NX;ix++){
			sprintf(filename,"xdist%d.dat",ix);
			xdfile=fopen(filename,"w");
			fprintf(xdfile,"%s\n",xname[ix].c_str());
			for(ibin=0;ibin<NBINS;ibin++){
				fprintf(xdfile,"%3d %g\n",ibin,NBINS*xdist[ix][ibin]/double(norm));
			}
			fclose(xdfile);
			xmcmc[ix]=oldrealx[ix];
		}

		delete [] oldx;
		delete [] oldrealx;
		delete [] x;
		delete [] bestx;
		delete [] xmcmcbar;
		delete [] rgauss;
		for(ix=0;ix<NX;ix++){
			delete [] spread[ix];
			delete [] xdist[ix];
		}
		delete [] spread;
		delete [] xdist;

	}

	double CRHICStat::GetLL(double *x){
		int ix,iz;
		double *z=new double[NZ];
		zgetter->GetZ(x,z);
		double ll=0.0;
		for(iz=0;iz<NZ;iz++){
			ll-=0.5*pow(z[iz]-expinfo->z[iz],2);
		}
		delete [] z;
		ll=ll/(1.0+SIGMA2_EMULATOR);
		return ll;
	}
	
	/*
	// This is just for testing purposes
	double CRHICStat::GetLL(double *x){
		double likelihood=1.0,da;
		const double root3=sqrt(3.0);
		int ix;
		for(ix=0;ix<NX;ix++){
			da=1.0-fabs(x[ix]/root3);
			likelihood*=da;
		}
		return log(likelihood);
	}
	*/

#endif