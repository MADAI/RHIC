#ifndef __RHICSTATMCMC_CC__
#define __RHICSTATMCMC_CC__
#include "rhicstat.h"
using namespace std;

void CRHICStat::Metropolis(long long int nburn,long long int nmcmc,double stepsize){
	NMCMC=nmcmc;
	NBURN=nburn;
	MCMC_STEPSIZE=stepsize;
	Metropolis();
}

void CRHICStat::Metropolis(){
	CRunInfo *bestrun=new CRunInfo(NX,NY);
	unsigned long long int isample,iburn,nsuccess=0,norm=0;
	int NBINS=50,nwrite,ixx,ix,ibin,iz;
	double ll,oldll,bestll=-1.0E99;
	double *x,*oldx,*realx,*oldrealx,*bestx,*xmcmcbar;
	double root12=sqrt(12),*bestz;
	double xpmax=1.0; // 1.0 corresponds to xmin and xmax at boundaries
	bool success;
	double **spread,**xdist;
	string command="mkdir -p figs";
	system(command.c_str());
	if(xmcmc==NULL){
		xmcmc=new double[NX];
		for(ix=0;ix<NX;ix++){
			xmcmc[ix]=xbar[ix];
		}
	}
	x=new double[NX];
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
	FILE *mcmc=fopen("figs/mcmctrace.csv","w");
	for(ix=0;ix<NX;ix++){
		fprintf(mcmc,"\"%s\"",xname[ix].c_str());
		if(ix!=NX-1)
			fprintf(mcmc,",");
	}
	fprintf(mcmc,"\n");
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
		for(ix=0;ix<NX;ix++){
			x[ix]=oldx[ix]+MCMC_STEPSIZE*randy->gauss();
		}
		for(ix=0;ix<NX;ix++){
			realx[ix]=xbar[ix]+x[ix]*(xmax[ix]-xmin[ix])/root12;
			while(realx[ix]>xmax[ix]||realx[ix]<xmin[ix]){
				if(realx[ix]>xmax[ix]) realx[ix]=xmax[ix]-(realx[ix]-xmax[ix]);
				if(realx[ix]<xmin[ix]) realx[ix]=xmin[ix]+(xmin[ix]-realx[ix]);
				x[ix]=(realx[ix]-xbar[ix])*root12/(xmax[ix]-xmin[ix]);
			}
		}
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
		for(ix=0;ix<NX;ix++){
			x[ix]=oldx[ix]+MCMC_STEPSIZE*randy->gauss();
		}
		for(ix=0;ix<NX;ix++){
			realx[ix]=xbar[ix]+x[ix]*(xmax[ix]-xmin[ix])/root12;
			while(realx[ix]>xmax[ix]||realx[ix]<xmin[ix]){
				if(realx[ix]>xmax[ix]) realx[ix]=xmax[ix]-(realx[ix]-xmax[ix]);
				if(realx[ix]<xmin[ix]) realx[ix]=xmin[ix]+(xmin[ix]-realx[ix]);
				x[ix]=(realx[ix]-xbar[ix])*root12/(xmax[ix]-xmin[ix]);
			}
		}
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
		if(nwrite==3){
			for(ix=0;ix<NX;ix++){
				//fprintf(mcmc,"%7.4f, ",(oldrealx[ix]-xmin[ix])/(xmax[ix]-xmin[ix]));
				fprintf(mcmc,"%7.4f, ",oldrealx[ix]);
			}
			fprintf(mcmc,"%g\n",oldll);
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


	printf("----- nsuccess=%lld, success rate=%g, bestll=%g\n",nsuccess,double(nsuccess)/double(NMCMC),bestll);
	bestz=new double[NZ];
	for(ix=0;ix<NX;ix++){
		bestrun->x[ix]=bestx[ix];
		//printf("%25s= %7.4f\n",xname[ix].c_str(),xbar[ix]+(xmax[ix]-xmin[ix])*bestx[ix]/root12);
	}
	zgetter->GetZ(bestx,bestz);
	for(iz=0;iz<NY;iz++){
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
		sprintf(filename,"figs/xdist%d.dat",ix);
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
	for(ix=0;ix<NX;ix++){
		delete [] spread[ix];
		delete [] xdist[ix];
	}
	delete [] spread;
	delete [] xdist;

}

double CRHICStat::GetLL(double *x){
	int ix,iz;
	double dll;
	double *z=new double[NZ];
	zgetter->GetZ(x,z);
	double ll=0.0;
	for(iz=0;iz<NZ;iz++){
		dll=0.5*pow(z[iz]-expinfo->z[iz],2);
		if(USE_COSH_LL){
			dll=log(cosh(2.0*dll));
		}
		ll-=dll;
	}
	delete [] z;
	ll=ll/(1.0+SIGMA2_EMULATOR);
	return ll;
}
	
	

#endif