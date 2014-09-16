#ifndef __LOCALLINEAR_CC__
#define __LOCALLINEAR_CC__
#include "rhicstat.h"
using namespace std;

CZGetter_LocalLinear::CZGetter_LocalLinear(CRHICStat *rsptr): CZGetter(){
	rhicstat=rsptr;
	int NX,NZ,NRUNS,ix,irun;
	NZ=rhicstat->NZ;
	NX=rhicstat->NX;
	NRUNS=rhicstat->NRUNS;
	weight=new double[NRUNS];
	ll_delx=parameter::getD(rhicstat->parmap,"LOCALLINEAR_DELX",0.1);
	ll_slope=new double[NX];
	ll_xz=new double[NX];
	xbar=new double[NX];
	ll_xx=new double *[NX];
	ll_xxinv=new double *[NX];
	for(ix=0;ix<NX;ix++){
		ll_xx[ix]=new double[NX];
		ll_xxinv[ix]=new double[NX];
	}
	//InitLocalLinearFit();
}

void CZGetter_LocalLinear::GetZ(double *x,double *z){
	int ix,jx,iz,irun,NRUNS,NX,NZ;
	CRunInfo **runinfo=rhicstat->runinfo;
	NX=rhicstat->NX;
	NZ=rhicstat->NZ;
	NRUNS=rhicstat->NRUNS;
	double wnorm,zbar;
	CRunInfo *ri;
	for(iz=0;iz<NZ;iz++){
		for(ix=0;ix<NX;ix++){
			ll_xz[ix]=0.0;
			for(jx=0;jx<NX;jx++){
				ll_xx[ix][jx]=ll_xxinv[ix][jx]=0.0;
			}
		}
		wnorm=0.0;
		for(irun=0;irun<NRUNS;irun++){
			ri=runinfo[irun];
			weight[irun]=0.0;
			for(ix=0;ix<NX;ix++){
				weight[irun]+=pow(fabs(x[ix]-ri->x[ix])/ll_delx,2);
			}
			if(weight[irun]<100)
				weight[irun]=1.0/cosh(sqrt(weight[irun]));
			else
				weight[irun]=0.0;
			wnorm+=weight[irun];
		}
		for(irun=0;irun<NRUNS;irun++)
			weight[irun]=weight[irun]/wnorm;

		zbar=0.0;
		for(irun=0;irun<NRUNS;irun++){
			ri=runinfo[irun];
			zbar+=weight[irun]*ri->z[iz];
		}
		for(ix=0;ix<NX;ix++){
			xbar[ix]=0.0;
			for(irun=0;irun<NRUNS;irun++){
				ri=runinfo[irun];
				xbar[ix]+=weight[irun]*ri->x[ix];
			}
		}
		
		for(irun=0;irun<NRUNS;irun++){
			ri=runinfo[irun];
			for(ix=0;ix<NX;ix++){
				ll_xz[ix]+=weight[irun]*(ri->x[ix]-xbar[ix])*(ri->z[iz]-zbar);
				for(jx=0;jx<NX;jx++){
					ll_xx[ix][jx]+=weight[irun]*(ri->x[ix]-xbar[ix])*(ri->x[jx]-xbar[jx]);
				}
			}
		}

		rhicstat->gslmatrix_NX->Invert(ll_xx,ll_xxinv);
		z[iz]=zbar;
		//printf("zbar=%g\n",zbar);
		//printf("  x=(%g,%g,%g,%g,%g,%g)\nxbar=(%g,%g,%g,%g,%g,%g)\n",x[0],x[1],x[2],x[3],x[4],x[5],xbar[0],xbar[1],xbar[2],xbar[3],xbar[4],xbar[5]);
		for(ix=0;ix<NX;ix++){
			for(jx=0;jx<NX;jx++){
				z[iz]+=ll_xz[ix]*ll_xxinv[ix][jx]*(x[jx]-xbar[jx]);
			}
		}
		//printf("zfit[%d]=%g\n",iz,z[iz]);
	}

}

#endif