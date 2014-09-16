#ifndef __GP_CC__
#define __GP_CC__
#include "coral.h"
#include "rhicstat.h"
using namespace std;

CZGetter_GP::CZGetter_GP(CRHICStat *rsptr) : CZGetter(){
	LINEAROFF=true;
	rhicstat=rsptr;
	hyperR_default=parameter::getD(rhicstat->parmap,"HYPERR_DEFAULT",1.0);
	hyperNugget_default=parameter::getD(rhicstat->parmap,"HYPERNUGGET_DEFAULT",0.2);
	READ_HYPERPARS=parameter::getB(rhicstat->parmap,"READ_HYPERPARS",false);
	LINEAROFF=parameter::getB(rhicstat->parmap,"GP_LINEAROFF",false);
	NORMALIZE=parameter::getB(rhicstat->parmap,"GP_NORMALIZE",true);
	hyperPower=parameter::getD(rhicstat->parmap,"HYPER_POWER",1.8);
	GETERROR=parameter::getB(rhicstat->parmap,"GP_GETERROR",false);
	int irun1,irun2,iz;
	NZ=rhicstat->NZ;
	NX=rhicstat->NX;
	NRUNS=rhicstat->NRUNS;
	NTESTRUNS=rhicstat->NTESTRUNS;
	Cov=new double **[NZ];
	CovInv=new double **[NZ];
	CovInvDotZ=new double *[NZ];
	hyperR=new double *[NZ];
	alphanorm=new double *[NZ];
	hyperTheta0=new double[NZ];
	hyperNugget=new double[NZ];
	Dzsquared=new double[NZ];
	for(iz=0;iz<NZ;iz++){
		CovInvDotZ[iz]=new double[NRUNS];
		alphanorm[iz]=new double[NRUNS];
		hyperR[iz]=new double[NX];
		Cov[iz]=new double*[NRUNS];
		CovInv[iz]=new double *[NRUNS];
		hyperTheta0[iz]=1.0;
		hyperNugget[iz]=0.0;
		for(irun1=0;irun1<NRUNS;irun1++){
			CovInvDotZ[iz][irun1]=0.0;
			alphanorm[iz][irun1]=0.0;
			Cov[iz][irun1]=new double[NRUNS];
			CovInv[iz][irun1]=new double[NRUNS];
			for(irun2=0;irun2<NRUNS;irun2++){
				Cov[iz][irun1][irun2]=0.0;
				CovInv[iz][irun1][irun2]=0.0;
			}
		}
	}
	gslmatrix_NRUNS=new CGSLMatrix_Real(NRUNS);
	InitInterpolator();
	for(iz=0;iz<NZ;iz++){
		for(irun1=0;irun1<NRUNS;irun1++){
			delete [] Cov[iz][irun1];
			delete [] CovInv[iz][irun1];
		}
		delete [] Cov[iz];
		delete [] CovInv[iz];
	}
	delete [] Cov;

	COVINVCALC=false;
}

void CZGetter_GP::LinearFit(){
	int ix,jx,iz,irun;
	int NTESTRUNS=rhicstat->NTESTRUNS;
	CRunInfo **runinfo=rhicstat->runinfo;
	CRunInfo **testinfo=rhicstat->testinfo;
	double *xz=new double[NX],**xx=new double *[NX],**xxinv=new double *[NX];
	double zbar,*xbar=new double[NX];
	for(ix=0;ix<NX;ix++){
		xx[ix]=new double[NX];
		xxinv[ix]=new double[NX];
	}
	if(mlinear==NULL){
		mlinear=new double *[NZ];
		blinear=new double [NZ];
		for(iz=0;iz<NZ;iz++){
			mlinear[iz]=new double[NX];
		}
	}
	for(iz=0;iz<NZ;iz++){
		printf("iz=%d\n",iz);
		zbar=0.0;
		Dzsquared[iz]=0.0;
		for(irun=0;irun<NRUNS;irun++){
			zbar+=runinfo[irun]->z[iz];
		}
		zbar=zbar/double(NRUNS);
		blinear[iz]=zbar;
		for(ix=0;ix<NX;ix++){
			xbar[ix]=0.0;
			for(irun=0;irun<NRUNS;irun++){
				xbar[ix]+=runinfo[irun]->x[ix];
			}
			xbar[ix]=xbar[ix]/double(NRUNS);
		}
		for(ix=0;ix<NX;ix++){
			xz[ix]=0.0;
			for(irun=0;irun<NRUNS;irun++){
				xz[ix]+=(runinfo[irun]->z[iz]-zbar)*(runinfo[irun]->x[ix]-xbar[ix])/double(NRUNS);
			}
			for(jx=0;jx<NX;jx++){
				xx[ix][jx]=0.0;
				for(irun=0;irun<NRUNS;irun++){
					xx[ix][jx]+=(runinfo[irun]->x[ix]-xbar[ix])*(runinfo[irun]->x[jx]-xbar[jx])/double(NRUNS);
				}
			}
		}
		rhicstat->gslmatrix_NX->Invert(xx,xxinv);
		for(ix=0;ix<NX;ix++){
			mlinear[iz][ix]=0.0;
			for(jx=0;jx<NX;jx++){
				mlinear[iz][ix]+=xz[jx]*xxinv[jx][ix];
			}
			blinear[iz]-=mlinear[iz][ix]*xbar[ix];
		}
		for(irun=0;irun<NRUNS;irun++){
			runinfo[irun]->zlinear[iz]=blinear[iz];
			for(ix=0;ix<NX;ix++)
				runinfo[irun]->zlinear[iz]+=mlinear[iz][ix]*runinfo[irun]->x[ix];
		}
		for(irun=0;irun<NTESTRUNS;irun++){
			testinfo[irun]->zlinear[iz]=blinear[iz];
			for(ix=0;ix<NX;ix++)
				testinfo[irun]->zlinear[iz]+=mlinear[iz][ix]*testinfo[irun]->x[ix];
		}
		printf("Linear fit: zbar[%d]=%g, (xbar,xx)=?(0,1) =",iz,zbar);
		for(ix=0;ix<NX;ix++){
			printf(" (%g,%g)",xbar[ix]/double(NRUNS),xx[ix][ix]);
		}
		printf("\n");
		printf("b[%d]=%g",iz,blinear[iz]);
		for(ix=0;ix<NX;ix++){
			printf(", m[%d][%d]=%g",iz,ix,mlinear[iz][ix]);
		}
		printf("\n");
	}
	for(ix=0;ix<NX;ix++){
		delete [] xx[ix];
		delete [] xxinv[ix];
	}
	delete [] xx;
	delete [] xxinv;
	delete [] xbar;
	delete [] xz;
}

void CZGetter_GP::CalcCovInv(){
	int iz,irun1,irun2,irun3;
	/*
	double **unitytest;
	unitytest=new double *[NRUNS];
	for(irun1=0;irun1<NRUNS;irun1++){
		unitytest[irun1]=new double[NRUNS];
	}
	*/
	printf("BEGIN COV[][] INVERSION\n");
		//gslmatrix_NRUNS->Cholesky_Invert(Cov[iz],CovInv[iz]);
	for(iz=0;iz<NZ;iz++)
		gslmatrix_NRUNS->Invert(Cov[iz],CovInv[iz]);
	/*
	for(irun1=0;irun1<NRUNS;irun1++){
		for(irun2=0;irun2<NRUNS;irun2++){
			unitytest[irun1][irun2]=0.0;
			for(irun3=0;irun3<NRUNS;irun3++){
				unitytest[irun1][irun2]+=Cov[0][irun1][irun3]*CovInv[iz][irun3][irun2];
			}
			printf("%7.5f ",unitytest[irun1][irun2]);
		}
		printf("\n");
		delete [] unitytest[irun1];
	}
	delete [] unitytest;
	*/
	printf("INVERSION FINISHED\n");
}

void CZGetter_GP::InitInterpolator(){
	int iz,ix,id,irun1,irun2;
	double dz;
	CRunInfo **runinfo=rhicstat->runinfo;

	if(mlinear==NULL) LinearFit();
	double *zz=new double[NRUNS],netcov;
	FILE *hyperfile;
	if(READ_HYPERPARS){
		ReadHyperPars();
	}
	else{
		for(iz=0;iz<NZ;iz++){
			for(ix=0;ix<NX;ix++){
				hyperR[iz][ix]=hyperR_default;
			}
			hyperNugget[iz]=hyperNugget_default;
		}
		for(iz=0;iz<NZ;iz++){
			for(irun1=0;irun1<NRUNS;irun1++){
				if(LINEAROFF){
					dz=runinfo[irun1]->z[iz];
				}
				else{
					dz=runinfo[irun1]->z[iz]-runinfo[irun1]->zlinear[iz];
				}
				Dzsquared[iz]+=dz*dz;
			}
			hyperTheta0[iz]=Dzsquared[iz]/double(NRUNS);
		}
	}
	for(iz=0;iz<NZ;iz++){
		for(irun1=0;irun1<NRUNS;irun1++){
			netcov=0.0;
			for(irun2=0;irun2<NRUNS;irun2++){
				Cov[iz][irun1][irun2]=GetCov(iz,runinfo[irun1],runinfo[irun2]);
				if(irun1!=irun2 && Cov[iz][irun1][irun2]>0.2){
					netcov+=Cov[iz][irun1][irun2];
					//printf("Cov[%d][%d][%d]=%g\n",iz,irun1,irun2,Cov[iz][irun1][irun2]);
				}
			}
			//if(netcov<1.0)
				//printf("netcov[%d][%d]=%g\n",iz,irun1,netcov);
		}
		
		for(irun1=0;irun1<NRUNS;irun1++)
			if(LINEAROFF) zz[irun1]=runinfo[irun1]->z[iz];
		else zz[irun1]=runinfo[irun1]->z[iz]-runinfo[irun1]->zlinear[iz];
		gslmatrix_NRUNS->SolveLinearEqs(zz,Cov[iz],CovInvDotZ[iz]);
		for(irun1=0;irun1<NRUNS;irun1++)
			zz[irun1]=1.0;
		gslmatrix_NRUNS->SolveLinearEqs(zz,Cov[iz],alphanorm[iz]);

			/*
		for(irun1=0;irun1<NRUNS;irun1++){
			CovInvDotZ[iz][irun1]=0.0;
			for(irun2=0;irun2<NRUNS;irun2++){
				CovInvDotZ[iz][irun1]+=CovInv[iz][irun1][irun2]*(runinfo[irun2]->z[iz]-runinfo[irun2]->zlinear[iz]);
			}
		}
		*/

	}
	delete [] zz;
}

double CZGetter_GP::GetCov(int iz,CRunInfo *runinfo1,CRunInfo *runinfo2){
	return GetCov(iz,runinfo1->x,runinfo2->x);
}

double CZGetter_GP::GetCov(int iz,double *x1,double *x2){
	double sum=0.0,dx,answer;
	int ix;
	for(ix=0;ix<rhicstat->NX;ix++){
		dx=x1[ix]-x2[ix];
		//sum+=fabs(dx);
		sum+=pow(fabs(dx)/hyperR[iz][ix],hyperPower);
	}
	if(sum>200) sum=200;
	answer=hyperTheta0[iz]*exp(-0.5*sum);
	if(x1==x2) answer+=hyperNugget[iz];
	//answer=1.0-sum/(hyperR_default);
	return answer;
}

void CZGetter_GP::GetCovDerivs(int iz,double *x1,double *x2,double &C,double *DC,double **DDC){
	double sum=0.0,dx1,dx2;
	int ix1,ix2,NX=rhicstat->NX;
	for(ix1=0;ix1<NX;ix1++){
		DC[ix1]=0.0;
		for(ix2=0;ix2<NX;ix2++){
			DDC[ix1][ix2]=0.0;
		}
	}
	for(ix1=0;ix1<rhicstat->NX;ix1++){
		dx1=x1[ix1]-x2[ix1];
		sum+=pow(fabs(dx1)/hyperR[iz][ix1],hyperPower);
	}
	if(sum>200) sum=200;
	C=exp(-0.5*sum);
	for(ix1=0;ix1<NX;ix1++){
		dx1=fabs(x1[ix1]-x2[ix1]);
		DC[ix1]=C*0.5*pow(dx1/hyperR[iz][ix1],hyperPower)*hyperPower/hyperR[iz][ix1];
		for(ix2=0;ix2<=ix1;ix2++){
			dx2=fabs(x2[ix1]-x2[ix2]);
			DDC[ix1][ix2]=DC[ix1]*0.5*pow(dx1/hyperR[iz][ix2],hyperPower)*hyperPower/hyperR[iz][ix2];
			if(ix2==ix1){
				DDC[ix1][ix2]+=DC[ix1]*(-hyperPower-1.0)/hyperR[iz][ix1];
			}
		}
	}

}

void CZGetter_GP::GetZ(double *x,double *z){
	int iz,izz,irun,ix;
	double norm,cov;
	for(iz=0;iz<NZ;iz++){
		z[iz]=blinear[iz];
		for(ix=0;ix<NX;ix++){
			z[iz]+=mlinear[iz][ix]*x[ix];
		}
		if(LINEAROFF) z[iz]=0.0;
		else{
			z[iz]=blinear[iz];
			for(ix=0;ix<NX;ix++)
				z[iz]+=mlinear[iz][ix]*x[ix];
		}
		//printf("In GetZ, zlinear=%g\n",z[iz]);
		norm=0.0;
		for(irun=0;irun<NRUNS;irun++){
			cov=GetCov(iz,x,rhicstat->runinfo[irun]->x);
			z[iz]+=cov*CovInvDotZ[iz][irun];
			if(NORMALIZE)
				norm+=cov*alphanorm[iz][irun];
		}
		if(NORMALIZE)
			z[iz]=z[iz]/norm;
	}
}

void CZGetter_GP::GetZError(double *x,double *error){
	if(COVINVCALC==false){
		CalcCovInv();
		COVINVCALC=true;
	}
	int iz,irun1,irun2;
	double *x1,*x2;
	for(iz=0;iz<NZ;iz++){
		error[iz]=GetCov(iz,x,x);
		for(irun1=0;irun1<NRUNS;irun1++){
			x1=rhicstat->runinfo[irun1]->x;
			for(irun2=0;irun2<NRUNS;irun2++){
				x2=rhicstat->runinfo[irun2]->x;
				error[iz]-=GetCov(iz,x,x1)*CovInv[iz][irun1][irun2]*GetCov(iz,x2,x);
			}
		}
	}
}

void CZGetter_GP::ReadHyperPars(){
	const double root3=sqrt(3.0);
	char dummy[150];
	int iz,ix,id;
	double scale=1.0,nugget;
	string filename=parameter::getS(rhicstat->parmap,"HYPERPARS_FILENAME","statinfo/hyperpars.dat");
	FILE *hyperfile=fopen(filename.c_str(),"r");
	fgets(dummy,150,hyperfile);		fgets(dummy,150,hyperfile);
	for(iz=0;iz<NZ;iz++){
		fscanf(hyperfile,"%d %lf %lf %lf",&id,&hyperTheta0[iz],&scale,&hyperNugget[iz]);
		//hyperTheta0[iz]*=30;
		for(ix=0;ix<NX;ix++){
			fscanf(hyperfile,"%lf",&hyperR[iz][ix]);
			hyperR[iz][ix]=hyperR[iz][ix]*(2*root3);
			printf("hyperR[%d][%d]=%g\n",iz,ix,hyperR[iz][ix]);
		}
	}
	fclose(hyperfile);
}

void CZGetter_GP::PrintHyperPars(){
	const double root3=sqrt(3.0);
	char dummy[150];
	int iz,ix,id;
	double scale,nugget;
	string filename=parameter::getS(rhicstat->parmap,"NEWHYPERPARS_FILENAME","statinfo/newhyperpars.dat");
	FILE *hyperfile=fopen(filename.c_str(),"w");
	//printf("##-- EMULATOR LENGTH SCALES (thetas) IN PCA SPACE -- ##\n");

	//printf("##id             pca-var         Scale             Nugget                length_0        length_1      length_2        length_3        length_4        length_5 -- ##\n");
	for(iz=0;iz<NZ;iz++){
		fprintf(hyperfile,"%d %lf %lf %lf",id,hyperTheta0[iz],scale,hyperNugget[iz]);
		for(ix=0;ix<NX;ix++){
			fprintf(hyperfile,"%lf",hyperR[iz][ix]);
			hyperR[iz][ix]=hyperR[iz][ix]/(2*root3);
		}
		fprintf(hyperfile,"\n");
	}
	fclose(hyperfile);
}

// Dead below
//JFN 11/30/12- Well then. lets just comment it out
/*void CZGetter_GP::CalcHyperPars(){
	//Use Newton's method
	int iz,ix1,ix2,irun1,irun2;
	// solve eqs for delxx, yy=MM*xx+b, delyy=MM*delxx=-yy, 
	double *yy,**MM,*xx,*delxx;
	double C,***DC,****DDC;
	DC=new double **[NRUNS];
	DDC=new double ***[NRUNS];
	for(irun1=0;irun1<NRUNS;irun1++){
		DC[irun1]=new double *[NRUNS];
		DDC[irun1]=new double **[NRUNS];
		for(irun2=0;irun2<NRUNS;irun2++){
			DC[irun1][irun2]=new double[NX];
			DDC[irun1][irun2]=new double *[NX];
			for(ix1=0;ix1<NX;ix1++){
				DDC[irun1][irun2][ix1]=new double[NX];
			}
		}
	}
	for(iz=0;iz<NZ;iz++){
		// Initial guess
		for(ix1=0;ix1<NX;ix1++){
			xx[ix1]=0.3;
		}
	}

	gslmatrix_NRUNS->Cholesky_Invert(Cov[iz],CovInv[iz]);


	for(irun1=0;irun1<NRUNS;irun1++){
		for(irun2=0;irun2<NRUNS;irun2++){
			for(ix1=0;ix1<NX;ix1++){
				delete [] DDC[irun1][irun2][ix1];
			}
			delete [] DC[irun1][irun2];
			delete [] DDC[irun1][irun2];
		}
		delete [] DC[irun1];
		delete [] DDC[irun1];
	}
	delete [] DC;
	delete [] DDC;


}*/

#endif