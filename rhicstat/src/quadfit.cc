#ifndef __QUADFIT_CC__
#define __QUADFIT_CC__
#include "coral.h"
#include "rhicstat.h"
using namespace std;

CZGetter_QuadFit::CZGetter_QuadFit(CRHICStat *rsptr) : CZGetter(){
	rhicstat=rsptr;
	InitQuadFit();
}

void CZGetter_QuadFit::InitQuadFit(){
	int NX=rhicstat->NX,NZ=rhicstat->NZ,NRUNS=rhicstat->NRUNS;
	CRunInfo **runinfo=rhicstat->runinfo;
	int i,j,i1,i2,j1,j2,iz,irun,n=(NX+1)*(NX+2)/2;
	double *w=new double[NRUNS];
	double besterror=1.0E99;
	int ibestrun=0;
	for(irun=0;irun<NRUNS;irun++){
		if(runinfo[irun]->good)	rhicstat->CalcNetDiffExp(runinfo[irun]);
		w[irun]=exp(-0.2*runinfo[irun]->netdiff_exp);
		if(runinfo[irun]->netdiff_exp<besterror){
			ibestrun=irun;
			besterror=runinfo[irun]->netdiff_exp;
		}
	}
	printf("best model run = run%d, besterror^2=%g\n",ibestrun+1,besterror);
	rhicstat->PrintX(runinfo[ibestrun]);
	rhicstat->PrintY(runinfo[ibestrun]);
	//Misc::Pause();
	CRunInfo *ri;
	double **M=new double *[n];
	double *error=new double[NRUNS];
	for(i=0;i<n;i++){
		M[i]=new double[n];
	}
	double *yy=new double[n];
	double *xx=new double[n];
	Aquad=new double **[NX];
	for(iz=0;iz<NX;iz++){
		Aquad[iz]=new double*[NX];
		for(i1=0;i1<NX;i1++){
			Aquad[iz][i1]=new double[NX];
			for(i2=0;i2<NX;i2++) Aquad[iz][i1][i2]=0.0;
		}
}
Bquad=new double*[NX];
for(iz=0;iz<NX;iz++){
	Bquad[iz]=new double[NX];
	for(i1=0;i1<NX;i1++) Bquad[iz][i1]=0.0;
}
Cquad=new double[NX];
for(iz=0;iz<NX;iz++) Cquad[iz]=0.0;
	CGSLMatrix_Real *gslmatrix=new CGSLMatrix_Real(n);

	/** for checking
	for(irun=0;irun<NRUNS;irun++){
		ri=runinfo[irun];
		for(iz=0;iz<NX;iz++){
			ri->z[iz]=-1.0+5*iz;
			for(i1=0;i1<NX;i1++){
				for(i2=0;i2<=i1;i2++) ri->z[iz]+=0.5*(iz+1+i1-0.12345*i2)*ri->x[i1]*ri->x[i2];
				ri->z[iz]+=0.6*(2.3352-iz)*(i1-0.25252*i2)*ri->x[i1];
			}
		}
	}
	*/
	
	for(iz=0;iz<NX;iz++){

		for(i=0;i<n;i++){
			yy[i]=xx[i]=0.0;
		}
		for(i=0;i<n;i++){
			for(j=0;j<n;j++) M[i][j]=0.0;
		}
			// from dM/dA_ij=0
	i=0;
	for(i1=0;i1<NX;i1++){
		for(i2=0;i2<=i1;i2++){
			j=0;
			for(j1=0;j1<NX;j1++){
				for(j2=0;j2<=j1;j2++){
					for(irun=0;irun<NRUNS;irun++){
						if(runinfo[irun]->good) M[i][j]+=w[irun]*runinfo[irun]->x[i1]*runinfo[irun]->x[i2]*runinfo[irun]->x[j1]*runinfo[irun]->x[j2];
					}
					j+=1;
				}
			}
			j=NX*(NX+1)/2;
			for(j1=0;j1<NX;j1++){
				for(irun=0;irun<NRUNS;irun++){
					if(runinfo[irun]->good) M[i][j]+=w[irun]*runinfo[irun]->x[i1]*runinfo[irun]->x[i2]*runinfo[irun]->x[j1];
				}
				j+=1;
			}
			j=-1+(NX+1)*(NX+2)/2;
			for(irun=0;irun<NRUNS;irun++){
				if(runinfo[irun]->good) M[i][j]+=w[irun]*runinfo[irun]->x[i1]*runinfo[irun]->x[i2];
			}
			yy[i]=0.0;
			for(irun=0;irun<NRUNS;irun++){
				if(runinfo[irun]->good) yy[i]+=w[irun]*runinfo[irun]->x[i1]*runinfo[irun]->x[i2]*runinfo[irun]->z[iz];
			}
			i+=1;
		}
	}
			// from dM/dB_i=0
	i=NX*(NX+1)/2;
	for(i1=0;i1<NX;i1++){
		j=0;
		for(j1=0;j1<NX;j1++){
			for(j2=0;j2<=j1;j2++){
				for(irun=0;irun<NRUNS;irun++){
					if(runinfo[irun]->good) M[i][j]+=w[irun]*runinfo[irun]->x[i1]*runinfo[irun]->x[j1]*runinfo[irun]->x[j2];
				}
				j+=1;
			}
		}
		j=NX*(NX+1)/2;
		for(j1=0;j1<NX;j1++){
			for(irun=0;irun<NRUNS;irun++){
				if(runinfo[irun]->good) M[i][j]+=w[irun]*runinfo[irun]->x[i1]*runinfo[irun]->x[j1];
			}
			j+=1;
		}
		j=-1+(NX+1)*(NX+2)/2;
		for(irun=0;irun<NRUNS;irun++){
			if(runinfo[irun]->good) M[i][j]+=w[irun]*runinfo[irun]->x[i1];
		}
		yy[i]=0.0;
		for(irun=0;irun<NRUNS;irun++){
			if(runinfo[irun]->good) yy[i]+=w[irun]*runinfo[irun]->x[i1]*runinfo[irun]->z[iz];
		}
		i+=1;
	}
			// from dM/dC=0
	i=-1+(NX+1)*(NX+2)/2;
	j=0;
	for(j1=0;j1<NX;j1++){
		for(j2=0;j2<=j1;j2++){
			for(irun=0;irun<NRUNS;irun++){
				if(runinfo[irun]->good) M[i][j]+=w[irun]*runinfo[irun]->x[j1]*runinfo[irun]->x[j2];
			}
			j+=1;
		}
	}
	j=NX*(NX+1)/2;
	for(j1=0;j1<NX;j1++){
		for(irun=0;irun<NRUNS;irun++){
			if(runinfo[irun]->good) M[i][j]+=w[irun]*runinfo[irun]->x[j1];
		}
		j+=1;
	}
	j=-1+(NX+1)*(NX+2)/2;
	for(irun=0;irun<NRUNS;irun++){
		if(runinfo[irun]->good) M[i][j]+=w[irun];
	}
	yy[i]=0.0;
	for(irun=0;irun<NRUNS;irun++){
		if(runinfo[irun]->good) yy[i]+=w[irun]*runinfo[irun]->z[iz];
	}

	gslmatrix->SolveLinearEqs(yy,M,xx);

	i=0;
	for(i1=0;i1<NX;i1++){
		for(i2=0;i2<=i1;i2++){
			Aquad[iz][i1][i2]=xx[i];
					//Aquad[iz][i2][i1]=xx[i];
			i+=1;
		}
	}
	i=(NX+1)*NX/2;
	for(i1=0;i1<NX;i1++){
		Bquad[iz][i1]=xx[i];
		i+=1;
	}
	i=-1+(NX+1)*(NX+2)/2;
	Cquad[iz]=xx[i];
}
	/**
	for(iz=0;iz<NX;iz++){
		printf("---------------- iz=%d ---------------------\n",iz);
		printf("Aquad\n");
		for(i1=0;i1<NX;i1++){
			for(i2=0;i2<NX;i2++){
				printf("%8.4f ",Aquad[iz][i1][i2]);
			}
			printf("\n");
		}
		printf("Bquad\n");
		for(i2=0;i2<NX;i2++){
			printf("%8.4f ",Bquad[iz][i2]);
		}
		for(i2=0;i2<NX;i2++){
			printf("%8.4f ",Bquad[iz][i2]);
		}
		printf("\n");
		printf("Cquad=%8.4f\n",Cquad[iz]);
	}
	 */
	
	
	for(irun=0;irun<NRUNS;irun++){
		ri=runinfo[irun];
		if(ri->good){
			error[irun]=0.0;
			GetZ(ri->x,ri->zfit);
			rhicstat->CalcNetDiffFit(ri);
			rhicstat->CalcNetDiffFitExp(ri);

			printf("irun=%d: netdiff_fit=%g, netdiff_exp=%g, netdiff_fitexp=%g\n",irun,ri->netdiff_fit,ri->netdiff_exp,ri->netdiff_fitexp);
			
			if(error[irun]>10.0){
				printf("----- irun=%3d: error^2=%g -----\n",irun,error[irun]);
				rhicstat->PrintX(ri);
				rhicstat->PrintY(ri);
				printf("zfit = %8.4f=?%8.4f %8.4f=?%8.4f %8.4f=?%8.4f, %8.4f=?%8.4f, %8.4f=?%8.4f, %8.4f=?%8.4f\n",ri->zfit[0],ri->z[0],ri->zfit[1],ri->z[1],ri->zfit[2],ri->z[2],ri->zfit[3],ri->z[3],ri->zfit[4],ri->z[4],ri->zfit[5],ri->z[5]);
			}

		}
	}
	
	delete gslmatrix;
	delete [] error;
	delete [] w;
}

void CZGetter_QuadFit::GetZ(double *x,double *z){
	int iz,i1,i2;
	for(iz=0;iz<rhicstat->NX;iz++){
		z[iz]=Cquad[iz];
		for(i1=0;i1<rhicstat->NX;i1++){
			z[iz]+=Bquad[iz][i1]*x[i1];
			for(i2=0;i2<=i1;i2++)	z[iz]+=Aquad[iz][i1][i2]*x[i1]*x[i2];			
		}
}
}

void CZGetter_QuadFit::PrintQuadCode(){
	int iz,i1,i2,iy;
	int NZ=rhicstat->NZ,NX=rhicstat->NX, NY=rhicstat->NY;
	CRunInfo **runinfo=rhicstat->runinfo;
	double ***A=new double **[NY];
	double **B=new double *[NY];
	double *C=new double[NY]();
	for(iy=0;iy<NY;iy++){
		A[iy]=new double *[NX];
		B[iy]=new double [NX]();
		for(i1=0;i1<NX;i1++) A[iy][i1]=new double[NX]();
	}
for(iy=0;iy<NY;iy++){
	for(i1=0;i1<NX;i1++){
		for(i2=0;i2<NX;i2++){
			A[iy][i1][i2]=0.0;
			for(iz=0;iz<NZ;iz++) A[iy][i1][i2]+=rhicstat->sigmaybar[iy]*rhicstat->Uytoz_inv[iz][iy]*Aquad[iz][i1][i2];
		}
	B[iy][i1]=0.0;
	for(iz=0;iz<NZ;iz++) B[iy][i1]+=rhicstat->sigmaybar[iy]*rhicstat->Uytoz_inv[iz][iy]*Bquad[iz][i1];
}
C[iy]=rhicstat->ybar[iy];
for(iz=0;iz<NZ;iz++) C[iy]+=rhicstat->sigmaybar[iy]*rhicstat->Uytoz_inv[iz][iy]*Cquad[iz];
}

double *ytest=new double[NY];
int itest=28;
for(iy=0;iy<NY;iy++){
	ytest[iy]=C[iy];
	for(i1=0;i1<NX;i1++){
		ytest[iy]+=B[iy][i1]*runinfo[itest]->x[i1];
		for(i2=0;i2<NX;i2++){
			ytest[iy]+=A[iy][i1][i2]*runinfo[itest]->x[i1]*runinfo[itest]->x[i2];
		}
	}
	printf("ytest[%d]=%g =? %g, sigmay=%g, ybar=%g\n",iy,ytest[iy],runinfo[itest]->y[iy]*rhicstat->sigmaybar[iy]+rhicstat->ybar[iy],rhicstat->sigmaybar[iy],rhicstat->ybar[iy]);
}
delete [] ytest;

FILE *qfile=fopen("quad.cc","w");
fprintf(qfile,"#include <cstdlib>\n");
fprintf(qfile,"#include <cstdio>\n");
fprintf(qfile,"#include <cmath>\n");
fprintf(qfile,"using namespace std;\n");

fprintf(qfile,"int main(int argc, char *argv[]){\n");
fprintf(qfile,"\tconst int NX=%d,NY=%d;\n",NX,NY);
fprintf(qfile,"\tint iy,i1,i2;\n");
fprintf(qfile,"\tdouble A[NY][NX][NX],B[NY][NX],C[NY];\n");
for(iy=0;iy<NY;iy++){
	fprintf(qfile,"\tC[%d]=%g;\n",iy,C[iy]);
	for(i1=0;i1<NX;i1++){
		fprintf(qfile,"\tB[%d][%d]=%g;\n",iy,i1,B[iy][i1]);
		for(i2=0;i2<NX;i2++)	fprintf(qfile,"\tA[%d][%d][%d]=%g;\n",iy,i1,i2,A[iy][i1][i2]);
	}
}

fprintf(qfile,"\tdouble xmin[NX],xmax[NX],xbar[NX],x[NX],xprime[NX];\n");
for(i1=0;i1<NX;i1++){
	fprintf(qfile,"\txmin[%d]=%g; xmax[%d]=%g; xbar[%d]=%g;\n",i1,rhicstat->xmin[i1],i1,rhicstat->xmax[i1],i1,rhicstat->xbar[i1]);
	fprintf(qfile,"\tx[%d]=atof(argv[%d]);\n",i1,i1+1);
}
fprintf(qfile,"\tfor(i1=0;i1<NX;i1++){\n");
fprintf(qfile,"\t\txprime[i1]=(x[i1]-xbar[i1])*sqrt(12.0)/(xmax[i1]-xmin[i1]);\n");
fprintf(qfile,"\t}\n");

fprintf(qfile,"\tdouble yexp[NY]={");
for(iy=0;iy<NY;iy++){
	if(iy!=0)
		fprintf(qfile,",");
	fprintf(qfile,"%g",rhicstat->ybar[iy]+rhicstat->sigmaybar[iy]*rhicstat->expinfo->y[iy]);
}
fprintf(qfile,"};\n");

fprintf(qfile,"\tdouble y[NY];\n");
fprintf(qfile,"\tfor(iy=0;iy<NY;iy++){\n");
fprintf(qfile,"\t\ty[iy]=C[iy];\n");
fprintf(qfile,"\t\tfor(i1=0;i1<NX;i1++){\n");
fprintf(qfile,"\t\t\ty[iy]+=B[iy][i1]*xprime[i1];\n");
fprintf(qfile,"\t\t\tfor(i2=0;i2<=i1;i2++) y[iy]+=A[iy][i1][i2]*xprime[i1]*xprime[i2];\n");
fprintf(qfile,"\t\t}\n\t};\n");

fprintf(qfile,"\tprintf(\"");
	for(iy=0;iy<NY;iy++){
		if(iy!=0) fprintf(qfile," ");
		fprintf(qfile,"%%g");
	}
	fprintf(qfile,"\\n\",");
	for(iy=0;iy<NY;iy++){
		if(iy!=0) fprintf(qfile,",");
		fprintf(qfile,"y[%d]",iy);
	}
	fprintf(qfile,");\n");
	
	fprintf(qfile,"\tprintf(\"");
		for(iy=0;iy<NY;iy++){
			if(iy!=0) fprintf(qfile," ");
			fprintf(qfile,"%g",rhicstat->sigmaybar[iy]);
		}
		fprintf(qfile,"\\n\");\n");

		fprintf(qfile,"\treturn 0;\n}\n");
		fclose(qfile);
		for(iy=0;iy<NY;iy++){
			for(i1=0;i1<NX;i1++) delete [] A[iy][i1];
				delete [] B[iy];
		}
		delete [] A;
		delete [] B;
		delete [] C;
	}

	void CZGetter_QuadFit::PrintCoefficients(){
		CRunInfo **runinfo=rhicstat->runinfo;
		int NX=rhicstat->NX,NY=rhicstat->NY;
		int iz,i1,i2,iy;
		const double root12=sqrt(12.0);
		int NZ=NX;
		double ***A=new double **[NY];
		double **B=new double *[NY];
		double *C=new double[NY]();
		for(iy=0;iy<NY;iy++){
			A[iy]=new double *[NX];
			B[iy]=new double [NX]();
			for(i1=0;i1<NX;i1++) A[iy][i1]=new double[NX]();
		}
	for(iy=0;iy<NY;iy++){
		for(i1=0;i1<NX;i1++){
			for(i2=0;i2<NX;i2++){
				A[iy][i1][i2]=0.0;
				for(iz=0;iz<NZ;iz++) A[iy][i1][i2]+=12.0*rhicstat->sigmaybar[iy]*rhicstat->Uytoz_inv[iz][iy]*Aquad[iz][i1][i2];
			}
		B[iy][i1]=0.0;
		for(iz=0;iz<NZ;iz++) B[iy][i1]+=root12*rhicstat->sigmaybar[iy]*rhicstat->Uytoz_inv[iz][iy]*Bquad[iz][i1];
	}
C[iy]=rhicstat->ybar[iy];
for(iz=0;iz<NZ;iz++) C[iy]+=rhicstat->sigmaybar[iy]*rhicstat->Uytoz_inv[iz][iy]*Cquad[iz];
}

double *ytest=new double[NY];
int itest=28;
for(iy=0;iy<NY;iy++){
	ytest[iy]=C[iy];
	for(i1=0;i1<NX;i1++){
		ytest[iy]+=B[iy][i1]*runinfo[itest]->x[i1]/root12;
		for(i2=0;i2<NX;i2++){
			ytest[iy]+=A[iy][i1][i2]*runinfo[itest]->x[i1]*runinfo[itest]->x[i2]/12.0;
		}
	}
	printf("ytest[%d]=%g =? %g, sigmay=%g, ybar=%g\n",iy,ytest[iy],runinfo[itest]->y[iy]*rhicstat->sigmaybar[iy]+rhicstat->ybar[iy],rhicstat->sigmaybar[iy],rhicstat->ybar[iy]);
}
delete [] ytest;

FILE *fptr=fopen("jeff.dat","w");

fprintf(fptr,"%d %d\n",NX,NY);
for(iy=0;iy<NY;iy++){
	for(i1=0;i1<NX;i1++){
		for(i2=i1;i2<NX;i2++){
			fprintf(fptr,"%g ",A[iy][i2][i1]);
		}
	}
	fprintf(fptr,"\n");
	for(i1=0;i1<NX;i1++) fprintf(fptr,"%g ",B[iy][i1]);
		fprintf(fptr,"\n%g\n",C[iy]);
	fprintf(fptr,"%g %g\n",rhicstat->expinfo->y[iy],rhicstat->sigmaybar[iy]);
}
fclose(fptr);
for(iy=0;iy<NY;iy++){
	for(i1=0;i1<NX;i1++) delete [] A[iy][i1];
		delete [] B[iy];
}
delete [] A;
delete [] B;
delete [] C;
}

#endif