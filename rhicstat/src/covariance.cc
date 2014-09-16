#ifndef __COVARIANCE_CC__
#define __COVARIANCE_CC__
#include "rhicstat.h"
using namespace std;

void CRHICStat::CalcCovariance(){
	int iz,ix,irun,jrun,NZ=NX,ia,NA=35;
	double xdiff,delx=0.1,zi=0.0,zj=0.0;
	double ***covarray=new double **[NZ];
	double ***covnorm=new double **[NZ];
	CZGetter_GP *zg;
	if(FIT_TYPE=="GP")
		zg=(CZGetter_GP *)zgetter;
	for(iz=0;iz<NZ;iz++){
		covarray[iz]=new double*[NX];
		covnorm[iz]=new double*[NX];
		for(ix=0;ix<NX;ix++){
			covarray[iz][ix]=new double[NA];
			covnorm[iz][ix]=new double[NA];
			for(ia=0;ia<NA;ia++){
				covarray[iz][ix][ia]=0.0;
				covnorm[iz][ix][ia]=0.0;
			}
		}
	}
	
	for(iz=0;iz<NZ;iz++){
		for(irun=1;irun<NRUNS;irun++){
			if(FIT_TYPE=="GP" && !(zg->LINEAROFF)){
				zi=runinfo[irun]->zlinear[iz];
			}
			for(jrun=0;jrun<irun;jrun++){
				if(FIT_TYPE=="GP" && !(zg->LINEAROFF)){
					zj=runinfo[jrun]->zlinear[iz];
				}
				for(ix=0;ix<NX;ix++){
					xdiff=fabs(runinfo[irun]->x[ix]-runinfo[jrun]->x[ix]);
					ia=lrint(floor(xdiff/delx));
					if(ia<NA){
						covarray[iz][ix][ia]+=(runinfo[irun]->z[iz]-zi)*(runinfo[jrun]->z[iz]-zj);
						covnorm[iz][ix][ia]+=1.0;
					}
				}
			}
		}
		for(ix=0;ix<NX;ix++){
			for(ia=0;ia<NA;ia++){
				covarray[iz][ix][ia]=covarray[iz][ix][ia]/covnorm[iz][ix][ia];
			}
		}
	}
	FILE *fptr;
	char filename[100];
	for(iz=0;iz<NZ;iz++){
		sprintf(filename,"statinfo/covariance/z%02d.dat",iz);
		fptr=fopen(filename,"w");
		printf("------ iz=%d ----------\n",iz);
		for(ia=0;ia<NA;ia++){
			printf("%5.3f ",(ia+0.5)*delx);
			fprintf(fptr,"%5.3f ",(ia+0.5)*delx);
			for(ix=0;ix<NX;ix++){
				printf("%7.3f ",covarray[iz][ix][ia]);
				fprintf(fptr,"%7.3f ",covarray[iz][ix][ia]);
			}
			printf("\n");
			fprintf(fptr,"\n");
		}
		fclose(fptr);
	}
	printf("-----------------------\n");

	for(iz=0;iz<NZ;iz++){
		for(ix=0;ix<NX;ix++){
			delete [] covarray[iz][ix];
			delete [] covnorm[iz][ix];
		}
		delete [] covarray[iz];
		delete [] covnorm[iz];
	}
	delete [] covarray;
}



#endif
