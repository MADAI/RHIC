#ifndef __PCA_CC__
#define __PCA_CC__
#include "pca.h"
using namespace std;

CPCA::CPCA(int nruns_set){
	string filename;
	string segment,dumbstring;
	char dummy[100],type[30];
	double dummyvalue;
	FILE *fptr;
	nruns=nruns_set;
	printf("nruns=%d\n",nruns);
	// qualifiers.Read("qualifiers.dat");
	// printf("qualifiers read\n");
	int iy,jy,iqual,iname;
	nnames=0;
	//for(iqual=0;iqual<qualifiers.nqualifiers;iqual++){
	fptr=fopen("pcanames.dat","r");
	while(!feof(fptr)){
		fscanf(fptr,"%s",dummy);
		if(dummy[0]!='#'){
			pcaname[nnames]=dummy;
			printf("pca name=%s\n",dummy);
			nnames+=1;
		}
	}
	fclose(fptr);
	//}
	printf("nnames=%d\n",nnames);
	ny=0;
	//for(iqual=0;iqual<qualifiers.nqualifiers;iqual++){
	filename="model_results/run1/results.dat";
	printf("reading %s\n",filename.c_str());
	fptr=fopen(filename.c_str(),"r");
	do{
		fscanf(fptr,"%s",type);
		if(!feof(fptr)){
			if(type[0]!='#'){
				fscanf(fptr,"%s",dummy);
				//printf("varname=%s\n",dummy);
				if(namecheck(dummy)){
					yname[ny]=string(dummy);
					fscanf(fptr,"%lf",&dummyvalue);
					fscanf(fptr,"%lf",&sigmay[ny]);
					dumbstring=dummy;
					segment=dumbstring.substr(dumbstring.length()-6,dumbstring.length()-1);
					//printf("segment=%s\n",segment.c_str());
					if(segment=="_YIELD") sigmay[ny]=0.075*dummyvalue;
					if(segment=="MEANPT")	sigmay[ny]=0.075*dummyvalue;
					if(segment=="PCAVAL") sigmay[ny]=0.1*dummyvalue;
					segment=dumbstring.substr(dumbstring.length()-8,dumbstring.length()-1);
					if(segment=="PTWEIGHT") sigmay[ny]=0.01+0.05*dummyvalue;
					segment=dumbstring.substr(dumbstring.length()-5,dumbstring.length()-1);
					if(segment=="_PION") sigmay[ny]=0.2+0.05*dummyvalue;
					//printf("yname[%d]=%s\n",ny,yname[ny].c_str());
					segment=dumbstring.substr(0,23);
					if(segment=="cent20to30_STAR_PION_V2") sigmay[ny]=0.05;//*dummyvalue;//0.07*dummyvalue;
					//printf("segment=%s\n",segment.c_str());
					ny+=1;
				}
				else fgets(dummy,100,fptr);
			}
			else fgets(dummy,100,fptr);
		}
	}while(!feof(fptr));
	fclose(fptr);
	printf("_______ ny=%d___________\n",ny);

	//}
	ybar=new double[ny];
	value=new double[ny];
	spread=new double*[ny];
	for(iy=0;iy<ny;iy++){
		ybar[iy]=0.0;
		spread[iy]=new double[ny];
		for(jy=0;jy<ny;jy++) spread[iy][jy]=0.0;
	}
	gslmatrix=new CGSLMatrix_Real(ny);
	//printf("CPCA initialized\n");
}

void CPCA::ReadResults(){
	FILE *fptr;
	string filename,varname;
	char type[30],dummy[100];
	double dummysigma,valsum;
	int iy,jy,irun,iqual;
	for(irun=0;irun<nruns;irun++){
		iy=0;
		//for(iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		sprintf(dummy,"%d",irun+1);
		filename="model_results/run"+string(dummy)+"/results.dat";
		//printf("reading %s\n",filename.c_str());
		fptr=fopen(filename.c_str(),"r");
		valsum=0.0;
		do{
			fscanf(fptr,"%s",type);
			if(!feof(fptr)){
				if(type[0]!='#'){
					fscanf(fptr,"%s",dummy);
					varname=dummy;
					if(namecheck(varname)){
						fscanf(fptr,"%lf",&value[iy]);
						fscanf(fptr,"%lf",&dummysigma);
						dummysigma=0.075*value[iy];
						if(irun==0) printf("%s[%d]= %g+/-%g\n",yname[iy].c_str(),irun,value[iy],sigmay[iy]);
						//if(fabs(value[iy])>1.0) Misc::Pause();
						ybar[iy]+=value[iy];
						valsum+=value[iy];
						for(jy=0;jy<=iy;jy++){
							spread[iy][jy]+=value[iy]*value[jy];
						}
						iy+=1;
					}
					else fgets(dummy,100,fptr);
				}
				else fgets(dummy,100,fptr);
			}
			if(iy>ny){
				printf("iy=%d, but ny=%d\n",iy,ny);
				exit(1);
			}
		}while(!feof(fptr));
		fclose(fptr);
		//printf("valsum=%g\n",valsum);
		//}
	}
	for(iy=0;iy<ny;iy++){
		ybar[iy]=ybar[iy]/double(nruns);
		for(jy=0;jy<=iy;jy++){
			spread[iy][jy]=spread[iy][jy]/double(nruns);
			spread[iy][jy]=spread[iy][jy]-ybar[iy]*ybar[jy];
			printf("spread[%d][%d]=%g\n",iy,jy,spread[iy][jy]);
			spread[iy][jy]=spread[iy][jy]/(sigmay[iy]*sigmay[jy]);
			if(iy!=jy) spread[jy][iy]=spread[iy][jy];
			//printf("%10.3e ",spread[iy][jy]);
		}
		//printf("\n");
	}
	double spreadtest=0.0;
	for(iy=0;iy<ny;iy++){
		for(jy=0;jy<ny;jy++) spreadtest+=spread[iy][jy];
	}
	printf("spreadtest=%g, ny=%d\n",spreadtest,ny);

	/* 
	printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	for(iy=0;iy<ny;iy++) printf("ybar[%d]=%g, name=%s\n",iy,ybar[iy],yname[iy].c_str());
	printf("--------------- SPREAD ----------------------\n");
	for(iy=0;iy<ny;iy++){
	for(jy=0;jy<ny;jy++) printf("%8.4f ",spread[iy][jy]);
	printf("\n");
	}
	printf("---------------------------------------------------------\n");
	*/
}

void CPCA::Calc(){
	int iy,jy;
	printf("CPCA::Calc(), ny=%d\n",ny);
	eigenval=new double[ny];
	evec=new double*[ny];
	for(iy=0;iy<ny;iy++) evec[iy]=new double[ny];
	double trace=0.0;
	for(iy=0;iy<ny;iy++) trace+=spread[iy][iy];

	CGSLMatrix_Real *gslmatrix=new CGSLMatrix_Real(ny);
	gslmatrix->EigenFind(spread,evec,eigenval);
	double etrace=0.0;
	for(iy=0;iy<ny;iy++){
		etrace+=eigenval[iy];
		printf("eigenval[%d]=%g\n",iy,eigenval[iy]);
	}
	printf("trace=%g=%g\n",trace,etrace);

	/*
	printf("0000000000000 checking 0000000000000000\n");
	double *yy=new double[ny];
	for(int ky=0;ky<ny;ky++){
	printf("------- checking eigenvector %d -----------------\n",ky);
	for(iy=0;iy<ny;iy++){
	yy[iy]=0.0;
	for(jy=0;jy<ny;jy++){
	yy[iy]+=spread[iy][jy]*evec[jy][ky];
	}
	printf("%g =? %g\n",yy[iy],eigenval[ky]*evec[iy][ky]);
	}
	}
	*/
	
	/*
	double *respower=new double[ny];
	for(iy=0;iy<ny;iy++){
	respower[iy]=0.0;
	for(jy=0;jy<ny;jy++) respower[iy]+=eigenval[jy]*pow(evec[iy][jy],2);
	printf("respower[%d]=%10.3e, %10.3e, %s\n",iy,respower[iy],spread[iy][iy],yname[iy].c_str());
	}
	delete [] respower;*/
	
	double sumcheck0=0,sumcheck1=0,sumcheck2=0,normcheck=0;
	double *enorm=new double[ny],*ynorm=new double[ny];
	for(jy=0;jy<ny;jy++){
		enorm[jy]=0.0;
		ynorm[jy]=0.0;
		for(iy=0;iy<ny;iy++){
			enorm[jy]+=pow(evec[iy][jy],2);
			ynorm[jy]+=pow(evec[iy][jy]*sigmay[iy],2);
		}
		enorm[jy]=sqrt(enorm[jy]);
		printf("enorm[%d]=%g\n",jy,enorm[jy]);
		ynorm[jy]=sqrt(ynorm[jy]);
	}
	printf("    ------------------- component --------------------   evec[0]  evec[1]  evec[2]  evec[3]\n");
	for(iy=0;iy<ny;iy++){
		printf("%3d %50s: %8.5f %8.5f %8.5f %8.5f\n", iy,yname[iy].c_str(),evec[iy][0]*sigmay[iy]/ynorm[iy],evec[iy][1]*sigmay[iy]/ynorm[iy],evec[iy][2]*sigmay[iy]/ynorm[iy],evec[iy][3]/ynorm[iy]*sigmay[iy]);
		sumcheck0+=evec[iy][0]*evec[iy][0];
		sumcheck1+=evec[iy][1]*evec[iy][1];
		sumcheck2+=evec[iy][2]*evec[iy][2];
		normcheck+=evec[iy][0]*evec[iy][1];
	}
	printf("sumchecks = %g, %g, %g, %g = ?0\n",sumcheck0,sumcheck1,sumcheck2,normcheck);
	normcheck=0.0;
	for(iy=0;iy<ny;iy++) normcheck+=evec[iy][0]*evec[iy][0];
	printf("normalizations=%g = ?1\n",normcheck);
	normcheck=0.0;
	for(iy=0;iy<ny;iy++) normcheck+=evec[iy][0]*evec[iy][0]/(enorm[0]*enorm[0]);
	printf("normalizations=%g = ?1\n",normcheck);
	normcheck=0.0;
	for(iy=0;iy<ny;iy++) normcheck+=(evec[iy][0]*evec[iy][0]*sigmay[iy]*sigmay[iy])/(ynorm[0]*ynorm[0]);
	printf("normalizations=%g = ?1\n",normcheck);

	delete [] enorm;
	delete [] ynorm;

}

bool CPCA::namecheck(string varname){
	bool answer=false;
	int iname=0;
	while(answer==false && iname<nnames){
		if(varname==pcaname[iname]) answer=true;
		iname+=1;
	}
	return answer;
}

#endif