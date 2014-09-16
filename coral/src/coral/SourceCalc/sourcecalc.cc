#ifndef INCLUDE_SOURCECALC_CC
#define INCLUDE_SOURCECALC_CC

#include "constants.h"
#include "sourcecalc.h"

using namespace std;

CSourceCalc::CSourceCalc(){
	randy=new CRandom(-1234);
}

void CSourceCalc::CalcS(CCHArray *A){
	printf("i'm a dummy CalcS(CCHArray)\n");
	// Dummy function to be over-written by inherited class
}

void CSourceCalc::CalcS(C3DArray *threed){
	printf("i'm a dummy CalcS(C3DArray)\n");
	// Dummy function to be over-written by inherited class
}

void CSourceCalc::GaussCFCalc(C3DArray *cf3d){
	printf("i'm a dummy CalcS(C3DArray *)\n");
	// Dummy function to be over-written by inherited class
}

void CSourceCalc::CalcS(int lx,int ly,int lz,CCHArray *A){
	printf("i'm a dummy CalcS(int,int,int,CCHArray *)\n");
	// Dummy function to be over-written by inherited class
}

void CSourceCalc::CalcS(CMCList *&lista,CMCList *&listb){
	printf("i'm a dummy CalcS(CMCList *a,CMCList *b)\n");
	// Dummy function to be over-written by inherited class
}

void CSourceCalc::CombineMCLists(CMCList *lista,CMCList *listb,CCHArray *A){
	A->ZeroArray();
	double rcm[4];
	long long int icalc=0,ncalc,icount=0;
	int ia,ib,nbmax,ir,na=lista->GetNMC(),nb=listb->GetNMC(),nrmax=A->GetNRADIAL();
	double *ra,*rb;
	double r,volume,delr=A->GetRADSTEP(),snorm;
	bool AEQUALB;
	if(lista==listb){
		AEQUALB=true;
		ncalc=na*(na-1)/2;
	}
	else{
		AEQUALB=false;
		ncalc=na*nb;
	}
	rcm[0]=0.0;
	printf("_______ In CombineMCLists: na=%d, nb=%d, AEQUALB=%d ______\n",na,nb,int(AEQUALB));
	A->ZeroArray();
	for(ia=0;ia<na;ia++){
		ra=lista->GetR(ia);
		if(AEQUALB) nbmax=ia-1;
		else nbmax=nb;
		for(ib=0;ib<=nbmax;ib++){
			rb=listb->GetR(ib);
			rcm[1]=ra[1]-rb[1];
			rcm[2]=ra[2]-rb[2];
			rcm[3]=ra[3]-rb[3];
			r=sqrt(rcm[1]*rcm[1]+rcm[2]*rcm[2]+rcm[3]*rcm[3]);
			if(r<nrmax*delr){
				ir=int(floor(r/delr));
				if(ir<nrmax){
					A->IncrementAExpArrayFromE(rcm[1]/r,rcm[2]/r,rcm[3]/r,1.0,ir);
				}
			}
			icalc+=1;
			icount+=1;
			if(icount*10>=ncalc){
				printf("finished %4.1f percent\n",100.0*double(icalc)/double(ncalc));
				icount=0;
			}
		}
	}
	A->FillRemainderX();
	
	if(AEQUALB) snorm=2.0/(double(na)*double(na-1));
	else snorm=1.0/(double(na)*double(nb));
	
	for(ir=0;ir<nrmax;ir++){
		volume=(4.0*PI/3)*(pow((ir+1)*delr,3)-pow(double(ir)*delr,3));
		A->ScaleArray(snorm/volume,ir);
	}
}

void CSourceCalc::CombineMCPRLists(CMCPRList *lista,CMCPRList *listb,CCHArray *A){
	A->ZeroArray();
	double rcm[4];
	long long int icalc=0,npairs=0,icount=0,ncalc;
	int increment;
	int ia,ib,nbmax,ir,na=lista->GetNMC(),nb=listb->GetNMC(),nrmax=A->GetNRADIAL();
	double *ra,*rb,*pa,*pb;
	double phia,phib,deltaphi;
	double r,volume,delr=A->GetRADSTEP(),snorm;
	double DELTAPHI_DEG=parameter::getD(spars,"DELTAPHI_DEG",30.0);
	bool XREFLECTIONSYMMETRY=parameter::getB(spars,"XREFLECTIONSYMMETRY",false);
	bool YREFLECTIONSYMMETRY=parameter::getB(spars,"YREFLECTIONSYMMETRY",false);
	bool AEQUALB;
	if(lista==listb){
		AEQUALB=true;
		ncalc=na*(na-1)/2;
	}
	else{
		AEQUALB=false;
		ncalc=na*nb;
	}
	rcm[0]=0.0;
	printf("_______ In CombineMCLists: na=%d, nb=%d, AEQUALB=%d ______\n",na,nb,int(AEQUALB));
	A->ZeroArray();
	for(ia=0;ia<na;ia++){
		ra=lista->GetR(ia);
		pa=lista->GetP(ia);
		phia=atan2(pa[2],pa[1]);
		if(AEQUALB) nbmax=ia-1;
		else nbmax=nb;
		for(ib=0;ib<=nbmax;ib++){
			rb=listb->GetR(ib);
			pb=listb->GetP(ib);
			//
			increment=0;
			phib=atan2(pb[2],pb[1]);
			deltaphi=(180.0/PI)*fabs(phia-phib);
			if(deltaphi>180.0) deltaphi=fabs(deltaphi-360.0);
			if(deltaphi<DELTAPHI_DEG){
				increment=1;
			}
			/*
			 if(XREFLECTIONSYMMETRY){
			 phib=atan2(pb[2],-pb[1]);
			 deltaphi=(180.0/PI)*fabs(phia-phib);
			 if(deltaphi>180.0) deltaphi=fabs(deltaphi-360.0);
			 if(deltaphi<DELTAPHI_DEG){
			 increment+=1;
			 }
			 }
			 if(YREFLECTIONSYMMETRY){
			 phib=atan2(-pb[2],pb[1]);
			 deltaphi=(180.0/PI)*fabs(phia-phib);
			 if(deltaphi>180.0) deltaphi=fabs(deltaphi-360.0);
			 if(deltaphi<DELTAPHI_DEG){
			 increment+=1;
			 }
			 }
			 if(XREFLECTIONSYMMETRY && YREFLECTIONSYMMETRY){
			 phib=atan2(-pb[2],-pb[1]);
			 deltaphi=(180.0/PI)*fabs(phia-phib);
			 if(deltaphi>180.0) deltaphi=fabs(deltaphi-360.0);
			 if(deltaphi<DELTAPHI_DEG){
			 increment+=1;
			 }
			 }
			 */
			if(increment!=0){
				npairs+=increment;
				rcm[1]=ra[1]-rb[1];
				rcm[2]=ra[2]-rb[2];
				rcm[3]=ra[3]-rb[3];
				r=sqrt(rcm[1]*rcm[1]+rcm[2]*rcm[2]+rcm[3]*rcm[3]);
				if(r<nrmax*delr){
					ir=int(lrint(floor(r/delr)));
					if(ir<nrmax){
						A->IncrementAExpArrayFromE(rcm[1]/r,rcm[2]/r,rcm[3]/r,double(increment),ir);
					}
				}
			}
			icalc+=1;
			icount+=1;
			if(icount*10>=ncalc){
				printf("finished %4.1f percent\n",100.0*double(icalc)/double(ncalc));
				icount=0;
			}
		}
	}
	A->FillRemainderX();
	
	snorm=1.0/double(npairs);
	for(ir=0;ir<nrmax;ir++){
		volume=(4.0*PI/3)*(pow((ir+1)*delr,3)-pow(double(ir)*delr,3));
		A->ScaleArray(snorm/volume,ir);
	}
}

void CSourceCalc::CombineMCLists(CMCList *lista,CMCList *listb,CCHArray *A,int NMC){
	A->ZeroArray();
	double rcm[4];
	int icalc=0,icount=0,ia,ib,imc,nbmax,ir,na=lista->GetNMC(),nb=listb->GetNMC(),nrmax=A->GetNRADIAL();
	double *ra,*rb;
	double r,volume,delr=A->GetRADSTEP(),snorm;
	bool AEQUALB;
	if(lista==listb){
		AEQUALB=true;
	}
	else{
		AEQUALB=false;
	}
	rcm[0]=0.0;
	printf("_______ In CombineMCLists: na=%d, nb=%d, AEQUALB=%d ______\n",na,nb,int(AEQUALB));
	A->ZeroArray();
	for(imc=0;imc<NMC;imc++){
		ia=lrint(floor(randy->ran()*na));
		do{
			ib=lrint(floor(randy->ran()*na));
		}while(AEQUALB && ia==ib);
		printf("ia=%d, ib=%d\n",ia,ib);
		ra=lista->GetR(ia);
		rb=listb->GetR(ib);
		rcm[1]=ra[1]-rb[1];
		rcm[2]=ra[2]-rb[2];
		rcm[3]=ra[3]-rb[3];
		r=sqrt(rcm[1]*rcm[1]+rcm[2]*rcm[2]+rcm[3]*rcm[3]);
		if(r<nrmax*delr){
			ir=int(floor(r/delr));
			if(ir<nrmax){
				A->IncrementAExpArrayFromE(rcm[1]/r,rcm[2]/r,rcm[3]/r,1.0,ir);
			}
		}
		icalc+=1;
		icount+=1;
		if(icount*10>=NMC){
			printf("finished %4.1f percent\n",100.0*double(icalc)/double(NMC));
			icount=0;
			
		}
	}
	A->FillRemainderX();
	snorm=1.0/double(NMC);
	
	for(ir=0;ir<nrmax;ir++){
		volume=(4.0*PI/3)*(pow((ir+1)*delr,3)-pow(double(ir)*delr,3));
		A->ScaleArray(snorm/volume,ir);
	}
}

void CSourceCalc::CombineMCPRLists(CMCPRList *lista,CMCPRList *listb,CCHArray *A,int NMC){
	A->ZeroArray();
	double rcm[4];
	double phia,phib,deltaphi;
	double DELTAPHI_DEG=parameter::getD(spars,"DELTAPHI_DEG",30.0);
	long long int icalc=0,icount=0,npairs=0;
	int ia,ib,imc,nbmax,ir,na=lista->GetNMC(),nb=listb->GetNMC(),nrmax=A->GetNRADIAL();
	double *ra,*rb,*pa,*pb;
	double r,volume,delr=A->GetRADSTEP(),snorm;
	bool AEQUALB;
	if(lista==listb){
		AEQUALB=true;
	}
	else{
		AEQUALB=false;
	}
	rcm[0]=0.0;
	printf("_______ In CombineMCLists: na=%d, nb=%d, AEQUALB=%d ______\n",na,nb,int(AEQUALB));
	A->ZeroArray();
	for(imc=0;imc<NMC;imc++){
		ia=lrint(floor(randy->ran()*na));
		do{
			ib=lrint(floor(randy->ran()*na));
		}while(AEQUALB && ia==ib);
		printf("ia=%d, ib=%d\n",ia,ib);
		ra=lista->GetR(ia);
		rb=listb->GetR(ib);
		pa=lista->GetP(ia);
		pb=listb->GetP(ib);
		phia=atan2(pa[2],pa[1]);
		phib=atan2(pb[2],pb[1]);
		deltaphi=(180.0/PI)*fabs(phia-phib);
		if(deltaphi<DELTAPHI_DEG){
			npairs+=1;
			rcm[1]=ra[1]-rb[1];
			rcm[2]=ra[2]-rb[2];
			rcm[3]=ra[3]-rb[3];
			r=sqrt(rcm[1]*rcm[1]+rcm[2]*rcm[2]+rcm[3]*rcm[3]);
			if(r<nrmax*delr){
				ir=int(floor(r/delr));
				if(ir<nrmax){
					A->IncrementAExpArrayFromE(rcm[1]/r,rcm[2]/r,rcm[3]/r,1.0,ir);
				}
			}
		}
		icalc+=1;
		icount+=1;
		if(icount*10>=NMC){
			printf("finished %4.1f percent\n",100.0*double(icalc)/double(NMC));
			icount=0;
			
		}
	}
	A->FillRemainderX();
	snorm=1.0/double(npairs);
	
	for(ir=0;ir<nrmax;ir++){
		volume=(4.0*PI/3)*(pow((ir+1)*delr,3)-pow(double(ir)*delr,3));
		A->ScaleArray(snorm/volume,ir);
	}
}

void CSourceCalc::CombineMCLists(CMCList *lista,CMCList *listb,C3DArray *threed){
	threed->ZeroArray();
	double rcm[4];
	int ia,ib,nbmax,na=lista->GetNMC(),nb=listb->GetNMC();
	double *ra,*rb;
	double volume,snorm;
	bool AEQUALB=false;
	if(lista==listb) AEQUALB=true;
	rcm[0]=0.0;
	for(ia=0;ia<na;ia++){
		ra=lista->GetR(ia);
		if(AEQUALB) nbmax=ia-1;
		else nbmax=nb;
		for(ib=0;ib<nbmax;ib++){
			rb=listb->GetR(ib);
			rcm[1]=ra[1]-rb[1];
			rcm[2]=ra[2]-rb[2];
			rcm[3]=ra[3]-rb[3];
			threed->IncrementElement(rcm[1],rcm[2],rcm[3],1.0);	
		}
		if(10*(ia+1)%(10*int(na/10))==0)
			printf("finished %g percent\n",100*double(ia+1)/(10*int(na/10)));
	}
	
	if(AEQUALB) snorm=2.0/double(na*(na-1));
	else snorm=1.0/double(na*nb);
	
	volume=threed->GetDELX()*threed->GetDELY()*threed->GetDELZ();
	threed->ScaleArray(snorm/volume);
}

void CSourceCalc::CombineMCPRLists(CMCPRList *lista,CMCPRList *listb,C3DArray *threed){
	threed->ZeroArray();
	double DELTAPHI_DEG=parameter::getD(spars,"DELTAPHI_DEG",30.0);
	double rcm[4],phia,phib,deltaphi;
	int ia,ib,nbmax,na=lista->GetNMC(),nb=listb->GetNMC();
	long long int npairs=0;
	double *ra,*rb,*pa,*pb;
	double volume,snorm;
	bool AEQUALB=false;
	if(lista==listb) AEQUALB=true;
	rcm[0]=0.0;
	for(ia=0;ia<na;ia++){
		ra=lista->GetR(ia);
		pa=lista->GetP(ia);
		phia=atan2(pa[2],pa[1]);
		if(AEQUALB) nbmax=ia-1;
		else nbmax=nb;
		for(ib=0;ib<nbmax;ib++){
			rb=listb->GetR(ib);
			pb=listb->GetP(ib);
			phib=atan2(pb[2],pb[1]);
			deltaphi=(180.0/PI)*fabs(phia-phib);
			if(deltaphi<DELTAPHI_DEG){
				npairs+=1;
				rcm[1]=ra[1]-rb[1];
				rcm[2]=ra[2]-rb[2];
				rcm[3]=ra[3]-rb[3];
				threed->IncrementElement(rcm[1],rcm[2],rcm[3],1.0);
			}
		}
		if(10*(ia+1)%(10*int(na/10))==0)
			printf("finished %g percent\n",100*double(ia+1)/(10*int(na/10)));
	}
	snorm=1.0/double(npairs);
	
	volume=threed->GetDELX()*threed->GetDELY()*threed->GetDELZ();
	threed->ScaleArray(snorm/volume);
}

void CSourceCalc::ReadSPars(char *sparsfilename){
	parameter::ReadParsFromFile(spars,sparsfilename);
}

double CSourceCalc::GetNorm(CCHArray *A){
	double check,delr,volume;
	int ir,NRMAX;
	NRMAX=A->GetNRADIAL();
	delr=A->GetRADSTEP();
	check=0;
	for(ir=0;ir<NRMAX;ir++){
		volume=4.0*PI*(pow(delr*(ir+1.0),3)-pow(delr*ir,3))/3.0;
		check+=volume*A->GetElement(0,0,0,ir);
	}
	return check;
}

void CSourceCalc::NormCheck(CCHArray *A){
	double check=GetNorm(A);
	printf("normalization check = %g\n",check);
}

double CSourceCalc::GetNorm(C3DArray *threed){
	int nsx,nsy,nsz,isx,isy,isz,ix,iy,iz;
	int nxmax=threed->GetNXMAX();
	int nymax=threed->GetNYMAX();
	int nzmax=threed->GetNZMAX();
	double prefactor=threed->GetDELX()*threed->GetDELY()*threed->GetDELZ();
	double norm=0.0;
	nsx=nsy=nsz=2;
	if(threed->GetXSYM()) nsx=1;
	if(threed->GetYSYM()) nsy=1;
	if(threed->GetZSYM()) nsz=1;
	if(nsx==1) prefactor*=2.0;
	if(nsy==1) prefactor*=2.0;
	if(nsz==1) prefactor*=2.0;
	for(ix=0;ix<nxmax;ix++){
		for(iy=0;iy<nymax;iy++){
			for(iz=0;iz<nzmax;iz++){
				for(isz=0;isz<nsz;isz++){
					for(isy=0;isy<nsy;isy++){
						for(isx=0;isx<nsx;isx++){
							norm+=threed->GetElement(isx,ix,isy,iy,isz,iz)*prefactor;
						}
					}
				}
			}
		}
	}
	return norm;
}

void CSourceCalc::NormCheck(C3DArray *threed){
	double norm=GetNorm(threed);
	printf("Norm Check of 3DArray = %g\n",norm);
}


#endif
