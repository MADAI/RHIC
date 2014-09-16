#ifndef __B3DIO_CC__
#define __B3DIO_CC__
#include "b3d.h"

double CB3D::WriteOSCAR(int ievent){
	CB3DBinaryPartInfo bpart;
	char dummy[100];
	if(oscarfile==NULL){
		if(BINARY_RW)
			oscarfile=fopen(oscarfilename.c_str(),"wb");
		else{
			oscarfile=fopen(oscarfilename.c_str(),"w");
			fprintf(oscarfile,":OSCAR1997a\n");
			fprintf(oscarfile,"ipart -- id -- p[4] -- m -- x[4]\n");
			fprintf(oscarfile,"b3d output\n");
		}
	}
	int nparts=PartMap.size()+FinalPartMap.size();
	if(BINARY_RW){
		fwrite(&ievent,sizeof(int),1,oscarfile);
		fwrite(&nparts,sizeof(int),1,oscarfile);
	}
	else
		fprintf(oscarfile,"%7d %6d    %8.5f     %8.5f\n",ievent,nparts,parameter::getD(parmap,"GLAUBER_B",0.0),
	parameter::getD(parmap,"GLAUBER_B",0.0));
	double v,pperp,eperp,dnchdeta=0.0,dnchdy=0,t,twrite,tauwrite,etawrite,eta,deleta,y;
	FourVector rwrite,pwrite;
	double mass;
	int ipart,nmesons=0,nch=0;
	CPart *part;
	CPartMap::iterator ppos;
	ipart=0;
	ppos=FinalPartMap.begin();
	do{
		if(ppos==FinalPartMap.end())
			ppos=PartMap.begin();
		part=ppos->second;
		part->Propagate(part->tau_lastint);
		if(BJORKEN && fabs(part->eta)>ETAMAX){
			part->CyclicReset();
		}
		if(part->resinfo->baryon!=0)
			nbaryons+=1;
		if(BINARY_RW){
			bpart.ID=part->resinfo->code;
			bpart.tau=part->tau0;
			bpart.x=part->r[1];
			bpart.y=part->r[2];
			bpart.eta=part->eta;
			bpart.px=part->p[1];
			bpart.py=part->p[2];
			bpart.rapidity=part->y;
			bpart.weight=part->weight;
			bpart.reality=part->reality;
			fwrite(&bpart,sizeof(bpart),1,oscarfile);
		}
		else
			fprintf(oscarfile,"%5d %5d %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %d %d\n",
		ipart,part->resinfo->code,part->p[1],part->p[2],part->p[3],part->p[0],sqrt(part->msquared),part->r[1],part->r[2],part->r[3],part->r[0],part->weight,int(part->reality));
		if(ppos==PartMap.end()){
			printf("ppos shouldn't be here\n");
		}
		++ppos;
		ipart+=1;
	} while(ipart<nparts);
	return dnchdy/(2.0*ETAMAX);
}

void CB3D::ReadOSCARHeader(){
	int ndead=3,idead;
	char dummy[200];
	oscarfilename="model_output/"+run_name+"/"+qualifier+"/oscar.dat";
	if(BINARY_RW)
		oscarfile=fopen(oscarfilename.c_str(),"rb");
	else{
		oscarfile=fopen(oscarfilename.c_str(),"r");
		for(idead=0;idead<ndead;idead++)
			fgets(dummy,200,oscarfile);
	}
}

int CB3D::ReadOSCAR(int ievent){
	CB3DBinaryPartInfo bpart;
	CResInfo *resinfo;
	double p[4],r[4],mass,rapidity,eta,tau0;
	int weight,reality_int,ID;
	int nparts,nparts_read,ipart=0;
	int ievent_read;
	bool reality;
	double bmin,bmax; // impact parameter
	double mtot,mothermass;
	CPart *mother;
	tau=0.0;
	if(oscarfile==NULL){
		ReadOSCARHeader();
	}
	if(BINARY_RW){
		fread(&ievent_read,sizeof(int),1,oscarfile);
		fread(&nparts_read,sizeof(int),1,oscarfile);
	}
	else{
		fscanf(oscarfile,"%d %d %lf %lf",&ievent_read,&nparts_read,&bmin,&bmax);
		if(!feof(oscarfile) && ievent_read!=ievent){
			printf("trying to read wrong event, ievent=%d, ievent_read=%d\n",ievent,ievent_read);
			exit(1);
		}
	}
	if(feof(oscarfile)){
		return 0;
	}
	for(ipart=0;ipart<nparts_read;ipart++){
		mother=GetDeadPart();
		if(BINARY_RW){
			fread(&bpart,sizeof(bpart),1,oscarfile);
			ID=bpart.ID;
			tau0=bpart.tau;
			r[1]=bpart.x;
			r[2]=bpart.y;
			eta=bpart.eta;
			p[1]=bpart.px;
			p[2]=bpart.py;
			rapidity=bpart.rapidity;
			resinfo=reslist->GetResInfoPtr(ID);
			mass=resinfo->mass;
			weight=bpart.weight;
			reality=bpart.reality;
		}
		else{
			fscanf(oscarfile,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d",
			&ipart,&ID,&p[1],&p[2],&p[3],&p[0],&mass,&r[1],&r[2],&r[3],&r[0],&weight,&reality_int);
			tau0=sqrt(r[0]*r[0]-r[3]*r[3]);
			eta=asinh(r[3]/tau0);
			rapidity=asinh(p[3]/p[0]);
			reality=true;
			if(reality_int==0)
				reality=false;
		}
		mother->Init(ID,r[1],r[2],tau0,eta,p[1],p[2],mass,rapidity,weight,reality);
	}
	return nparts_read;
}

void CB3D::WriteDens(){
	string densfilename="model_output/"+run_name+"/"+qualifier+"/dens.dat";
	FILE *densfile = fopen(densfilename.c_str(),"w");
	fprintf(densfile,"#ix iy  dens[itau=0] dens[itau=1]...\n");
	double dxy;
	int ix,iy,ieta,itau;
	for(ix=0;ix<2*NXY;ix++){
		for(iy=0;iy<2*NXY;iy++){
			fprintf(densfile,"%3d %3d",ix,iy);
			for(itau=0;itau<DENSWRITE_NTAU;itau++){
				dxy=0.0;
				for(ieta=0;ieta<2*NETA;ieta++){
					dxy+=cell[ix][iy][ieta]->dens[itau];
				}
				fprintf(densfile," %6.0f",dxy);
			}
			fprintf(densfile,"\n");
		}
	}
	fclose(densfile);
}

#endif
