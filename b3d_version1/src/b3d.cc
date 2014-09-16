#ifndef __B3D_CC__
#define __B3D_CC__

#include "b3d.h"
#include "hydrotob3d.h"
#ifdef __B3D_USE_HDF5__
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif
#endif
using namespace std;

CB3D::CB3D(){
	randy=new CRandom(-1234);
	BJORKEN=false;
};

CB3D::CB3D(string run_name_set){
	run_name=run_name_set;
	string parsfilename,dirname;
	dirname="parameters/"+run_name;
	parsfilename=dirname+"/fixed.param";
	printf("reading %s\n",parsfilename.c_str());
	parameter::ReadParsFromFile(parmap,parsfilename);
	parsfilename=dirname+"/stats.param";
	printf("reading %s\n",parsfilename.c_str());
	parameter::ReadParsFromFile(parmap,parsfilename);
	NACTIONSMAX=parameter::getI(parmap,"B3D_NACTIONSMAX",5000);
	NPARTSMAX=parameter::getI(parmap,"B3D_NPARTSMAX",2000);
	TAUCOLLMAX=parameter::getD(parmap," B3D_TAUCOLLMAX",50.0);
	DENSWRITE=parameter::getB(parmap,"B3D_DENSWRITE",false);
	DENSWRITE_NTAU=parameter::getI(parmap,"B3D_DENSWRITE_NTAU",20);
	DENSWRITE_DELTAU=parameter::getD(parmap,"B3D_DENSWRITE_DELTAU",1.0);
	input_dataroot=parameter::getS(parmap,"B3D_INPUT_DATAROOT","data/b3d");
	output_dataroot=parameter::getS(parmap,"B3D_OUTPUT_DATAROOT","data/b3d");
	//NxExVxExNxTxS=parameter::getI(parmap,"B3D_NACTIONS",1);
	NSAMPLE=parameter::getI(parmap,"B3D_NSAMPLE",1);
	ERROR_PRINT=parameter::getB(parmap,"B3D_ERROR_PRINT",true);
	XYMAX=parameter::getD(parmap,"B3D_XYMAX",15);
	ETAMAX=parameter::getD(parmap,"B3D_ETAMAX",1.0);
	NETA=parameter::getI(parmap,"B3D_NETA",10);
	NXY=parameter::getI(parmap,"B3D_NXY",10);
	SIGMAMAX=parameter::getD(parmap,"B3D_SIGMAMAX",30);
	NRINGSMAX=parameter::getI(parmap,"B3D_NRINGSMAX",200);
	NPRCELLSMAX=parameter::getI(parmap,"B3D_PR_NPRCELLSMAX",100000);
	COLLISIONS=parameter::getB(parmap,"B3D_COLLISIONS",true);
	BJORKEN=parameter::getB(parmap,"B3D_BJORKEN",false);
	SIGMADEFAULT=parameter::getD(parmap,"B3D_SIGMADEFAULT",1.5);
	SIGMAINELASTIC=parameter::getD(parmap, "B3D_SIGMAINELASTIC", 7);
	INELASTIC=parameter::getB(parmap, "B3D_INELASTIC",true);
	Q0=parameter::getD(parmap, "B3D_INELASTIC_Q0", 0);
	OSUHYDRO=parameter::getB(parmap,"B3D_OSUHYDRO",false);
	PRHYDRO=parameter::getB(parmap,"B3D_PRHYDRO",false);
	HYDRO_OCTANT_SYMMETRY=parameter::getI(parmap,"HYDRO_OCTANT",2);
	HYDRO_PURE_BJORKEN=parameter::getB(parmap,"HYDRO_PURE_BJORKEN",false);
	ANNIHILATION_CHECK=parameter::getB(parmap,"B3D_ANNIHILATION_CHECK",false);
	ANNIHILATION_SREDUCTION=parameter::getD(parmap,"B3D_ANNIHILATION_SREDUCTION",1.0);

	SIGMAMAX=SIGMAMAX/double(NSAMPLE);
	NPARTSMAX*=NSAMPLE;
	NACTIONSMAX*=NSAMPLE;
	string command="mkdir -p output/"+run_name;
	system(command.c_str());
	double xmin,xmax,ymin,ymax,etamin,etamax;
	DXY=XYMAX/double(NXY);
	DETA=ETAMAX/double(NETA);
	int jx,jy,jeta;
	CB3DCell *c;
	CB3DCell::b3d=this;
	CPart::b3d=this;
	if(OSUHYDRO){
		osuhydrotob3d=new COSUHydrotoB3D();
		osuhydrotob3d->b3d=this;
		osuhydrotob3d->initialization=false;
	}
	else{
		hydrotob3d=new CHYDROtoB3D();
		hydrotob3d->b3d=this;
		hydrotob3d->initialization=false;
	}
	bjmaker.b3d=this;
	CResList::b3d=this;
	CInelasticList::b3d=this;
	CInelasticList::UseFile = false;
	CInelasticList::UseInelasticArray = false;
	tau=0.0;
	itau=0;
	DeadPartMap.clear();
	PartMap.clear();
	FinalPartMap.clear();
	ActionMap.clear();
	DeadActionMap.clear();
	cell=new CB3DCell***[2*NXY];
	int ix,iy,ieta;
	for(ix=0;ix<2*NXY;ix++){
		xmin=-XYMAX+ix*DXY;
		xmax=xmin+DXY;
		cell[ix]=new CB3DCell**[2*NXY];
		for(iy=0;iy<2*NXY;iy++){
			ymin=-XYMAX+iy*DXY;
			ymax=ymin+DXY;
			cell[ix][iy]=new CB3DCell*[2*NETA];
			for(ieta=0;ieta<2*NETA;ieta++){
				etamin=-ETAMAX+DETA*ieta;
				etamax=etamin+DETA;
				cell[ix][iy][ieta]=new CB3DCell(xmin,xmax,ymin,ymax,etamin,etamax);
			}
		}
	}
	for(ix=0;ix<2*NXY;ix++){
		for(iy=0;iy<2*NXY;iy++){
			for(ieta=0;ieta<2*NETA;ieta++){
				c=cell[ix][iy][ieta];
				c->ix=ix; c->iy=iy; c->ieta=ieta;
				c->creflection=NULL;
				for(jx=ix-1;jx<=ix+1;jx++){
					for(jy=iy-1;jy<=iy+1;jy++){
						for(jeta=ieta-1;jeta<=ieta+1;jeta++){
							if(jx>=0 && jy>=0 && jeta>=0 && jx<2*NXY && jy<2*NXY && jeta<2*NETA){
								c->neighbor[1+jx-ix][1+jy-iy][1+jeta-ieta]=cell[jx][jy][jeta];
							}
							else c->neighbor[1+jx-ix][1+jy-iy][1+jeta-ieta]=NULL;
						}
					}
				}
			}
		}
	}
	if(BJORKEN){
		for(ix=0;ix<2*NXY;ix++){
			for(iy=0;iy<2*NXY;iy++){
				c=cell[ix][iy][0];
				c->ireflection=-1;
				c->creflection=cell[ix][iy][2*NETA-1];
				for(jx=ix-1;jx<=ix+1;jx++){
					for(jy=iy-1;jy<=iy+1;jy++){
						if(jx>=0 && jy>=0 && jx<2*NXY && jy<2*NXY){
							c->neighbor[1+jx-ix][1+jy-iy][0]=cell[jx][jy][2*NETA-1];
						}
					}
				}				
				c=cell[ix][iy][2*NETA-1];
				c->ireflection=1;
				c->creflection=cell[ix][iy][0];
				for(jx=ix-1;jx<=ix+1;jx++){
					for(jy=iy-1;jy<=iy+1;jy++){
						if(jx>=0 && jy>=0 && jx<2*NXY && jy<2*NXY){
							c->neighbor[1+jx-ix][1+jy-iy][2]=cell[jx][jy][0];
						}
					}
				}
			}
		}
	}
	randy=new CRandom(-1234);
	// create particle objects, and put into deadpartmap
	int ipart;
	partarray=new CPart *[NPARTSMAX];
	for(ipart=0;ipart<NPARTSMAX;ipart++){
		partarray[ipart]=new CPart();
	}
	// create array of action objects
	int iaction;
	CAction::b3d=this;
	actionarray=new CAction *[NACTIONSMAX];
	for(iaction=0;iaction<NACTIONSMAX;iaction++){
		actionarray[iaction]=new CAction();
	}

	reslist=new CResList();
	if(INELASTIC) inelasticlist = new CInelasticList();
	oscarfile=NULL;

#ifdef __B3D_USE_HDF5__
	h5infile=NULL;
	h5outfile=NULL;
	
	ptype=new CompType(sizeof(CPartH5));
	ptype->insertMember("listid", HOFFSET(CPartH5,listid), PredType::NATIVE_INT);
	ptype->insertMember("ID", HOFFSET(CPartH5,ID), PredType::NATIVE_INT);
	ptype->insertMember("x", HOFFSET(CPartH5,x), PredType::NATIVE_DOUBLE);
	ptype->insertMember("y", HOFFSET(CPartH5,y), PredType::NATIVE_DOUBLE);
	ptype->insertMember("eta", HOFFSET(CPartH5,eta), PredType::NATIVE_DOUBLE);
	ptype->insertMember("tau", HOFFSET(CPartH5,tau), PredType::NATIVE_DOUBLE);
	ptype->insertMember("px", HOFFSET(CPartH5,px), PredType::NATIVE_DOUBLE);
	ptype->insertMember("py", HOFFSET(CPartH5,py), PredType::NATIVE_DOUBLE);
	ptype->insertMember("rapidity", HOFFSET(CPartH5,rapidity), PredType::NATIVE_DOUBLE);
	ptype->insertMember("mass", HOFFSET(CPartH5,mass), PredType::NATIVE_DOUBLE);
#endif

	if(DENSWRITE){
		CB3DCell *cptr;
		int itau,ix,iy,ieta;
		for(ix=0;ix<2*NXY;ix++){
			for(iy=0;iy<2*NXY;iy++){
				for(ieta=0;ieta<2*NETA;ieta++){
					cptr=cell[ix][iy][ieta];
					cptr->dens=new double[DENSWRITE_NTAU];
					for(itau=0;itau<DENSWRITE_NTAU;itau++){
						cptr->dens[itau]=0.0;
					}
				}
			}
		}
	}
	InitArrays();
	if(ANNIHILATION_CHECK){
		int iarray,imax;
		imax=lrint(TAUCOLLMAX);
		annihilation_array=new double[imax];
		for(int i=0;i<imax;i++){
			annihilation_array[i]=0.0;
		}
	}
}

void CB3D::InitArrays(){
	int ipart,alpha;
	PartMap.clear();
	DeadPartMap.clear();
	FinalPartMap.clear();
	for(ipart=0;ipart<NPARTSMAX;ipart++){
		partarray[ipart]->listid=ipart;
		partarray[ipart]->key=ipart;
		partarray[ipart]->tau0=0.0;
		partarray[ipart]->tauexit=0.0;
		partarray[ipart]->tau_lastint=0.0;
		partarray[ipart]->currentmap=&DeadPartMap;
		partarray[ipart]->cell=NULL;
		partarray[ipart]->actionmap.clear();
		partarray[ipart]->active=false;
		partarray[ipart]->taudecay=0.0;
		reslist->GetResInfoptr(211,partarray[ipart]->resinfo);
		for(alpha=0;alpha<4;alpha++){
			partarray[ipart]->r[alpha]=partarray[ipart]->p[alpha]=0.0;
		}
		partarray[ipart]->r[0]=1.0;
		partarray[ipart]->tau0=1.0;
		partarray[ipart]->p[0]=139.58;
		partarray[ipart]->y=partarray[ipart]->eta=0.0;
		DeadPartMap.insert(CPartPair(partarray[ipart]->key,partarray[ipart]));
	}
	// create array of action objects
	ActionMap.clear();
	DeadActionMap.clear();
	int iaction;
	for(iaction=0;iaction<NACTIONSMAX;iaction++){
		actionarray[iaction]->listid=iaction;
		actionarray[iaction]->key=iaction;
		actionarray[iaction]->tau=0.0;
		actionarray[iaction]->type=-1;
		actionarray[iaction]->currentmap=&DeadActionMap;
		DeadActionMap.insert(CActionPair(actionarray[iaction]->key,actionarray[iaction]));
	}
}

void CB3D::SetQualifier(string qualifier_set){
	ievent_write=ievent_read=0;
	qualifier=qualifier_set;
#ifdef __B3D_USE_HDF5__
	if(h5outfile!=NULL) delete h5outfile;
	if(h5infile!=NULL) delete h5infile;
	string command="mkdir -p output/"+run_name+"/"+qualifier;
	system(command.c_str());
	outfilename="output/"+run_name+"/"+qualifier+"/b3d.h5";
	h5outfile = new H5File(outfilename,H5F_ACC_TRUNC);
	h5infile=NULL;
#else
	if(oscarfile!=NULL){
		fclose(oscarfile);
		oscarfile=NULL;
	}
	//if(oscarinfile!=NULL) fclose(oscarinfile);
	string command="mkdir -p output/"+run_name+"/"+qualifier;
	system(command.c_str());
	oscarfilename="output/"+run_name+"/"+qualifier+"/oscar.dat";
	//oscarfile = fopen(oscarfilename.c_str(),"w");
	//oscarinfile=NULL;
#endif

}

#ifdef __B3D_USE_HDF5__
int CB3D::ReadDataH5(int ievent){
	KillAllActions();
	nactions=0;
	KillAllParts();
	if(h5infile==NULL){
		string infilename="output/"+run_name+"/"+qualifier+"/hydro.h5";
		h5infile = new H5File(infilename,H5F_ACC_RDONLY);
	}

	int nparts,ipart;
	char eventno[20];
	H5D_space_status_t status;
	sprintf(eventno,"event%d",ievent);
	hsize_t dim[1];
	DataSet *dataset = new DataSet (h5infile->openDataSet(eventno));
	//dataset->getSpaceStatus(status);
	//hsize_t datasetsize=dataset->getStorageSize();
	//printf("ievent=%d, status=%d, size=%d\n",ievent,int(status),int(datasetsize));
	DataSpace filespace = dataset->getSpace();
	int rank=filespace.getSimpleExtentDims(dim);
	nparts=dim[0];
	CPartH5 *partH5=new CPartH5[nparts];
	if(nparts>NPARTSMAX){
		printf("Increase NPARTSMAX, nparts=%d\n",nparts);
		exit(1);
	}
	dataset->read(partH5,*ptype);
	delete dataset;
	for(ipart=0;ipart<nparts;ipart++){
		partarray[ipart]->Init(partH5[ipart].ID,partH5[ipart].x,partH5[ipart].y,partH5[ipart].tau,partH5[ipart].eta,
			partH5[ipart].px,partH5[ipart].py,partH5[ipart].mass,partH5[ipart].rapidity);
	}
	delete [] partH5;
	//printf("READ IN %d PARTS\n",nparts);
	return nparts;
}
#endif

void CB3D::MovePartsToFinalMap(){
	CPart *part;
	CPartMap::iterator ppos,pppos;
	CB3DCell *c;
	int ix,iy,ieta,iaction;
	for(ix=0;ix<2*NXY;ix++){
		for(iy=0;iy<2*NXY;iy++){
			for(ieta=0;ieta<2*NETA;ieta++){
				c=cell[ix][iy][ieta];
				while(c->partmap.size()>0){
					//printf("partmap size=%d\n",int(c->partmap.size()));
					ppos=c->partmap.begin();
					part=ppos->second;
					if(part->cell!=c){
						printf("cells don't match\n");
					}
					if(part->resinfo->decay){
						printf("all particles should have decayed\n");
						part->Print();
						exit(1);
					}
					part->ChangeMap(&FinalPartMap);
					part->DeleteFromMap(&c->partmap);
					part->cell=NULL;
				}
				c->partmap.clear();
			}
		}
	}
	DeadActionMap.clear();
	for(iaction=0;iaction<NACTIONSMAX;iaction++){
		actionarray[iaction]->listid=iaction;
		actionarray[iaction]->key=iaction;
		DeadActionMap.insert(CActionPair(actionarray[iaction]->key,actionarray[iaction]));
	}
	if(DeadActionMap.size()!=NACTIONSMAX){
		printf("DeadActionMap wrong length = %d\n",int(DeadActionMap.size()));
	}
}

#ifdef __B3D_USE_HDF5__
double CB3D::WriteDataH5(){
	double v,pperp,eperp,e,dnchdeta=0.0,dnchdy=0,twrite,t,tauwrite,eta,etawrite,deleta,y;
	int ipart,nmesons=0,nch=0;
	CPart *part;
	CPartMap::iterator ppos;
	ievent_write+=1;
	int nparts=int(FinalPartMap.size());
	CPartH5 *partH5=new CPartH5[nparts];
	ppos=FinalPartMap.begin();
	ipart=0;
	while(ppos!=FinalPartMap.end()){
		part=ppos->second;
		if(abs(part->resinfo->baryon)!=0) nbaryons+=1;
		else nmesons+=1;
		if(abs(part->resinfo->charge)!=0) nch+=1;
		pperp=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
		eperp=sqrt(pperp*pperp+part->resinfo->mass*part->resinfo->mass);
		v=pperp/eperp;
		if(part->resinfo->charge!=0) dnchdy+=1.0;
		if(part->resinfo->charge!=0) dnchdeta+=v;
		tauwrite=part->tau_lastint;
		eta=part->eta;
		etawrite=part->GetEta(tauwrite);
		y=part->y;
		part->resinfo->count+=1;
		
		while(fabs(eta)>ETAMAX && BJORKEN){
			if(etawrite<-ETAMAX) deleta=2.0*ETAMAX*ceil((-ETAMAX-eta)/(2.0*ETAMAX));
			if(etawrite>ETAMAX) deleta=-2.0*ETAMAX*ceil((eta-ETAMAX)/(2.0*ETAMAX));
			eta+=deleta;
			etawrite+=deleta;
			y+=deleta;
		}
		e=eperp*cosh(y);
		
		twrite=tauwrite*cosh(etawrite);
		t=part->tau0*cosh(eta);
		
		partH5[ipart].listid=part->listid;
		partH5[ipart].ID=part->resinfo->code;
		partH5[ipart].x=part->r[1]+(part->p[1]/e)*(twrite-t);
		partH5[ipart].y=part->r[2]+(part->p[2]/e)*(twrite-t);
		partH5[ipart].tau=tauwrite;
		partH5[ipart].eta=etawrite;
		partH5[ipart].px=part->p[1];
		partH5[ipart].py=part->p[2];
		partH5[ipart].rapidity=y;
		partH5[ipart].mass=part->resinfo->mass;
		++ppos;
		ipart+=1;
		part->Kill();
	}
	if(int(DeadPartMap.size())!=NPARTSMAX){
		printf("some particles still out there\n");
		exit(1);
	}

	hsize_t dim[] = {nparts};   /* Dataspace dimensions */
	DataSpace space(1,dim );
	DataSet* dataset;
	char event_number[20];
	sprintf(event_number,"event%d",ievent_write);
	printf("writing  for event %s, nparts=%d\n",event_number,ipart);
	dataset = new DataSet(h5outfile->createDataSet(event_number,*ptype, space));
	dataset->write(partH5,*ptype);

	delete dataset;
	delete [] partH5;
	return dnchdeta/(2.0*ETAMAX);
}
#endif

double CB3D::WriteOSCAR(){
	ievent_write+=1;
	if(oscarfile==NULL){
		oscarfile=fopen(oscarfilename.c_str(),"w");
		fprintf(oscarfile,":OSCAR1997a\n");
		fprintf(oscarfile,"ipart -- id -- p[4] -- m -- x[4]\n");
		fprintf(oscarfile,"b3d output\n");
	}
	else{
		oscarfile=fopen(oscarfilename.c_str(),"a");
	}
	int nparts=PartMap.size()+FinalPartMap.size();
	fprintf(oscarfile,"%7d %6d    %8.5f     %8.5f\n",ievent_write,nparts,parameter::getD(parmap,"GLAUBER_B",0.0),parameter::getD(parmap,"GLAUBER_B",0.0));
	double v,pperp,eperp,dnchdeta=0.0,dnchdy=0,t,twrite,tauwrite,etawrite,eta,deleta,y,rwrite[4],pwrite[4],mass;
	int ipart,nmesons=0,nch=0;
	CPart *part;
	CPartMap::iterator ppos;
	ipart=0;
	ppos=FinalPartMap.begin();
	do{
		if(ppos==FinalPartMap.end()) ppos=PartMap.begin();
		part=ppos->second;
		t=part->tau0*cosh(part->tau0);
		tauwrite=part->tau_lastint;
		pperp=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
		eperp=sqrt(pperp*pperp+part->resinfo->mass*part->resinfo->mass);
		v=pperp/eperp;
		if(part->resinfo->charge!=0) dnchdy+=1.0;
		if(part->resinfo->charge!=0) dnchdeta+=v;
		eta=part->eta;
		etawrite=part->GetEta(tauwrite);
		y=part->y;
		while(fabs(etawrite)>ETAMAX && BJORKEN){
			if(etawrite<-ETAMAX) deleta=2.0*ETAMAX*ceil((-ETAMAX-etawrite)/(2.0*ETAMAX));
			if(etawrite>ETAMAX) deleta=-2.0*ETAMAX*ceil((etawrite-ETAMAX)/(2.0*ETAMAX));
			etawrite+=deleta;
			y+=deleta;
			eta+=deleta;
		}
		pwrite[1]=part->p[1];
		pwrite[2]=part->p[2];
		pwrite[3]=eperp*sinh(y);
		mass=part->GetMass();
		pwrite[0]=sqrt(mass*mass+pwrite[1]*pwrite[1]+pwrite[2]*pwrite[2]+pwrite[3]*pwrite[3]);

		t=part->tau0*cosh(eta);
		twrite=tauwrite*cosh(etawrite);
		rwrite[0]=twrite;
		rwrite[1]=part->r[1]+(pwrite[1]/pwrite[0])*(twrite-t);
		rwrite[2]=part->r[2]+(pwrite[2]/pwrite[0])*(twrite-t);
		rwrite[3]=tauwrite*sinh(etawrite);

		fprintf(oscarfile,"%5d %5d %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e %15.10e\n",
			ipart,part->resinfo->code,pwrite[1],pwrite[2],pwrite[3],pwrite[0],mass,rwrite[1],rwrite[2],rwrite[3],rwrite[0]);
		if(ppos==PartMap.end()){
			printf("ppos shouldn't be here\n");
		}
		++ppos;
		ipart+=1;
	} while(ipart<nparts);
	fclose(oscarfile);
	return dnchdy/(2.0*ETAMAX);
}

void CB3D::WriteDens(){
	string densfilename="output/"+run_name+"/"+qualifier+"/dens.dat";
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
void CB3D::WriteAnnihilationData(){
	if(ANNIHILATION_CHECK){
		double total=0.0;
		int itau,imax=lrint(TAUCOLLMAX);
		for(int itau=0;itau<imax;itau++){
			printf("%6.2f %g\n",(itau+0.5),annihilation_array[itau]);
			total+=annihilation_array[itau];
		}
		printf("%g total annihilations, nbaryons=%d, annihilation fraction=%g\n",total,nbaryons,2.0*total/double(nbaryons));
	}
}

void CB3D::PerformAllActions(){
	if(DENSWRITE){
		for(int itau=0;itau<DENSWRITE_NTAU;itau++){
			AddAction_DensCalc((itau+1.0)*DENSWRITE_DELTAU);
		}
		//PrintActionMap(&ActionMap);
	}
	CActionMap::iterator epos=ActionMap.begin();
	//cout << ActionMap.size() << " actions to perform" << endl;
	CAction *action;
	nscatter=ndecay=npass=nmerge=nswallow=npass=nexit=nactivate=ninelastic=ncheck=nactionkills=0;
	ncollisions=nannihilate=0;
	tau=0.0;
	nactions=0;	
	while(epos!=ActionMap.end()){
		action=epos->second;
		action->Perform();
		epos=ActionMap.begin();
	}
	MovePartsToFinalMap();
	/*
	printf("Actions Finished, nactions=%lld\n",nactions);
	printf("%d actions in ActionMap\n",int(ActionMap.size()));
	printf("%d actions in DeadActionMap\n",int(DeadActionMap.size()));
	printf("%d particles in FinalPartMap\n",int(FinalPartMap.size()));
	printf("%d particles in PartMap\n",int(PartMap.size()));
	printf("%d particles in DeadPartMap\n",int(DeadPartMap.size()));
	*/

}

void CB3D::KillAllActions(){
	CAction *action;
	CActionMap::iterator epos=ActionMap.begin();
	while(epos!=ActionMap.end()){
		action=epos->second;
		action->Kill();
		epos=ActionMap.begin();
	}
	nactions=0;
}

void CB3D::KillAllParts(){
	CPartMap *partmap=&PartMap;
	CPart *part;
	CPartMap::iterator ppos=partmap->begin();
	while(ppos!=partmap->end()){
		part=ppos->second;
		part->Kill();
		ppos=partmap->begin();
	}
	partmap=&FinalPartMap;
	ppos=partmap->begin();
	while(ppos!=partmap->end()){
		part=ppos->second;
		part->Kill();
		ppos=partmap->begin();
	}
	int ix,iy,ieta;
	for(ix=0;ix<2*NXY;ix++){
		for(iy=0;iy<2*NXY;iy++){
			for(ieta=0;ieta<2*NETA;ieta++){
				partmap=&(cell[ix][iy][ieta]->partmap);
				ppos=partmap->begin();
				while(ppos!=partmap->end()){
					part=ppos->second;
					part->Kill();
					ppos=partmap->begin();
				}
			}
		}
	}
}


void CB3D::PrintActionMap(CActionMap *actionmap){
	CActionMap::iterator epos;
	CAction *action;
	int iaction=0;
	printf("_________________ ACTIONMAP %d actions _________________________\n",int(actionmap->size()));
	for(epos=actionmap->begin();epos!=actionmap->end();++epos){
		iaction+=1;
		action=epos->second;
		printf("iaction=%d : ",iaction);
		action->Print();
	}
}

void CB3D::FindAllCollisions(){
	double taucoll;
	CPartMap::iterator ppos1,ppos2;
	CPart *part1,*part2;
	CActionMap::iterator epos;
	CAction *action;
	int nbefore=ActionMap.size();
	//printf("CB3D::FindAllCollisions, Resetting Collisions\n");

	for(ppos1=PartMap.begin();ppos1!=PartMap.end();++ppos1){
		part1=ppos1->second;
		part1->KillActions();
	}

	for(epos=ActionMap.begin();epos!=ActionMap.end();++epos){
		action=epos->second;
		if(action->type==2){
			printf("CB3D::FindAllCollisions, expected all type-2 actions to be dead\n");
			exit(1);
		}
	}

	ppos1=PartMap.begin();
	part1=ppos1->second;
	ppos2=ppos1; ++ppos2;
	while(ppos2!=PartMap.end()){
		part1=ppos1->second; part2=ppos2->second;
		FindCollision(part1,part2,taucoll);
		ppos1=ppos2;
		++ppos2;
	}
	int nafter=ActionMap.size();
	if(nbefore!=nafter){
		printf("CB3D::FindAllCollisions, nbefore=%d, nafter=%d\n",nbefore,nafter);
		exit(1);
		//Misc::Pause();
	}

}

void CB3D::PrintPartList(){
	CPartMap::iterator ppos2,ppos1=FinalPartMap.begin();
	while(ppos1!=FinalPartMap.end()){
		printf("%d ",ppos1->second->listid);
		ppos2=ppos1; ++ppos2;
		if(ppos2!=FinalPartMap.end()){
			if(ppos1->second->actionmother!=ppos2->second->actionmother) printf("| ");
		}
		++ppos1;
	}
	printf("\n");
}

void CB3D::ListFutureCollisions(){
	CActionMap::iterator epos=ActionMap.begin();
	CAction *action;
	CPartMap::iterator p1,p2;
	printf("------------------- LIST OF FUTURE COLLISIONS ---------------------\n");
	while(epos!=ActionMap.end()){
		action=epos->second;
		if(action->type==2){
			p1=action->partmap.begin();
			p2=p1; ++p2;
			printf("%d  %d  will collide at %g\n",p1->second->listid,p2->second->listid,double(action->tau));
		}
		epos++;
	}
}

double CB3D::GetPiBsquared(CPart *part1,CPart *part2){
	ofstream outputfile;
	
	int alpha;
	double r[4],P[4],q[4],Pdotq=0.0,Pdotr=0.0,P2=0.0,qdotr=0.0,q2=0.0,rsquared=0.0,g[4]={1,-1,-1,-1},y1,mt,eta1;
	double *p1=part1->p,*p2=part2->p,*r1=part1->r,*r2=part2->r;
	bool flip=false;
	if(BJORKEN && ((part1->cell->ieta==0 && part2->cell->ieta==2*NETA-1) || (part2->cell->ieta==0 && part1->cell->ieta==2*NETA-1))){
		flip=true;
		p1=new double[4];
		r1=new double[4];
		for(alpha=0;alpha<4;alpha++){
			p1[alpha]=part1->p[alpha];
			r1[alpha]=part1->r[alpha];
		}
		mt=sqrt(p1[0]*p1[0]-p1[3]*p1[3]);
		y1=part1->y;
		eta1=part1->eta;
		if(part1->cell->ieta==0){
			y1+=2.0*ETAMAX;
			eta1+=2.0*ETAMAX;
		}
		else{
			y1-=2.0*ETAMAX;
			eta1-=2.0*ETAMAX;
		}
		p1[0]=mt*cosh(y1);
		p1[3]=mt*sinh(y1);
		r1[0]=part1->tau0*cosh(eta1);
		r1[3]=part1->tau0*sinh(eta1);
		//printf("p1=(%g,%g,%g,%g), r1=(%g,%g,%g,%g)\n",p1[0],p1[1],p1[2],p1[3],r1[0],r1[1],r1[2],r1[3]);
		//printf("p2=(%g,%g,%g,%g), r2=(%g,%g,%g,%g)\n",p2[0],p2[1],p2[2],p2[3],r2[0],r2[1],r2[2],r2[3]);
	}
	for(alpha=0;alpha<4;alpha++){
		P[alpha]=p1[alpha]+p2[alpha];
		q[alpha]=p1[alpha]-p2[alpha];
		r[alpha]=r1[alpha]-r2[alpha];
		P2+=g[alpha]*P[alpha]*P[alpha];
		Pdotq+=g[alpha]*P[alpha]*q[alpha];
		Pdotr+=g[alpha]*P[alpha]*r[alpha];
	}
	for(alpha=0;alpha<4;alpha++){
		q[alpha]-=Pdotq*P[alpha]/P2;
		r[alpha]-=Pdotr*P[alpha]/P2;
		q2+=g[alpha]*q[alpha]*q[alpha];
		qdotr+=g[alpha]*q[alpha]*r[alpha];
	}
	for(alpha=0;alpha<4;alpha++){
		r[alpha]-=qdotr*q[alpha]/q2;
		rsquared-=g[alpha]*r[alpha]*r[alpha];
	}
	if(flip){
		delete [] p1;
		delete [] r1;
	}
	return PI*rsquared;
}

void CB3D::freegascalc_onespecies(double m,double t,double &p,double &e,double &dens,double &sigma2,double &dedt){
	const double prefactor=1.0/(2.0*PI*PI*pow(HBARC,3));
	double k0,k1,z,k0prime,k1prime,m2,m3,m4,t2,t3,I1,I2,Iomega;
	m2=m*m;
	m3=m2*m;
	m4=m2*m2;
	t2=t*t;
	t3=t2*t;
	z=m/t;
	if(z>1000.0){
		p=e=dens=dedt=0.0;
		printf("z is huge=%g, m=%g, t=%g\n",z,m,t);
	}
	else{
		if(z<0.0){
			printf("___z=%g,m=%g,T=%g ___\n",z,m,t);
			exit(1);
		}
		k0=Bessel::K0(z);
		k1=Bessel::K1(z);
		p=prefactor*(m2*t2*k0+2.0*m*t3*k1);
		e=prefactor*(3.0*m2*t2*k0+(m3*t+6.0*m*t3)*k1);
		dens=p/t;
		k0prime=-k1;
		k1prime=-k0-k1/z;
		dedt=prefactor*(6.0*m2*t*k0+(m3+18.0*m*t2)*k1-3.0*m3*k0prime-((m4/t)+6.0*m2*t)*k1prime);
		Iomega=exp(-m/t)/(30.0*PI*PI*HBARC*HBARC*HBARC);
		I1=pow(m,1.5)*pow(t,3.5)*7.5*sqrt(2.0*PI);
		I2=24.0*pow(t,5);
		sigma2=Iomega*(I1+I2+0.5*sqrt(I1*I2));  // this is an approximation (+/-2%) to messy integral
	}
}

double CB3D::CalcSigma(CPart *part1,CPart *part2){
	double sigma=0.0,Gamma,G,G2,MR,M,m1,m2,b,q2=0.0,q3=0.0,q4=0.0,qR2,tan2delta;
	double inel_d=0.0, q_prime;
	double sigma_elastic, sigma_merge, sigma_inelastic,temp;
	int ir1,ir2,irflip,alpha,G_Value,L_merge,netq, netb, nets, pmq=0, pmb=0, pms=0;
	CMerge *merge;
	list<CInelasticInfo>::iterator inel;
	list<CInelasticInfo> inel_list;
	bool G_Parity = false;
	double jR,j1,j2,j1_i,j2_i;
	
	if(part1->resinfo->G_Parity && part2->resinfo->G_Parity){
		//cout << "Using G Parity" << endl;
		G_Parity = true;
		G_Value = part1->resinfo->G_Parity * part2->resinfo->G_Parity;
	}
	if(!G_Parity){
		//cout << "Not using G Parity" << endl;
	}
	
	ir1=part1->resinfo->ires; ir2=part2->resinfo->ires;
	if(ir1>ir2){
		irflip=ir1; ir1=ir2; ir2=irflip;
	}
	//cout << "SIGMADEFAULT " << SIGMADEFAULT << " NSAMPLE " << NSAMPLE << endl;
	sigma=SIGMADEFAULT/double(NSAMPLE);
	sigma_elastic = sigma;

	merge=reslist->MergeArray[ir1][ir2];
	
	if(INELASTIC){
		if(inelasticlist->UseInelasticArray){
			inel_list = inelasticlist->InelasticArray[ir1][ir2];
			//cout << "Getting inelastic data from InelasticArray[" << ir1 << "][" << ir2 << "]" << endl;
		}else{
			netq = part1->resinfo->charge+part2->resinfo->charge;
			netb = part1->resinfo->baryon+part2->resinfo->baryon;
			nets = part1->resinfo->strange+part2->resinfo->strange;
			inel_list = inelasticlist->ThermalArray[abs(netq)][abs(nets)][abs(netb)][Misc::Sign(netq)][Misc::Sign(nets)][Misc::Sign(netb)];
			//cout << "Getting inelastic data from ThermalArray["<< abs(netq) << "][" << abs(nets) << "][" << abs(netb) << "][" << Misc::Sign(netq) \
			<< "][" << Misc::Sign(nets) << "][" << Misc::Sign(netb) << "]" << endl;
}
}
//	cout << "And it has " << inel_list.size() << " elements." << endl;


inel = inel_list.begin();

if(merge!=NULL || inel!=inel_list.end()){
	j1=part1->resinfo->spin;
	j2=part2->resinfo->spin;
	m1=part1->GetMass();
	m2=part2->GetMass();
	M=pow(part1->p[0]+part2->p[0],2);
	for(alpha=1;alpha<4;alpha++) M-=pow(part1->p[alpha]+part2->p[alpha],2);
		M=sqrt(M);
	if(M<m1+m2){
		printf("SCREWY MASSES: m1=%g, m2=%g, M=%g\n",m1,m2,M);
		part1->Print();
		part2->Print();
		q2=1.0E-6;
	}
	else q2=Misc::triangle(M,m1,m2);
}
if(INELASTIC){

	inel_d = (2.0*j1+1.0)*(2.0*j2+1.0)*q2;
		//inel_d = 0;


	while(inel!=inel_list.end()){
		if((G_Parity && (inel->resinfo_1->G_Parity * inel->resinfo_2->G_Parity == G_Value)) || (!G_Parity)){
			if(M > inel->min_mass){
				j1_i=inel->resinfo_1->spin;
				j2_i=inel->resinfo_2->spin;
				q_prime = Misc::triangle(M,inel->resinfo_1->mass, inel->resinfo_2->mass);
				temp = (2.0*j1_i+1.0)*(2.0*j2_i+1.0)*q_prime;
				inel_d += temp;
			}
		}
		inel++;
	}
	inel = inel_list.begin();
}

while(merge!=NULL){
	Gamma=merge->resinfo->width;
	b=merge->branching;
	jR=merge->resinfo->spin;
	MR=merge->resinfo->mass;
	L_merge = merge->L;

	
	if(MR > m1 + m2){
		qR2=Misc::triangle(MR,m1,m2);
		q3=pow(q2/qR2,(2*L_merge + 1)/2);
		q4=pow(q2/qR2,(2*L_merge)/2);
		G=Gamma*(MR/M)*q3*1.2/(1.0+0.2*q4);
		tan2delta=0.25*G*G/((M-MR)*(M-MR));
		sigma+=b*((4.0*PI*HBARC*HBARC/q2)*(tan2delta/(1.0+tan2delta))
			*((2.0*jR+1.0)/((2.0*j1+1.0)*(2.0*j2+1.0))))/NSAMPLE;

	}
	merge=merge->next;
}


sigma_merge = sigma-sigma_elastic;
if(INELASTIC){
	while(inel!=inel_list.end()){
		if((G_Parity && (inel->resinfo_1->G_Parity * inel->resinfo_2->G_Parity == G_Value)) || !G_Parity){
			if(M > inel->min_mass){
				j1_i=inel->resinfo_1->spin;
				j2_i=inel->resinfo_2->spin;

				sigma += SIGMAINELASTIC*(2.0*j1_i+1.0)*(2.0*j2_i+1.0)*Misc::triangle(M,inel->resinfo_1->mass,inel->resinfo_2->mass)/(inel_d + Q0);
			}
		}
		inel++;
	}
		//cout << "Done considering inelastic exit channels" << endl;
	sigma_inelastic = sigma-sigma_merge-sigma_elastic;
}	
	//cout << sigma_elastic << "\t" << sigma_merge << "\t" << sigma_inelastic << endl;
printf("%6.4f\t%6.4f\t%6.4f\t\n",sigma_elastic, sigma_merge, sigma_inelastic);
if(part1->r[0]!=part1->r[0] || part2->r[0]!=part2->r[0]){
	printf("NaN problem in CB3D::Collide");
	exit(1);
}
return 0;
}

void CB3D::ReadHydroInput(){
	if(OSUHYDRO) osuhydrotob3d->ReadInput();
	else hydrotob3d->ReadInput();
}

int CB3D::HydrotoB3D(){
	int n0;
	if(OSUHYDRO) n0=osuhydrotob3d->MakeEvent();
	else n0=hydrotob3d->MakeEvent();
	return n0;
}

void CB3D::Reset(){
	tau=0.0;
	KillAllActions();
	nactions=0;
	KillAllParts();
}

CB3D::~CB3D(){
#ifdef __B3D_USE_HDF5__
	delete h5outfile;
#endif
}

#endif
