#ifndef __B3D_CC__
#define __B3D_CC__
#include <functional>
#include "b3d.h"

CB3D::CB3D(){
};

CB3D::CB3D(string run_name_set){
	run_name=run_name_set;
	string parsfilename,dirname;
	dirname="model_output/"+run_name;
	parsfilename=dirname+"/fixed_parameters.dat";
	printf("reading %s\n",parsfilename.c_str());
	parameter::ReadParsFromFile(parmap,parsfilename);
	parsfilename=dirname+"/parameters.dat";
	printf("reading %s\n",parsfilename.c_str());
	parameter::ReadParsFromFile(parmap,parsfilename);
	
	NACTIONSMAX=parameter::getI(parmap,"B3D_NACTIONSMAX",100000);
	NPARTSMAX=parameter::getI(parmap,"B3D_NPARTSMAX",200000);
	TAUCOLLMAX=parameter::getD(parmap,"B3D_TAUCOLLMAX",50.0);
	DENSWRITE=parameter::getB(parmap,"B3D_DENSWRITE",false);
	DENSWRITE_NTAU=parameter::getI(parmap,"B3D_DENSWRITE_NTAU",20);
	DENSWRITE_DELTAU=parameter::getD(parmap,"B3D_DENSWRITE_DELTAU",1.0);
	input_dataroot=parameter::getS(parmap,"B3D_INPUT_DATAROOT","data/b3d");
	output_dataroot=parameter::getS(parmap,"B3D_OUTPUT_DATAROOT","data/b3d");
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
	SIGMADEFAULT=parameter::getD(parmap,"B3D_SIGMADEFAULT",1.0);
	SIGMAINELASTIC=parameter::getD(parmap, "B3D_SIGMAINELASTIC",1.0);
	INELASTIC=parameter::getB(parmap, "B3D_INELASTIC",false);
	NBOSE=parameter::getI(parmap,"B3D_NBOSE",1);
	Q0=parameter::getD(parmap, "B3D_INELASTIC_Q0", 0);
	COOPERFRYE_CREATENEGPARTS=parameter::getB(parmap,"B3D_COOPERFRYE_CREATENEGPARTS","false");
	COOPERFRYE_WMAX=parameter::getI(parmap,"B3D_COOPERFRYE_WMAX",1);
	HYDRO_OCTANT_SYMMETRY=parameter::getI(parmap,"HYDRO_OCTANT",2);
	HYDRO_PURE_BJORKEN=parameter::getB(parmap,"HYDRO_PURE_BJORKEN",false);
	BARYON_ANNIHILATION=parameter::getB(parmap,"B3D_BARYON_ANNIHILATION",true);
	DELNPARTSTOT=parameter::getD(parmap,"B3D_DELNPARTSTOT",1000);
	DELNACTIONSTOT=parameter::getD(parmap,"B3D_DELNACTIONSTOT",2000);
	NSCATT_MAX=parameter::getI(parmap,"B3D_NSCATT_MAX",0);
	USE_OLD_SAMPLER=parameter::getB(parmap,"B3D_USE_OLD_SAMPLER",false);
	BINARY_RW=parameter::getB(parmap,"B3D_BINARY_RW",false);
	MUTCALC=parameter::getB(parmap,"B3D_MUTCALC",false);
	MUTCALC_DELTAU=parameter::getD(parmap,"B3D_MUTCALC_DELTAU",0.5);
	ANNIHILATION_SREDUCTION=parameter::getD(parmap,"B3D_ANNIHILATION_SREDUCTION",1.0);
	SECALC=parameter::getB(parmap,"B3D_SECALC",false);
	SIGMAMAX=SIGMAMAX/double(NSAMPLE);
	NPARTSMAX*=NSAMPLE;
	NACTIONSMAX*=NSAMPLE;
	npartstot=nactionstot=0;
	DXY=XYMAX/double(NXY);
	DETA=ETAMAX/double(NETA);

	ibalmax=0;
	BALANCE_DECAY=parameter::getB(parmap,"B3D_BALANCE_DECAY",false);
	BALANCE_CALC=parameter::getB(parmap,"B3D_BALANCE_CALC",false);
	
	CResList::b3d=this;
	CPart::b3d=this;
	tau=0.0;
	npartstot=nactionstot=0;
	PartMap.clear();
	DeadPartMap.clear();
	FinalPartMap.clear();
	ActionMap.clear();
	DeadActionMap.clear();
	randy=new CRandom(std::hash<string>()(run_name));
	CAction::b3d=this;
	oscarfile=NULL;
	sampler=NULL;
	reslist=new CResList();
}

void CB3D::InitCascade(){
	// First initialize cells
	int ix,iy,ieta,jx,jy,jeta,itau,iarray,imax;
	double xmin,xmax,ymin,ymax,etamin,etamax;
	CB3DCell *c;
	string command="mkdir -p model_output/"+run_name;
	system(command.c_str());
	CB3DCell::b3d=this;
	if(USE_OLD_SAMPLER){
		hydrotob3d=new CHYDROtoB3D();
		hydrotob3d->b3d=this;
		hydrotob3d->initialization=false;
	}
	else
		sampler=new Csampler(this);
	CInelasticList::b3d=this;
	CInelasticList::UseFile = false;
	CInelasticList::UseInelasticArray = false;
	if(INELASTIC)
		inelasticlist = new CInelasticList();
	
	cell.resize(2*NXY);
	for(ix=0;ix<2*NXY;ix++){
		xmin=-XYMAX+ix*DXY;
		xmax=xmin+DXY;
		cell[ix].resize(2*NXY);
		for(iy=0;iy<2*NXY;iy++){
			ymin=-XYMAX+iy*DXY;
			ymax=ymin+DXY;
			cell[ix][iy].resize(2*NETA);
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
	if(DENSWRITE){
		for(ix=0;ix<2*NXY;ix++){
			for(iy=0;iy<2*NXY;iy++){
				for(ieta=0;ieta<2*NETA;ieta++){
					c=cell[ix][iy][ieta];
					c->dens.resize(DENSWRITE_NTAU);
					for(itau=0;itau<DENSWRITE_NTAU;itau++){
						c->dens[itau]=0.0;
					}
				}
			}
		}
	}
	if(BARYON_ANNIHILATION){
		imax=lrint(TAUCOLLMAX);
		annihilation_array.resize(imax);
		for(int i=0;i<imax;i++){
			annihilation_array[i]=0.0;
		}
	}
	if(MUTCALC){
		CMuTInfo::b3d=this;
		CMuTInfo::DELTAU=MUTCALC_DELTAU;
		CMuTInfo::NTAU=TAUCOLLMAX/MUTCALC_DELTAU;
		CMuTInfo::NETEVENTS=0;
		CRegenerate::b3d=this;
		regen=new CRegenerate();
		muTinfo.resize(2*NXY);
		for(ix=0;ix<2*NXY;ix++)
			muTinfo[ix].resize(2*NXY);
		for(ix=0;ix<2*NXY;ix++){
			for(iy=0;iy<2*NXY;iy++)
				muTinfo[ix][iy]=new CMuTInfo();
		}
	}
	if(SECALC){
		SEinfo=new CSEInfo(this);
	}
}

void CB3D::SetQualifier(string qualifier_set){
	qualifier=qualifier_set;
	string command="mkdir -p model_output/"+run_name+"/"+qualifier;
	system(command.c_str());
	if(sampler!=NULL)
		sampler->nevents=0;
	if(oscarfile!=NULL){
		fclose(oscarfile);
		oscarfile=NULL;
	}
	if(SECALC)
		SEinfo->Zero();
	oscarfilename="model_output/"+run_name+"/"+qualifier+"/oscar.dat";
}

void CB3D::MovePartsToFinalMap(){
	CPartMap::iterator ppos;
	CPart *part;
	CB3DCell *c;
	while(PartMap.size()!=0){
		ppos=PartMap.begin();
		part=ppos->second;
		c=part->cell;
		if(c!=NULL){
			part->DeleteFromMap(&c->partmap);
			part->cell=NULL;
		}
		part->ChangeMap(&FinalPartMap);
	}
}

void CB3D::Reset(){
	int ix,iy,itau,ntau;
	double taucalc;
	CB3DCell *c;
	KillAllParts();
	KillAllActions();
	tau=0.0;
	nactions=0;
	if(MUTCALC){
		itau=0;
		ntau=lrint(TAUCOLLMAX/MUTCALC_DELTAU);
		for(itau=0;itau<ntau;itau++){
			taucalc=(1+itau)*MUTCALC_DELTAU;
			AddAction_MuTCalc(taucalc);
			for(ix=0;ix<2*NXY;ix++){
				for(iy=0;iy<2*NXY;iy++){
					muTinfo[ix][iy]->Tpi=muTinfo[ix][iy]->TK=150.0;
					muTinfo[ix][iy]->muK=muTinfo[ix][iy]->mupi=0.0;
				}
			}
		}
		CMuTInfo::NETEVENTS+=NSAMPLE;
	}
	if(SECALC){
		SEinfo->NETEVENTS+=NSAMPLE;
		SEinfo->ETAOVERS=parameter::getD(parmap,"SEINFO_ETAOVERS",0.3);
		int itau0;
		itau0=lrint(SEinfo->TAU0/SEinfo->DELTAU)-1;
		for(itau=itau0;itau<SEinfo->NTAU;itau++){
			AddAction_SECalc(SEinfo->DELTAU*(itau+1));
		}
	}
}

CB3D::~CB3D(){
	if(oscarfile!=NULL)
		fclose(oscarfile);
}

#endif
