#ifndef __INCLUDE_SOURCECALC_HDF5__
#define __INCLUDE_SOURCECALC_HDF5__
#include "sourcecalc.h"

using namespace std;

CSourceCalc_HDF5_MultiBin::CSourceCalc_HDF5_MultiBin(){
	InitSPars();	
	randy=new CRandom(1234);
#ifdef __B3D_USE_HDF5__
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
}

CSourceCalc_HDF5_MultiBin::CSourceCalc_HDF5_MultiBin(string sparsfilename){
	InitSPars();
	parameter::ReadParsFromFile(spars,sparsfilename);
	NPHIBINS=parameter::getI(spars,"NPHIBINS",6);
	NPTBINS=parameter::getI(spars,"NPTBINS",16);
	DELPT=parameter::getD(spars,"DELPT",25.0);
	PTMIN=parameter::getD(spars,"PTMIN",150.0);
	PTMAX=PTMIN+NPTBINS*DELPT;
	DELPHI=90.0/double(NPHIBINS);
	parameter::PrintPars(spars);
	randy=new CRandom(1234);
}

void CSourceCalc_HDF5_MultiBin::InitSPars(){
		// DEFAULT VALUES
	parameter::set(spars,"DELPT",25.0);
	parameter::set(spars,"NPTBINS",16);
	parameter::set(spars,"NPHIBINS",6);
	parameter::set(spars,"PTMIN",150.0);
	parameter::set(spars,"PHIMIN_DEG",0.0);
	parameter::set(spars,"PHIMAX_DEG",360.0);
	parameter::set(spars,"YMIN",-1.0);
	parameter::set(spars,"YMAX",1.0);
	parameter::set(spars,"NPARTSMAX",10000);
	parameter::set(spars,"HDF5filename",string("HDF5_UNDEFINED_FILENAME"));
	parameter::set(spars,"NEVENTSMAX",10000);
	parameter::set(spars,"ETA_GAUSS",1.2);
}

void CSourceCalc_HDF5_MultiBin::SetSPars(double PT_set,double DELPT_set,double PHIMIN_DEG_set,double PHIMAX_DEG_set,double YMIN_set,double YMAX_set){
	parameter::set(spars,"PT",PT_set);
	parameter::set(spars,"DELPT",DELPT_set);
	parameter::set(spars,"PHIMIN_DEG",PHIMIN_DEG_set);
	parameter::set(spars,"PHIMAX_DEG",PHIMAX_DEG_set);
	parameter::set(spars,"YMIN",YMIN_set);
	parameter::set(spars,"YMAX",YMAX_set);
}

void CSourceCalc_HDF5_MultiBin::SetIDs(int *idlista,int nida,int *idlistb,int nidb){
	idlist_a=idlista;
	nid_a=nida;
	nid_b=nidb;
	idlist_b=idlistb;
}

void CSourceCalc_HDF5_MultiBin::CalcS(CMCPRList ***&lista,CMCPRList ***&listb){
	double ****ra,****rb,****pa,****pb;
	int ia,ib,**na,**nb,ipt,iphi;
	bool AEQUALB=parameter::getB(spars,"AEQUALB",true);
	int NPARTSMAX=parameter::getI(spars,"NPARTSMAX",10000);
	
	if(lista!=NULL){
		for(ipt=0;ipt<NPTBINS;ipt++){
			for(iphi=0;iphi<NPHIBINS;iphi++){
				if(lista[ipt][iphi]!=NULL){
					delete lista[ipt][iphi];
					lista[ipt][iphi]=NULL;
				}
				delete [] lista[ipt];
			}
			delete [] lista;
		}
	}
	lista=NULL;
	if(!AEQUALB && listb!=NULL){
		for(ipt=0;ipt<NPTBINS;ipt++){
			for(iphi=0;iphi<NPHIBINS;iphi++){
				if(listb[ipt][iphi]!=NULL){
					delete listb[ipt][iphi];
					listb[ipt][iphi]=NULL;
				}
				delete [] listb[ipt];
			}
			delete [] listb;
		}
		listb=NULL;
	}
	
	lista=new CMCPRList**[NPTBINS];
	for(ipt=0;ipt<NPTBINS;ipt++){
		lista[ipt]=new CMCPRList*[NPHIBINS];
		for(iphi=0;iphi<NPHIBINS;iphi++) lista[ipt][iphi]=NULL;
	}
	if(!AEQUALB){
		listb=new CMCPRList**[NPTBINS];
		for(ipt=0;ipt<NPTBINS;ipt++){
			listb[ipt]=new CMCPRList*[NPHIBINS];
			for(iphi=0;iphi<NPHIBINS;iphi++) listb[ipt][iphi]=NULL;
		}
	}
	
	ra=new double ***[NPTBINS];
	pa=new double ***[NPTBINS];
	na=new int *[NPTBINS];
	for(ipt=0;ipt<NPTBINS;ipt++){
		ra[ipt]=new double **[NPHIBINS];
		pa[ipt]=new double **[NPHIBINS];
		na[ipt]=new int[NPHIBINS];
		for(iphi=0;iphi<NPHIBINS;iphi++){
			ra[ipt][iphi]=new double *[NPARTSMAX];
			pa[ipt][iphi]=new double*[NPARTSMAX];
			for(ia=0;ia<NPARTSMAX;ia++){
				ra[ipt][iphi][ia]=new double[4];
				pa[ipt][iphi][ia]=new double[4];
			}
		}
	}
	
	
	if(AEQUALB){
		rb=ra;
		pb=pa;
		nb=na;
	}
	else{
		rb=new double ***[NPTBINS];
		pb=new double ***[NPTBINS];
		nb=new int *[NPTBINS];
		for(ipt=0;ipt<NPTBINS;ipt++){
			rb[ipt]=new double **[NPHIBINS];
			pb[ipt]=new double **[NPHIBINS];
			nb[ipt]=new int[NPHIBINS];
			for(iphi=0;iphi<NPHIBINS;iphi++){
				rb[ipt][iphi]=new double *[NPARTSMAX];
				pb[ipt][iphi]=new double *[NPARTSMAX];
				for(ia=0;ia<NPARTSMAX;ia++){
					rb[ipt][iphi][ia]=new double[4];
					pb[ipt][iphi][ia]=new double[4];
				}
			}
		}
	}
	
	ReadPR(pa,ra,na,pb,rb,nb);

	for(ipt=0;ipt<NPTBINS;ipt++){
		for(iphi=0;iphi<NPHIBINS;iphi++){	
			lista[ipt][iphi]=new CMCPRList(na[ipt][iphi]);
			if(!AEQUALB){
				listb[ipt][iphi]=new CMCPRList(nb[iphi][iphi]);
			}
			for(int ia=0;ia<na[ipt][iphi];ia++)
				lista[ipt][iphi]->SetPR(ia,pa[ipt][iphi][ia],ra[ipt][iphi][ia]);
			if(lista!=listb){
				for(int ib=0;ib<nb[ipt][iphi];ib++)
					listb[ipt][iphi]->SetPR(ib,pb[ipt][iphi][ib],rb[ipt][iphi][ib]);
			}
		}
	}
	
	for(ipt=0;ipt<NPTBINS;ipt++){
		for(iphi=0;iphi<NPHIBINS;iphi++){
			for(ia=0;ia<NPARTSMAX;ia++){
				delete [] ra[ipt][iphi][ia];
				delete [] pa[ipt][iphi][ia];
			}
			delete [] ra[ipt][iphi];
			delete [] pa[ipt][iphi];
		}
		delete [] ra[ipt];
		delete [] pa[ipt];
		delete [] na[ipt];
	}
	delete [] ra;
	delete [] pa;
	delete [] na;
	if(!AEQUALB){
		for(ipt=0;ipt<NPTBINS;ipt++){
			for(iphi=0;iphi<NPHIBINS;iphi++){
				for(ia=0;ia<NPARTSMAX;ia++){
					delete [] rb[ipt][iphi][ia];
					delete [] pb[ipt][iphi][ia];
				}
				delete [] rb[ipt][iphi];
				delete [] pb[ipt][iphi];
			}
			delete [] rb[ipt];
			delete [] pb[ipt];
			delete [] nb[ipt];
		}
	}
}

void CSourceCalc_HDF5_MultiBin::ReadPR(double ****pa,double ****ra,int **&na,double ****pb,double ****rb,int **&nb){
	int ipt,iphi,i;
	double r[4],p[4],phi,pt;
	double MA=parameter::getD(spars,"MA",139.57);
	double MB=parameter::getD(spars,"MB",139.57);
	double mass,mt,tau,eta,rdummy1,rdummy2;
	int ID,idummy,alpha,iipart;
	int ipart,nparts,ievent,nevents;
	int NEVENTSMAX=parameter::getI(spars,"NEVENTSMAX",100);
	bool AEQUALB=parameter::getB(spars,"AEQUALB",false);
	int NPARTSMAX=parameter::getI(spars,"NPARTSMAX",10000);
#ifdef __B3D_USE_HDF5__
	CPartH5 *partH5;
	string h5_infilename=parameter::getS(spars,"HDF5filename","HDF5filename_UNDEFINED");
	printf("CSourceCalc_HDF5_MultiBin::ReadPR, opening %s\n",h5_infilename.c_str());
	H5File *h5file=new H5File(h5_infilename,H5F_ACC_RDONLY);
#else
	string oscarfilename=parameter::getS(spars,"OSCARfilename","OSCARfilename_UNDEFINED");
	FILE *oscarfile=fopen(oscarfilename.c_str(),"r");
#endif
	char dummy[160];
	long long int nsigmay=0;
		//double y,eta,sigmay=0.0,sigmaz=0.0,z;
		//double taucompare=parameter::getD(spars,"TAUCOMPARE",12.0);
		//taucompare=12.0;
	for(ipt=0;ipt<NPTBINS;ipt++){
		for(iphi=0;iphi<NPHIBINS;iphi++){
			na[ipt][iphi]=0;
			if(!AEQUALB) nb[ipt][iphi]=0;
		}
	}
#ifdef __B3D_USE_HDF5__
	nevents=int(h5file->getNumObjs());
	char eventno[20];
	H5D_space_status_t status;
	hsize_t dim[1];
	if(nevents>NEVENTSMAX) nevents=NEVENTSMAX;
	partH5=new CPartH5[NPARTSMAX];
#else
	nevents=NEVENTSMAX;
#endif
	printf("nevents=%d\n",nevents);
	for(ievent=1;ievent<=nevents;ievent++){
		randy->reset(ievent);
#ifdef __B3D_USE_HDF5__
		sprintf(eventno,"event%d",ievent);
		DataSet *dataset = new DataSet (h5file->openDataSet(eventno));
			//DataSet dataset(h5file->openDataSet(eventno));
			//dataset->getSpaceStatus(status);
			//hsize_t datasetsize=dataset->getStorageSize();
			//printf("ievent=%d, status=%d, size=%d\n",ievent,int(status),int(datasetsize));
		DataSpace filespace = dataset->getSpace();
		int rank=filespace.getSimpleExtentDims(dim);
		nparts=dim[0];
		dataset->read(partH5,*ptype);
#else
		if(ievent==0){
			for(i=0;i<3;i++)
				fgets(dummy,160,oscarfile);
		}
		fscanf(oscarfile,"%d %d %lf %lf",&idummy,&nparts,&rdummy1,&rdummy2);
		fgets(dummy,160,oscarfile);
		//printf("ievent=%d, idummy=%d, npartmax=%d\n",ievent,idummy,npartmax);
#endif
		for(ipart=0;ipart<nparts;ipart++){
#ifdef __B3D_USE_HDF5__
			ID=partH5[ipart].ID;
			p[1]=partH5[ipart].px;
			p[2]=partH5[ipart].py;
#else
			fscanf(oscarfile,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
				&idummy,&ID,&p[1],&p[2],&p[3],&p[0],&mass,&r[1],&r[2],&r[3],&r[0]);
			//for(alpha=0;alpha<4;alpha++) p[alpha]*=1000.0; mass*=1000;
			if(idummy!=ipart){
				printf("In CSourceCalc_HDF5_MultiBin::ReadPR, idummy=%d != ipart=%d\n",idummy,ipart);
				exit(1);
			}
#endif
			pt=sqrt(p[1]*p[1]+p[2]*p[2]);
			if(pt>PTMIN && pt<PTMAX){
#ifdef __B3D_USE_HDF5__
				r[1]=partH5[ipart].x;
				r[2]=partH5[ipart].y;
				tau=partH5[ipart].tau;
				eta=partH5[ipart].eta;
				r[0]=tau*cosh(eta);
				r[3]=tau*sinh(eta);
				mass=partH5[ipart].mass;
				mt=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]);
				p[3]=mt*sinh(partH5[ipart].rapidity);
#endif
				mt=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]);
				p[0]=sqrt(mt*mt+p[3]*p[3]);				
				phi=fabs(atan(p[2]/p[1])*180.0/PI);
				iphi=int(lrint(floor(phi/DELPHI)));
				ipt=int(lrint(floor((pt-PTMIN)/DELPT)));
				if(ipt>=NPTBINS || ipt<0){
					printf("ipt=%d is out of range\n",ipt);
					exit(1);
				}
				if(iphi<0 || iphi>=NPHIBINS){
					printf("iphi=%d is out of range\n",iphi);
					exit(1);
				}
				if(IDMatch(ID,idlist_a,nid_a) && na[ipt][iphi]<NPARTSMAX)
					Check(p,r,MA,pa[ipt][iphi],ra[ipt][iphi],na[ipt][iphi]);
				if(!AEQUALB){
					if(IDMatch(ID,idlist_b,nid_b) && nb[ipt][iphi]<NPARTSMAX){
						Check(p,r,MB,pb[ipt][iphi],rb[ipt][iphi],nb[ipt][iphi]);
					}
				}
			}
		}
#ifdef __B3D_USE_HDF5__
		dataset->close();
		delete dataset;
#endif
	}
#ifdef __B3D_USE_HDF5__
	delete [] partH5;
	h5file->close();
	delete h5file;
#else
	fclose(oscarfile);
#endif
	if(AEQUALB) nb=na;
}

bool CSourceCalc_HDF5_MultiBin::Check(double *p,double *r,double m,double **pa,double **ra,int &n){
	double YMIN=parameter::getD(spars,"YMIN",-1.0);
	double YMAX=parameter::getD(spars,"YMAX",1.0);
	double PHIMIN=2.0*PI*parameter::getD(spars,"PHIMIN_DEG",0.0)/360.0;
	double PHIMAX=2.0*PI*parameter::getD(spars,"PHIMAX_DEG",0.0)/360.0;
	double ETA_GAUSS=parameter::getD(spars,"ETA_GAUSS",1.2);
	double MA=parameter::getD(spars,"MA",139.57);
	double MB=parameter::getD(spars,"MB",139.57);
	double phi,eta;
	double rout,rlong,rside,sinhy,coshy,tau,vperp,y;
	const double TAUCOMPARE=12.0;
	int NPARTSMAX=parameter::getI(spars,"NPARTSMAX",20000);
	bool XREFLECTIONSYMMETRY=parameter::getB(spars,"XREFLECTIONSYMMETRY",false);
	bool YREFLECTIONSYMMETRY=parameter::getB(spars,"YREFLECTIONSYMMETRY",false);
	bool success=false;
	double pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	double gammav=pt/MA;
	double gamma=sqrt(1.0+gammav*gammav);
	if(n<NPARTSMAX){
	
		if(p[1]!=p[1] || p[2]!=p[2] || p[3]!=p[3]){
			printf("bad particle has nan, p=(%g,%g,%g)\n",p[1],p[2],p[3]);
			return false;
		}
		if(XREFLECTIONSYMMETRY){
			if(p[1]<0.0){
				p[1]=-p[1];
				r[1]=-r[1];
			}
		}
		if(YREFLECTIONSYMMETRY){
			if(p[2]<0.0){
				p[2]=-p[2];
				r[2]=-r[2];
			}
		}
		phi=atan2(p[1],p[2]);
		if(phi>PHIMIN && phi<PHIMAX){
			y=atanh(p[3]/p[0]);
			if(y>YMIN && y<YMAX){
				p[0]=sqrt(pt*pt+p[3]*p[3]+m*m);
				rout=(p[1]*r[1]+p[2]*r[2])/pt;
				rside=(p[1]*r[2]-p[2]*r[1])/pt;
				sinhy=sinh(y);
				coshy=cosh(y);
				rlong=coshy*r[3]-sinhy*r[0];
				tau=coshy*r[0]-sinhy*r[3];
				eta=asinh(rlong/tau);
				if(randy->ran()<exp(-0.5*eta*eta/(ETA_GAUSS*ETA_GAUSS))){
					if(n==NPARTSMAX){
						printf("TOO MANY PARTICLES FIT CRITERIA, increase parameter NPARTSMAX=%d if you want more\n",NPARTSMAX);
							//exit(1);
					}
					vperp=pt/sqrt(m*m+pt*pt);
					rout=rout-vperp*(tau-TAUCOMPARE);
					tau=TAUCOMPARE;
					rout=gamma*rout-gammav*tau;
					ra[n][0]=0.0;
					ra[n][1]=rout;
					ra[n][2]=rside;
					ra[n][3]=rlong;
					pa[n][0]=p[0];
					pa[n][1]=p[1];
					pa[n][2]=p[2];
					pa[n][3]=p[3];
					n+=1;
					success=true;
				}
			}
		}
	}
	return success;
}

bool CSourceCalc_HDF5_MultiBin::IDMatch(int ident,int *idlist,int nid){
	int i=0;
	bool answer=false;
	while(answer==false && i<nid){
		if(ident==idlist[i]) answer=true;
		i+=1;
	}
	return answer;
}

#endif
