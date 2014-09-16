#ifndef __INCLUDE_SOURCECALC_HDF5__
#define __INCLUDE_SOURCECALC_HDF5__
#ifdef __B3D_USE_HDF5__
#include "sourcecalc.h"

using namespace std;

CSourceCalc_HDF5::CSourceCalc_HDF5(){
	InitSPars();	
	parameter::PrintPars(spars);
	randy=new CRandom(1234);
}

CSourceCalc_HDF5::CSourceCalc_HDF5(string sparsfilename){
	InitSPars();
	parameter::ReadParsFromFile(spars,sparsfilename);
	parameter::PrintPars(spars);
	randy=new CRandom(1234);
}

void CSourceCalc_HDF5::InitSPars(){
	// DEFAULT VALUES
	parameter::set(spars,"PT",600.0);
	parameter::set(spars,"DELPT",20.0);
	parameter::set(spars,"PHIMIN_DEG",0.0);
	parameter::set(spars,"PHIMAX_DEG",360.0);
	parameter::set(spars,"YMIN",-1.0);
	parameter::set(spars,"YMAX",1.0);
	parameter::set(spars,"NPARTSMAX",20000);
	parameter::set(spars,"HDF5filename",string("HDF5_UNDEFINED_FILENAME"));
	parameter::set(spars,"NEVENTSMAX",10000);
	parameter::set(spars,"ETA_GAUSS",1.2);
}

void CSourceCalc_HDF5::SetSPars(double PT_set,double DELPT_set,double PHIMIN_DEG_set,double PHIMAX_DEG_set,double YMIN_set,double YMAX_set){
	parameter::set(spars,"PT",PT_set);
	parameter::set(spars,"DELPT",DELPT_set);
	parameter::set(spars,"PHIMIN_DEG",PHIMIN_DEG_set);
	parameter::set(spars,"PHIMAX_DEG",PHIMAX_DEG_set);
	parameter::set(spars,"YMIN",YMIN_set);
	parameter::set(spars,"YMAX",YMAX_set);
}

void CSourceCalc_HDF5::SetIDs(int *idlista,int nida,int *idlistb,int nidb){
	idlist_a=idlista;
	nid_a=nida;
	nid_b=nidb;
	idlist_b=idlistb;
}

void CSourceCalc_HDF5::CalcS(CMCPRList *&lista,CMCPRList *&listb){
	double **ra,**rb,**pa,**pb;
	int ia,ib,na=0,nb=0;
	bool AEQUALB=parameter::getB(spars,"AEQUALB",true);
	int NPARTSMAX=parameter::getI(spars,"NPARTSMAX",50000);
	if(lista!=NULL){
		delete lista;
		lista=NULL;
	}
	if(listb!=NULL){
		delete listb;
		listb=NULL;
	}

	ra=new double *[NPARTSMAX];
	pa=new double *[NPARTSMAX];
	for(ia=0;ia<NPARTSMAX;ia++){
		ra[ia]=new double[4];
		pa[ia]=new double[4];
	}
	if(AEQUALB){
		rb=ra;
		pb=pa;
	}
	else{
		rb=new double *[NPARTSMAX];
		pb=new double *[NPARTSMAX];
		for(ib=0;ib<NPARTSMAX;ib++){
			rb[ib]=new double[4];
			pb[ib]=new double[4];
		}
	}
	ReadPR(pa,ra,na,pb,rb,nb);

	lista=new CMCPRList(na);
	if(!AEQUALB){
		listb=new CMCPRList(nb);
	}

	for(int i=0;i<na;i++) lista->SetPR(i,pa[i],ra[i]);
	if(lista!=listb){
		for(int i=0;i<nb;i++) listb->SetPR(i,pb[i],rb[i]);
	}

	for(ia=0;ia<NPARTSMAX;ia++){
		delete [] ra[ia];
		delete [] pa[ia];
	}
	delete [] ra;
	delete [] pa;
	if(!AEQUALB){
		for(ib=0;ib<NPARTSMAX;ib++){
			delete [] rb[ib];
			delete [] pb[ib];
		}
		delete [] rb;
		delete [] pb;
	}
	printf("______ FINISHED CREATING MCLISTS na=%d, nb=%d___________\n",na,nb);

}

void CSourceCalc_HDF5::ReadPR(double **pa,double **ra,int &na,double **pb,double **rb,int &nb){
	double r[4],p[4];
	double MA=parameter::getD(spars,"MA",139.57);
	double MB=parameter::getD(spars,"MB",139.57);
	double mass,mt,tau,eta;
	int ID,idummy,alpha,iipart;
	int ipart,nparts,ievent,nevents;
	int NEVENTSMAX=parameter::getI(spars,"NEVENTSMAX",100);
	bool AEQUALB=parameter::getB(spars,"AEQUALB",false);
	CPartH5 *partH5;
	CompType *ptype;
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
	
	string h5_infilename=parameter::getS(spars,"HDF5filename","HDF5filename_UNDEFINED");
	printf("CAnalyze::CalcSpectra, opening %s, ",h5_infilename.c_str());
	H5File *h5file=new H5File(h5_infilename,H5F_ACC_RDONLY);
	char dummy[160];
	long long int nsigmay=0;
	//double y,eta,sigmay=0.0,sigmaz=0.0,z;
	//double taucompare=parameter::getD(spars,"TAUCOMPARE",12.0);
	//taucompare=12.0;
	na=nb=0;
	nevents=int(h5file->getNumObjs());
	printf("nevents=%d\n",nevents);
	char eventno[10];
	H5D_space_status_t status;
	hsize_t dim[1];
	if(nevents>NEVENTSMAX) nevents=NEVENTSMAX;
	for(ievent=1;ievent<=nevents;ievent++){
		randy->reset(ievent);
		sprintf(eventno,"event%d",ievent);
		DataSet *dataset = new DataSet (h5file->openDataSet(eventno));
		//dataset->getSpaceStatus(status);
		//hsize_t datasetsize=dataset->getStorageSize();
		//printf("ievent=%d, status=%d, size=%d\n",ievent,int(status),int(datasetsize));
		DataSpace filespace = dataset->getSpace();
		int rank=filespace.getSimpleExtentDims(dim);
		nparts=dim[0];
		partH5=new CPartH5[nparts];
		dataset->read(partH5,*ptype);
		for(ipart=0;ipart<nparts;ipart++){
			ID=partH5[ipart].ID;
			p[1]=partH5[ipart].px;
			p[2]=partH5[ipart].py;
			r[1]=partH5[ipart].x;
			r[2]=partH5[ipart].y;
			mass=partH5[ipart].mass;
			mt=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]);
			p[3]=mt*sinh(partH5[ipart].rapidity);
			p[0]=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
			tau=partH5[ipart].tau;
			eta=partH5[ipart].eta;
			r[0]=tau*cosh(eta);
			r[3]=tau*sinh(eta);
			//if(ipart<10) printf("ipart=%d, p=(%g,%g,%g,%g) r=(%g,%g,%g,%g)\n",ipart,p[0],p[1],p[2],p[3],r[0],r[1],r[2],r[3]);
			if(IDMatch(ID,idlist_a,nid_a)) Check(p,r,MA,pa,ra,na);
			if(!AEQUALB) if(IDMatch(ID,idlist_b,nid_b)) Check(p,r,MB,pb,rb,nb);
		}
		dataset->close();
		delete dataset;
		delete [] partH5;
	}	
	h5file->close();
	delete h5file;
	delete ptype;
	if(AEQUALB) nb=na;
}

bool CSourceCalc_HDF5::Check(double *p,double *r,double m,double **pa,double **ra,int &n){
	double PT=parameter::getD(spars,"PT",600.0);
	double DELPT=parameter::getD(spars,"DELPT",20.0);
	double YMIN=parameter::getD(spars,"YMIN",-1.0);
	double YMAX=parameter::getD(spars,"YMAX",1.0);
	double PHIMIN=2.0*PI*parameter::getD(spars,"PHIMIN_DEG",0.0)/360.0;
	double PHIMAX=2.0*PI*parameter::getD(spars,"PHIMAX_DEG",0.0)/360.0;
	double ETA_GAUSS=parameter::getD(spars,"ETA_GAUSS",1.2);
	double MA=parameter::getD(spars,"MA",139.57);
	double MB=parameter::getD(spars,"MB",139.57);
	double pt,phi,eta;
	double gammav=PT/(MA+MB);
	double gamma=sqrt(1.0+gammav*gammav);
	double pttarget=PT*m/(MA+MB);
	double rout,rlong,rside,sinhy,coshy,tau,vperp,y;
	const double TAUCOMPARE=12.0;
	int NPARTSMAX=parameter::getI(spars,"NPARTSMAX",50000);
	bool XREFLECTIONSYMMETRY=parameter::getB(spars,"XREFLECTIONSYMMETRY",false);
	bool YREFLECTIONSYMMETRY=parameter::getB(spars,"YREFLECTIONSYMMETRY",false);
	bool success=false;
	
	pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	if(fabs(pt-pttarget)<DELPT){
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
		if( ((PHIMIN<PHIMAX)&&(phi>PHIMIN && phi<PHIMAX)) || ((PHIMIN>PHIMAX)&&(phi>PHIMIN||phi<PHIMAX)) ){
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
						printf("TOO MANY PARTICLES FIT CRITERIA, increase parameter NPARTSMAX=%d\n",NPARTSMAX);
						exit(1);
					}
					//printf("%g %g %g %g\n",tau,rout,rside,rlong);
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

bool CSourceCalc_HDF5::IDMatch(int ident,int *idlist,int nid){
	int i=0;
	bool answer=false;
	while(answer==false && i<nid){
		if(ident==idlist[i]) answer=true;
		i+=1;
	}
	return answer;
}

#endif
#endif
