#ifndef __ANALYZE_REALITYDIFF_CC__
#define __ANALYZE_REALITYDIFF_CC__
#include "b3d.h"
//
void CB3D::CalcRealityDiff(){
	const int NPT=32;
	double weight,et,pt,Et,NetEt,NetEtsquared,vr,rperp;
	int ipt,nevents,ievent,ipart,nparts;
	double spectradiff[NPT]={0.0};
	long long int netfakeparts;
	CPartMap::iterator ppos;
	CPart *part;
	FILE *fptr;
	parameter::set(parmap,string("B3D_RESONANCES_DECAYS_FILE"),string("progdata/madai/resinfo/decays_none.dat"));
	parameter::set(parmap,string("B3D_RESONANCES_INFO_FILE"),string("progdata/madai/resinfo/resonances_nodecays.dat"));
	//parameter::set(parmap,string("B3D_RESONANCES_DECAYS_FILE"),string("progdata/madai/resinfo/decays_pdg_weak.dat"));
	//parameter::set(parmap,string("B3D_RESONANCES_INFO_FILE"),string("progdata/madai/resinfo/resonances_pdg_weak.dat"));
	reslist=new CResList();
	reslist->ReadResInfo();
	COLLISIONS=false;
	nevents=parameter::getI(parmap,"B3D_NEVENTSMAX",4000);
	netfakeparts=0;
	NetEtsquared=NetEt=0.0;
	fptr=fopen("vrtau.dat","w");
	for(ievent=1;ievent<=nevents;ievent++){
		KillAllParts();
		nparts=ReadOSCAR(ievent);
		if(nparts>0){
			PerformAllActions(); // Decays unstable particles
			Et=0.0;
			for(ppos=FinalPartMap.begin();ppos!=FinalPartMap.end();ppos++){
				part=ppos->second;
				if(!part->reality){
					if(part->tau_lastint<0.8){
						part->Print();
					}
					part->Propagate(part->tau_lastint);
					netfakeparts+=1;
					weight=part->weight;
					pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
					et=sqrt(pt*pt+part->msquared)*weight;
					rperp=sqrt(part->r[1]*part->r[1]+part->r[2]*part->r[2]);
					vr=(part->p[1]*part->r[1]+part->p[2]*part->r[2])/(fabs(et)*rperp);
					//if(fabs(part->resinfo->code)==211 || part->resinfo->code==111){
						Et+=et;
						//}
					if(weight==-1)
						fprintf(fptr,"%8.4f %5d %7.3f %7.5f %7.3f\n",part->tau0,part->resinfo->code,pt,vr,rperp);
					//if(fabs(part->resinfo->code)==211 || part->resinfo->code==111){
						ipt=lrint(floor(pt/50.0));
						if(ipt<NPT){
							spectradiff[ipt]+=weight;
						}
						//}
				}
			}
			NetEtsquared+=Et*Et;
			NetEt+=Et;
		}
	}
	fclose(fptr);
	NetEt=NetEt/double(nevents);
	NetEtsquared=NetEtsquared/double(nevents);
	printf("<Nfake>=%g, fake Et/event=%g += %g\n",double(netfakeparts)/double(nevents),NetEt,sqrt((NetEtsquared-NetEt*NetEt)/double(nevents)));
	
	string filename="analysis/reality.dat";
	fptr=fopen(filename.c_str(),"a");
	fprintf(fptr,"%2d %g %g %g\n",NSCATT_MAX,double(netfakeparts)/double(nevents),NetEt,sqrt((NetEtsquared-NetEt*NetEt)/double(nevents)));
	fclose(fptr);
	filename="analysis/"+run_name+"/spectradiff_pion.dat";
	fptr=fopen(filename.c_str(),"w");
	for(ipt=0;ipt<NPT;ipt++){
		fprintf(fptr,"%6.1f %g\n",(0.5+ipt)*50.0,spectradiff[ipt]);
	}
	fclose(fptr);
	delete reslist;
}
#endif
