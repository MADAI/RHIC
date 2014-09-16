#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <cstring>
#include <iostream>
#include <sstream>

const double PI=3.14159265358979323844;
const double HBARC=197.3269602;

using namespace std;

int main(){
	FILE *input=fopen("pdg05.dat","r");
	FILE *output=fopen("resonances_pdg.dat","w");
	fprintf(output,"319\n");
	FILE *decay=fopen("decays_pdg.dat","w");
	string line,pidstring,namestring,massstring,widthstring;
	stringstream ss;
	int i,ichar,istring,pid,n,nread=0,oldpid=0,degen,baryon,strangeness,charge,whoknows1,whoknows2,nchannels,multiplet,idecay;
	int product[20][5],iproduct,ichannel,nproducts[20],ichannelprint,pp;
	bool ppcheck;
	int ires,nres=0,pidlist[330];
	double branching[20],btot;
	double mass,width;
	char wholeline[100];
	pidstring.clear();
	namestring.clear();
	massstring.clear();
	widthstring.clear();

	while(!feof(input)&&nread<1000){
		fgets(wholeline,100,input);
		line=wholeline;
		n=0;
		idecay=0;
		pidstring.assign(&line[0],&line[9]);
		pidstring[10]='\0';
		n+=strlen(pidstring.c_str());
		pid=atoi(pidstring.c_str());
		if(pid!=oldpid){
			btot=0.0;
			ichannel=0;
			//printf("-----------------------------------------\n");
			//printf("%s -> %d\n",pidstring.c_str(),pid);
			namestring.assign(&line[10],&line[27]);
			n+=strlen(namestring.c_str());
			//printf("namestring=%s\n",namestring.c_str());
			for(ichar=0;ichar<27;ichar++){
				line[ichar]=' ';
			}
			ss.clear();
			//printf("line=%s",line.c_str());
			ss << line;
			ss >> mass >> width >> degen >> baryon >> strangeness >> whoknows1 >> whoknows2 >> multiplet >> charge >> nchannels;
			if(width>0.00001) idecay=1;
			//printf("mass=%g, width=%g, degen=%d, baryon=%d, strangeness=%d, multiplet=%d, charge=%d, nchannels=%d\n",mass, width,degen,baryon,strangeness,multiplet,charge,nchannels);
			//printf("namestring=%s\n",namestring.c_str());
			if(nchannels>20){
				printf("nchannels=%d\n",nchannels);
				exit(1);
			}
			pidlist[nres]=pid;
			nres+=1;
			fprintf(output,"%d %g %d %d %d %g %d %g %s\n",pid,1000.0*mass,charge,baryon,strangeness,(degen-1.0)/2.0,idecay,1000.0*width,namestring.c_str());
			if(baryon!=0){
				fprintf(output,"%d %g %d %d %d %g %d %g %s\n",-pid,1000.0*mass,-charge,-baryon,-strangeness,(degen-1.0)/2.0,idecay,1000.0*width,namestring.c_str());
				pidlist[nres]=-pid;
				nres+=1;
			}
			nread+=1;
		}
		else{
			//printf("width=%g, %s -> %d, namestring=%s\n",width,pidstring.c_str(),pid,namestring.c_str());
			idecay=0;
			if(width>0.00001){
				idecay=1;
				//printf("%s -> %d\n",pidstring.c_str(),pid);
				ss << line;
				//printf("XXXXXXXXXX %s\n",line.c_str());
				ss >> pid;
				ss >> nproducts[ichannel];
				if(nproducts[ichannel]<0) nproducts[ichannel]*=(-1);
				ss >> branching[ichannel];
				btot+=branching[ichannel];
				if(nproducts[ichannel]>3) printf("%s: ichannel=%d, pid=%d, nproducts[%d]=%d, branching[%d]=%g\n", namestring.c_str(),ichannel,pid,ichannel,nproducts[ichannel],ichannel,branching[ichannel]);
				for(iproduct=0;iproduct<5;iproduct++)	ss >> product[ichannel][iproduct];
			}
			ichannel+=1;
			if(ichannel==nchannels){
				if(idecay==1 && fabs(btot-1.0)>0.001){
					printf("btot=%g is screwed up\n",btot);
					exit(1);
				}
				//printf("%d %g %d %d %d spin=%g idecay=%d namestring=%s\n",pid,mass,charge,baryon,strangeness,(degen-1.0)/2.0,idecay,namestring.c_str());
				if(idecay==1){
					//
					fprintf(decay,"%d %g %d %d %d %g %d %g %s\n",pid,1000.0*mass,charge,baryon,strangeness,(degen-1.0)/2.0,idecay,1000.0*width,namestring.c_str());
					fprintf(decay,"%d %d ",pid,nchannels);
					for(ichannelprint=0;ichannelprint<nchannels;ichannelprint++){
						fprintf(decay,"%d ",nproducts[ichannelprint]);
						for(iproduct=0;iproduct<nproducts[ichannelprint];iproduct++) fprintf(decay,"%d ",product[ichannelprint][iproduct]);
						fprintf(decay,"%g ",branching[ichannelprint]);
					}
					fprintf(decay,"\n");
					//
					if(baryon!=0){
						fprintf(decay,"%d %g %d %d %d %g %d %s\n",-pid,mass,-charge,-baryon,-strangeness,(degen-1.0)/2.0,idecay,namestring.c_str());
						fprintf(decay,"%d %d ",-pid,nchannels);
						for(ichannelprint=0;ichannelprint<nchannels;ichannelprint++){
							fprintf(decay,"%d ",nproducts[ichannelprint]);
							for(iproduct=0;iproduct<nproducts[ichannelprint];iproduct++){
								pp=product[ichannelprint][iproduct];
								ppcheck=false;
								for(ires=0;ires<nres;ires++){
									if(pidlist[ires]==-pp){
										ppcheck=true;
									}
								}
								if(ppcheck==true) pp=-pp;
								fprintf(decay,"%d ",pp);
							}
							fprintf(decay,"%g ",branching[ichannelprint]);
						}
						fprintf(decay,"\n");
					}
				}
				pidstring.clear();
				namestring.clear();
				massstring.clear();
				widthstring.clear();
			}
		}
		oldpid=pid;
	}
	fclose(decay);
	fclose(output);
	fclose(input);
	printf("successfully finishing, nres=%d\n",nres);
	return 0;
}


