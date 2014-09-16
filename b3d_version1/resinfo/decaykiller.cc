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
	FILE *input=fopen("resonances_pdg.dat","r");
	FILE *output=fopen("resonances_nodecays.dat","w");
	fprintf(output,"319\n");
	int baryon,strangeness,charge,idecay,pid;
	double mass,width,spin;
	char name[100];


	while(!feof(input)){
		fscanf(input,"%d %lf %d %d %d %lf %d %lf",&pid,&mass,&charge,&baryon,&strangeness,&spin,&idecay,&width);
		idecay=0;
		fgets(name,100,input);
		fprintf(output,"%d %g %d %d %d %g %d %g %s\n",pid,mass,charge,baryon,strangeness,spin,idecay,width,name);
		fclose(output);
		fclose(input);
		return 0;
	}
}


