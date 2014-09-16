#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>

const double PI=3.14159265358979323844;
const double HBARC=197.3269602;

using namespace std;

int main(int argc, char *argv[]){
	if (argc != 2) {
		printf("Usage: bmaker NBINS\n");
		exit(-1);
  }
	int NB=atoi(argv[1]),ib;
	double bb2,oldbb2,b;
	double bb2max=216.45; // = sigma_tot/pi, sigmatot=6.8 barns for AuAu at 100 GeV
	
	printf("FOR AUAU at 100A GeV\n");
	oldbb2=0.0;
	for(ib=0;ib<NB;ib++){
		bb2=(double(ib+1)/double(NB))*bb2max;
		b=sqrt(0.5*(bb2+oldbb2));
		oldbb2=bb2;
		printf("%3.0f-%3.0f%%: %g\n",100.0*ib/double(NB),100.0*(ib+1)/double(NB),b);
	}
	
	printf("FOR PbPb at 7? TeV (the total X-section needs to be fixed)\n");
	printf("Right now, I'm using a guess\n");
	bb2max=248.282; // using sigma_tot=7.66 barns (http://dde.web.cern.ch/dde/glauber_lhc.htm)
	oldbb2=0.0;
	for(ib=0;ib<NB;ib++){
		bb2=(double(ib+1)/double(NB))*bb2max;
		b=sqrt(0.5*(bb2+oldbb2));
		oldbb2=bb2;
		printf("%3.0f-%3.0f%%: %g\n",100.0*ib/double(NB),100.0*(ib+1)/double(NB),b);
	}

  return 0;
}


