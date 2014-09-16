#ifndef __QUALIFIER_CC__
#define __QUALIFIER_CC__

#include "qualifier.h"
using namespace std;

void CQualifiers::Read(string qfilename){
	char cread[120],dummy[120];
	int ipar;
	string sread;
	nqualifiers=-1;
	FILE *fptr=fopen(qfilename.c_str(),"r");
	while(!feof(fptr)){
		fscanf(fptr,"%s",cread); sread=cread;
		if(sread=="qualifier"){
			nqualifiers+=1;
			npars[nqualifiers]=0;
			fscanf(fptr,"%s",cread); qualifier[nqualifiers]=cread;
		}
		else if(sread=="int" || sread=="double" || sread=="bool"){
			type[nqualifiers][npars[nqualifiers]]=sread;
			fscanf(fptr,"%s",cread); parname[nqualifiers][npars[nqualifiers]]=cread;
			fscanf(fptr,"%s",cread); value[nqualifiers][npars[nqualifiers]]=cread;
			npars[nqualifiers]+=1;
			fgets(dummy,120,fptr);
		}
		else{
			if(sread[0]=='#'){
				fgets(dummy,120,fptr);
			}
			else{
				strcpy(cread,sread.c_str());
				type[nqualifiers][npars[nqualifiers]]="unknown";
				parname[nqualifiers][npars[nqualifiers]]=cread;
				fscanf(fptr,"%s",cread); value[nqualifiers][npars[nqualifiers]]=cread;
				npars[nqualifiers]+=1;
				fgets(dummy,120,fptr);
			}
		}
	}
	nqualifiers+=1;
	fclose(fptr);
}

void CQualifiers::SetPars(parameterMap *pmap,int iqual){
	for(int ipar=0;ipar<npars[iqual];ipar++){
		parameter::set(*pmap,parname[iqual][ipar],value[iqual][ipar]);
	}
}

#endif