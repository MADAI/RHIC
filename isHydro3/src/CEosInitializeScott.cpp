/*
*  CEos.cpp
*  Created by Joshua Vredevoogd on 3/4/09.
*/

#include "CEos.h"
#include "b3d.h"

//constructor - read from EOS interpolation
void CEosCalculator::InitializeEosScott(parameterMap* pM){
	double epsilonmax=1000000,dx=0.05;
	int NT; // number of points for hadronic EoS below Th
	double Th=parameter::getD(*pM,"EQOFST_SCOTT_TH",165.0);
	NT=lrint(Th-20.0);
	int iT,iX,iT0;
	double epsilonh,Ph,sh,nhadrons,cs2h,cs2,cs20,T0;
	double ff,x0min,x,x0,P,epsilon,oldcs2,oldepsilon,oldP,logToverTh,T;
	vector<double> density;
	double xprime=parameter::getD(*pM,"EQOFST_SCOTT_XPRIME",1.9);
	double x0ratio=parameter::getD(*pM,"EQOFST_SCOTT_X0RATIO",0.45);
	//printf("XPRIME=%g, X0RATIO=%g, Th=%g\n",xprime,x0ratio,Th);
	bExpTail   = false;
	bGaussTail = false;
	bPowerTail = true; 
	
	CResList *reslist=new CResList(pM);
	
	temp[0]=ed[0]=pr[0]=sd[0]=0.0;
	for(iT=0;iT<=NT;iT++){
		T=20.0+double(iT)*(Th-20.0)/double(NT);
		reslist->CalcEoS(T,epsilon,P,nhadrons,cs2,density);
		temp[iT]=T/1000.0;
		ed[iT]=epsilon/1000.0;
		pr[iT]=P/1000.0;
		sd[iT]=(P+epsilon)/T;
	}
	cs2h=cs2;
	epsilonh=epsilon;
	Ph=P;
	
	//
	x0min=-xprime*sqrt(12.0*cs2h);
	x0=x0ratio*fabs(x0min);
	P=Ph;
	epsilon=epsilonh;
	cs2=cs2h;
	x=logToverTh=0.0;
	iT=NT;
	aSize=NT+1;
	while(epsilon<epsilonmax){
		oldepsilon=epsilon;
		oldcs2=cs2;
		oldP=P;
		x+=dx;
		cs2=cs2h+((1.0/3.0)-cs2h)*(x*x+x0*x)/(x*x+x*x0+xprime*xprime);
		epsilon*=exp(dx);
		P+=0.5*(oldcs2+cs2)*(epsilon-oldepsilon);
		logToverTh+=(epsilon-oldepsilon)*(cs2+oldcs2)/(P+oldP+epsilon+oldepsilon);
		T=Th*exp(logToverTh);
		iT+=1;
		//printf("iT=%d, T=%g, epsilon=%g, P=%g, cs2=%g\n",iT,T,epsilon/1000.0,P/1000.0,cs2);
		temp[iT]=T/1000.0;
		sd[iT]=(epsilon+P)/T;
		ed[iT]=epsilon/1000.0;
		pr[iT]=P/1000.0;
		aSize+=1;
		if(aSize>1000){
			printf("make arrays bigger!\n");
			exit(1);
		}
	}
	/*
	for(iT=0;iT<aSize;iT++){
		printf("%4d: %6.4f %10.6f %10.6f %10.6f\n",iT,temp[iT],ed[iT],pr[iT],sd[iT]);
	}
	*/
	
	eAtT = getEGivenT(svSwitchTemp);
	sAtT = getSGivenT(svSwitchTemp);
}
