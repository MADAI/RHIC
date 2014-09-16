#ifndef __SCATTER_CC__
#define __SCATTER_CC__
#include "b3d.h"
using namespace std;

/*
 InelasticScatter: Method to implement inelastic scattering. Right now, it magically
 converts the two incoming particles to outgoing particles and then scatters them
 elastically.
 */

 void CB3D::InelasticScatter(CPart *part1, CPart *part2, CInelasticInfo inelinfo){
 	double mass[3], *p[3], mtot, cthet, sthet, phi, q, u[4], p0tot[4],*p1=part1->p,*p2=part2->p,y1,y2,mt;
 	int ibody, alpha;
 	CPart * dptr;

 	part1->resinfo = inelinfo.resinfo_1;
 	part2->resinfo = inelinfo.resinfo_2;

 	CPart * daughter[2] = {part1, part2 };

 	for(int b = 0; b< 4; b++){
 		p0tot[b] = part1->p[b]+part2->p[b];
 	}
 	p[0] = p0tot;

 	mass[0] = pow(p0tot[0],2);
 	for(alpha = 1; alpha<4; alpha++){
 		mass[0] -= pow(p0tot[alpha],2);
 	}
 	mass[0] = sqrt(mass[0]);

 	mtot=0.0;
 	for(ibody=0;ibody<2;ibody++){
 		mass[ibody+1]=daughter[ibody]->resinfo->mass;
 		mtot+=mass[ibody+1];
 		p[ibody+1]=daughter[ibody]->p;
 	}
 	if(mtot<mass[0]){
		//pick the angles
 		cthet=1.0-2.0*randy->ran();
 		sthet=sqrt(1.0-cthet*cthet);
 		phi=2.0*PI*randy->ran();
 		q=sqrt(Misc::triangle(mass[0],mass[1],mass[2]));
 		p[1][3]=q*cthet;
 		p[1][1]=q*sthet*cos(phi);
 		p[1][2]=q*sthet*sin(phi);
 		p[2][3]=-p[1][3];
 		p[2][2]=-p[1][2];
 		p[2][1]=-p[1][1];
 		p[1][0]=sqrt(mass[1]*mass[1]+p[1][1]*p[1][1]+p[1][2]*p[1][2]+p[1][3]*p[1][3]);
 		p[2][0]=sqrt(mass[2]*mass[2]+p[2][1]*p[2][1]+p[2][2]*p[2][2]+p[2][3]*p[2][3]);
		/* Boost the new particles */
 		for(alpha=0;alpha<4;alpha++) u[alpha]=p0tot[alpha]/mass[0];
 			double pprime[4];
 		for(ibody=0;ibody<2;ibody++){
 			dptr=daughter[ibody];
 			Misc::lorentz(u,p[ibody+1],pprime);
 			for(alpha=0;alpha<4;alpha++) dptr->p[alpha]=pprime[alpha];
 				dptr->SetY();
 			dptr->Setp0();
 		}
 	}

 }

/*
 Scatter: Elastically scatters two CPart objects with s-wave angular distribution
 */

 void CB3D::Scatter(CPart *part1,CPart *part2){
 	double ptot[4],u[4],q[4],qprime[4],roots=0.0,newroots,taucoll,g[4]={1.0,-1.0,-1.0,-1.0};
 	double ctheta,stheta,phi,qmag,vr1,vr2,m1,m2,*p1=part1->p,*p2=part2->p;
 	double y1,mt, mass1,mass2;
 	int nparts,ipart,alpha;

 	m1=part1->GetMass();
 	m2=part2->GetMass();
	//else printf("b=%g, delrperp=%g, deleta=%g\n",sqrt(GetPiBsquared(part1,part2))/PI,
	//sqrt(pow(part1->r[1]-part2->r[1],2)+pow(part1->r[2]-part2->r[2],2)),fabs(part1->eta-part2->eta));

 	for(alpha=0;alpha<4;alpha++){
 		ptot[alpha]=p1[alpha]+p2[alpha];
 		q[alpha]=0.5*(p1[alpha]-p2[alpha]);
 		roots+=g[alpha]*ptot[alpha]*ptot[alpha];
 	}
 	roots=sqrt(roots);
	//printf("BEFORE: roots=%g, ptot=(%g,%g,%g,%g)\n",roots,ptot[0],ptot[1],ptot[2],ptot[3]);
 	for(alpha=0;alpha<4;alpha++) u[alpha]=ptot[alpha]/roots;
 		Misc::BoostToCM(u,q,qprime);
 	qmag=sqrt(qprime[1]*qprime[1]+qprime[2]*qprime[2]+qprime[3]*qprime[3]);

 	ctheta=1.0-2.0*randy->ran();
 	phi=2.0*PI*randy->ran();
 	qprime[3]=qmag*ctheta;
 	stheta=sqrt(1.0-ctheta*ctheta);
 	qprime[1]=qmag*stheta*cos(phi);
 	qprime[2]=qmag*stheta*sin(phi);
	//printf("qprime[0]=%g =? %g\n",qprime[0],0.5*sqrt(m1*m1+qmag*qmag)-0.5*sqrt(m2*m2+qmag*qmag));

 	qprime[0]=sqrt(m1*m1+qmag*qmag);
 	Misc::Boost(u,qprime,p1);
 	part1->SetY();
 	for(alpha=1;alpha<4;alpha++) qprime[alpha]=-qprime[alpha];
 		qprime[0]=sqrt(m2*m2+qmag*qmag);
 	Misc::Boost(u,qprime,p2);
 	part2->SetY();

 	mass1 = part1->GetMass();
 	mass2 = part2->GetMass();
 	if(mass1 < part1->resinfo->minmass){
 		cout << "Error in Scatter: Created a particle with wrong mass." <<endl;
 		cout << "Mass: " << mass1 << " minmass: " << part1->resinfo->minmass;
 		part1->Print();
 		part1->resinfo->Print();
 		exit(-1);
 	}
 	if(mass2 < part2->resinfo->minmass){
 		cout << "Error in Scatter: Created a particle with wrong mass." <<endl;
 		part2->Print();
 		exit(-1);
 	}
 }

/*
 If resonance can decay, A->B+C, then if B+C collide, this will merge B+C->A
 */
  bool CB3D::Merge(CPart *part1,CPart *part2,CResInfo *resinfo){
 	int alpha;
 	double ptot[4],minv,mt,y1,y2;
 	CPart *pswitch;
 	double *p1=part1->p,*p2=part2->p;
 	bool success=false;

 	for(alpha=0;alpha<4;alpha++)
 		ptot[alpha]=part1->p[alpha]+part2->p[alpha];
 	minv=sqrt(ptot[0]*ptot[0]-ptot[1]*ptot[1]-ptot[2]*ptot[2]-ptot[3]*ptot[3]);
 	if(minv>resinfo->minmass){
 		for(alpha=0;alpha<4;alpha++)
 			part1->p[alpha]=part1->p[alpha]+part2->p[alpha];
 		part1->resinfo=resinfo;
 		part1->SetY();
 		part2->KillActions();
 		part2->Kill();
 		success=true;
 	}
 	//else{
 		//printf("merge failed, minv=%g, minmass=%g\n",minv,resinfo->minmass);
 	//}
 	return success;
 }

#endif