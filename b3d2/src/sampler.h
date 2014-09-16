#ifndef __SAMPLER_H__
#define __SAMPLER_H__
#include "b3d.h"
//#define __SAMPLER_WRITE_XY__

using namespace std;
class Csampler;
class Cvertex2D;

class CvolumeElement2D{
public:
	double T,muB,muE,muS,XUD,XS,Xscale;
	double nhadrons,nbaryonsS0,nbaryonsS1,nbaryonsS2,nbaryonsS3;
	double nmesonsS0,nmesonsS1,nmesonsS2;
	vector<double> *density; // densities of various species
	FourVector Omega; // volume elements (four-vectors)
	double Omegamax;
	Cvertex2D *vertex[3]; //vertices of triangle
	double pitilde[4][4]; // shear tensor
	void Initialize();
	void CalcEquilibriumQuantities();
	void CopyBulkQuantities(CvolumeElement2D *element);
	void CopyEquilibriumQuantities(CvolumeElement2D *element);
	double P,epsilon,lambda; // prefactor used for viscous corrections
	void FillOutShearTensor(double &pixx,double &pixy,double &pixz,double &piyy,double &piyz,double &pizz);
	void Print();
	void CalcOmegamax();
	int MakeParts();
	int MakeParts_UniformXY();
	static Csampler *sampler;
};

class Cvertex2D{
public:
	double r[3],ux,uy;
	void Print(){
		printf("vertex: r=(%g,%g,%g), u=(%g,%g,%g)\n",r[0],r[1],r[2],sqrt(1.0+ux*ux+uy*uy),ux,uy);
	};
};

class Csampler{
public:
	CRandom *randy;
	CResList *reslist;
	vector<CvolumeElement2D> volume_element;
	vector<Cvertex2D> vertex;
	bool VISCOUSCORRECTIONS;
	FILE *xyfptr;

	Csampler(CB3D *b3d); // Constructor
	~Csampler();

	void ReadVolumeElements2D();
	void ReadVolumeElements3D();
	vector<CvolumeElement2D> element;
	int MakeB3DEvent();
	int MakeB3DEvent_UniformXY();
	double cummulative_N,cummulative_random,DELETA;
	int NBOSE;
	void SetVolumeElementsByHand(double T,int nelements_set,double elementvolume);
	void SetPiByHand(double pixx,double pixy,double pixz,double piyy,double piyz,double pizz);

	double GetLambda(); // Calculates lambda in terms of T and densities
	void CalcPiFromParts();
	double Tf,epsilonf,Pf,sf,lambdaf,nhadronsf;
	vector<double> densityf,boseweight;
	double GetLambda(double T,double P,double epsilon);
	CB3D *b3d;

	int nevents,ievent,nparts,nelements,nvertices;
};

#endif
