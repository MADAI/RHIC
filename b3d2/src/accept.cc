#include <math.h>
#include <cassert>

#include "accept.h"
#include "pow.h"

/// Translate the general particle ID to the specific STAR PID.
/// Returns -1 if there is not corresponding STAR PID.
/// Exits if the PID is not in the STAR list.
int get_starpid(int pid) {
    switch (pid) {
        case 211: case -211:
            return 1;
        case 321: case -321:
            return 2;
        case -2212:
            return 3;
        case 2212:
            return 4;
        default:
            if (abs(pid) != 2112 && abs(pid) != 311 && abs(pid) != 111 && abs(pid) != 22) {
                printf("pid=%d isn't in STAR list\n", pid);
                exit(1);
            }
            return -1;
    }
}

/// Routine to return acceptance and efficiency
/// Gary D. Westfall
/// January, 2014
/// Uses results published by STAR PRC 79 034909 2009, page 034909-13 for functional forms
/// Uses embedding results from Hui Wang for Run 10, 200 GeV
/// Determined TOF matching efficiency looking at pion spectra
///  and double checked with kaons and protons
/// Input
///  pid=1 pion (either sign)
///  pid=2 kaon (with sign)
///  pid=3 anti-proton
///  pid=4 proton
///  cen=0 0-5%
///  cen=1 5-10%
///  cen=2 10-20%
///  cen=3 20-30%
///  cen=4 30-40%
///  cen=5 40-50%
///  cen=6 50-60%
///  cen=7 60-70%
///  cen=8 70-80%
///  basic acceptance cuts are
///   |eta| < 1.0 for TPC only
///   |eta| < 0.9 for TPC+TOF
///   full phi acceptance
///  pt dependence is calculated in detail
Acceptance star_acc_eff(int pid,double pt,double eta,double phi,int cen){
	int i;
	float f1,f2,f3,f4,f5;
	float P0,P1,P2,P3,P4,P5,P6;
	float p,theta;
	// Pion parameters
	float P0_pion[9]={0.714,0.745,0.779,0.810,0.830,0.842,0.851,0.862,0.874};
	float P1_pion[9]={0.083,0.082,0.079,0.087,0.105,0.085,0.121,0.109,0.121};
	float P2_pion[9]={1.906,2.021,2.031,2.331,2.915,2.498,3.813,3.052,3.287};
	// Kaon parameters
	float P0_kaon[9]={0.634,0.646,0.691,0.711,0.748,0.742,0.783,0.758,0.748};
	float P1_kaon[9]={0.226,0.220,0.215,0.205,0.204,0.194,0.202,0.198,0.218};
	float P2_kaon[9]={1.387,1.490,1.447,1.533,1.458,1.485,1.372,1.559,1.665};
	float P3_kaon[9]={0.021,0.026,0.021,0.025,0.019,0.027,0.016,0.029,0.042};
	// pbar parameters
	float P0_pbar[9]={0.707,0.735,0.770,0.803,0.821,0.836,0.849,0.891,0.889};
	float P1_pbar[9]={0.191,0.185,0.196,0.194,0.191,0.200,0.199,0.156,0.197};
	float P2_pbar[9]={2.440,2.695,3.336,3.684,4.021,5.037,4.850,2.701,4.538};
	// p parameters
	float P0_p[9]={0.720,0.750,0.785,0.818,0.837,0.852,0.865,0.885,0.921};
	float P1_p[9]={0.201,0.197,0.201,0.193,0.202,0.202,0.192,0.192,0.178};
	float P2_p[9]={3.247,3.467,4.124,4.267,5.414,5.621,5.080,4.950,3.957};


	// Check the input data

	if ( (pid < 1) || (pid > 4)
	     || (pt < 0.2) || (pt > 3.0)
	     || (fabs(eta) > 1.0)
	     || (cen < 0) || (cen > 8) ) {
        return Acceptance(false, 0);
    }

	// The basic parameters are OK, check the details

    bool accept=false;
	double eff=0.0;

    switch (pid) {
        case 1: // Pions
            theta=2.0*atan(exp(-eta));
            p=pt/sin(theta);
            if((pt > 0.2) && (p < 1.6)){
                accept=true;
                P0=P0_pion[cen];
                P1=P1_pion[cen];
                P2=P2_pion[cen];
                f1=P1/pt;
                eff=P0*exp(-pow(f1,P2));
                if(pt > 0.60){
                    if(fabs(eta) < 0.9){
                        eff=0.59*eff;
                    } else {
                        eff=0.0;
                    }
                }
            }
            break;

        case 2: // Kaons
            theta=2.0*atan(exp(-eta));
            p=pt/sin(theta);
            if((pt > 0.20) && (p < 1.6)){
                accept=true;
                P0=P0_kaon[cen];
                P1=P1_kaon[cen];
                P2=P2_kaon[cen];
                P3=P3_kaon[cen];
                f1=P1/pt;
                eff=P0*exp(-pow(f1,P2))+P3*pt;
                if(pt > 0.6){
                    if(fabs(eta) < 0.9){
                        eff=0.59*eff;
                    } else {
                        eff=0.0;
                    }
                }
            }
            break;

        case 3: // pbar
            theta=2.0*atan(exp(-eta));
            p=pt/sin(theta);
            if((pt > 0.4) && (p < 3.0)){
                accept=true;
                P0=P0_pbar[cen];
                P1=P1_pbar[cen];
                P2=P2_pbar[cen];
                f1=P1/pt;
                eff=P0*exp(-pow(f1,P2));
                if(pt > 1.0){
                    if(fabs(eta) < 0.9){
                        eff=0.59*eff;
                    } else {
                        eff=0.0;
                    }
                }
            }
            break;

        case 4: // p
            theta=2.0*atan(exp(-eta));
            p=pt/sin(theta);
            if((pt > 0.4) && (p < 3.0)){
                accept=true;
                P0=P0_p[cen];
                P1=P1_p[cen];
                P2=P2_p[cen];
                f1=P1/pt;
                eff=P0*exp(-pow(f1,P2));
                if(pt > 1.0){
                    if(fabs(eta) < 0.9){
                        eff=0.59*eff;
                    } else {
                        eff=0.0;
                    }
                }
            }
            break;
        default:
            assert(0);
            // unreachable
    }

    return Acceptance(accept, eff);
}

Acceptance calc_acceptance_STAR(const CPart& part, int centrality, const double* dca) {
    // Those are the cuts independent of the STAR acceptance.
    // They might be more restrictive.
    const double ETAMIN = -0.9;
    const double ETAMAX = 0.9;
    const double PTMIN = 200;
    const double PTMAX = 1600;
    //const int CENTRALITY = 4;

    const int starpid = get_starpid(part.resinfo->code);
    if (starpid == -1)
        return Acceptance(false, 0);
    const auto& p = part.p;
    const double pt = sqrt(square(p[1]) + square(p[2]));
    const double pmag = sqrt(square(pt) + square(p[3]));
    const double eta = atanh(p[3]/pmag);

    if (ETAMIN < eta && eta < ETAMAX
        && dca[0] < 2.0
        && PTMIN < pt && pt < PTMAX) {
        const double phi = atan2(p[2], p[1]);
        return star_acc_eff(starpid, 0.001*pt, eta, phi, centrality);
    }

    return Acceptance(false, 0);
}

