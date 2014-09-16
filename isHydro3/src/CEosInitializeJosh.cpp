/*
*  CEos.cpp
*  Created by Joshua Vredevoogd on 3/4/09.
*/

#include "CEos.h"

//constructor - read from EOS interpolation
void CEosCalculator::InitializeEosJosh(parameterMap* pM){
		
	if (parameter::getS(*pM,"EQOFST_LATDATA","null") == string("null")){
		printf("\n** CEos requires parameterMap variable EQOFST_LATDATA if EQOFST_LATEOS == true **\nAborting...\n\n");
		exit(1);
	}
		
	string latFileName = parameter::getS(*pM,"EQOFST_LATDATA","none");
	if (!checkFile(latFileName)) {
		printf("\n** No file found at %s (EQOFST_LATDATA) **\nAborting...\n\n",latFileName.c_str());
		exit(1);
	}

	std::ifstream latFile(latFileName.c_str());
		
	double lowT  = parameter::getD(*pM,"HYDRO_FOTEMP",0.15);
	double highT = parameter::getD(*pM,"EQOFST_HIGH_MERGE_T",0.2);
		
	if (highT < lowT) {
		printf("Inconsistent values of HYDRO_FOTEMP and EQOFST_HIGH_CROSS_T.\n Aborting.\n");
		exit(1);
	}
	const int size=1000;
	double lT[size], lE[size], lP[size], lS[size];
	double hT[size], hE[size], hP[size], hS[size];
	int lSize=0, hSize=0;
		
	char test;
	latFile >> test;
	latFile.ignore(1000,'\n');
		
	// re-used counter -- gross!
	int i;
		
	// check first line to determine file format
	if (test == 'e') { // EOS from Pasi
		double junk;
		for(i=0; !latFile.eof() && i<size; i++) {
			latFile >> lE[i] >> lP[i] >> lS[i] >> junk;
			lT[i] = (lE[i] + lP[i])/lS[i];
		}
		lSize = i-1;
		latFile.close();
	}
	else { // EOS from Ron/Fodor
		double junk;
		for(i=0; !latFile.eof() && i<size; i++) {
			latFile >> lT[i] >>  junk >> lE[i] >> lP[i] >> lS[i] >> junk >> junk;
				
			lE[i]   *= pow(lT[i],4.)/pow(JHBARC,3.);
			lP[i]   *= pow(lT[i],4.)/pow(JHBARC,3.);
			lS[i]   *= pow(lT[i],3.)/pow(JHBARC,3.);
		}	
		lSize = i-1;
	}
	latFile.close();
		
	// get hrg file name from parameterMap
	string hrgFileName = parameter::getS(*pM,"EQOFST_HRGDATA","none");
		
	// not merging lattice with HRG data
	if (hrgFileName == string("none")){
		// move variables to permanent homes
		aSize = lSize;
		for (int j=0;j<aSize;j++) {
			temp[j] = lT[j]; 
			sd[j] = lS[j];
			ed[j] = lE[j];
			pr[j] = lP[j];
		} 
	}
	else if (!checkFile(hrgFileName)) {
		printf("Unable to find eos data at %s....\n Continuing with lattice data only....\n",
		hrgFileName.c_str());
			
		// move variables to permanent homes
		aSize = lSize;
		for (int j=0;j<aSize;j++) {
			temp[j] = lT[j]; 
			sd[j] = lS[j];
			ed[j] = lE[j];
			pr[j] = lP[j];
		} 
	}
	else {
		// Open hrg eos file for merging
		std::ifstream hrgFile;
		hrgFile.open(hrgFileName.c_str());
			
		// skip first lines of hrg file
		for (int j=0;j<8;j++) 
			hrgFile.ignore(1000,'\n');
			
		// read hrg data
		for(i=0; !hrgFile.eof() && i<size; i++) {
			hrgFile >> hT[i] >> hS[i] >> hP[i] >> hE[i];
			hT[i] /= 1000.; hP[i] /= 1000.; hE[i] /= 1000.;
			if (hT[i] < 0.02) i--;
		}
		hSize = i-1;	
			
		// merge will fail if 
		if (hT[hSize] < highT){
			if (hT[hSize] < lowT){
				printf("No hadron resonance gas data up to %g GeV (HYDRO_FOTEMP).\n",lowT);
				printf("Unclear how to proceed.  Aborting.\n");
			}
			else {
				printf("No hadron resonance gas data up to %g GeV (EQOFST_HIGH_CROSS_T).\n",highT);
				printf("Resetting to %g GeV, maximum value in %s (EQOFST_HRGDATA).\n",
				hT[hSize],parameter::getS(*pM,"EQOFST_LATDATA","none").c_str());
				highT = hT[hSize];
			}
		}
			
		// offset by tenth of step eliminating round-off error effects
		lowT  += (lT[1] - lT[0])/10.;
		highT += (lT[1] - lT[0])/10.;
			
		// find boundaries in each array
		int hlowTIndex=0, llowTIndex=0, highTIndex=0;
		while (hT[hlowTIndex] < lowT )
			hlowTIndex++;
		while (lT[llowTIndex] < lowT )
			llowTIndex++; 
		while (lT[highTIndex] < highT)
			highTIndex++;
		hlowTIndex--;
		llowTIndex--;
			
		// copy hrg data below TLow
		aSize=0;
		for (int j=0; j<=hlowTIndex; j++) {
			sd[j] = hS[j];
			temp[j] = hT[j];
			aSize++;
		}
			
		// merge between TLow and THigh
		for (int j=1;;j++) {
			temp[aSize] = lT[llowTIndex+j];
			if (temp[aSize] >= highT) break;
			double w = (tanh( tan((PI/2.)*( 2.*(highT-temp[aSize])/(highT-lowT) - 1.)))+1.)/2. ;
			sd[aSize] = w*hS[hlowTIndex+j] + (1.-w)*lS[llowTIndex+j];
			aSize++;
		}
		aSize--;

		// take lattice values above THigh
		for (int j=highTIndex-1;j<lSize;j++) {
			sd[aSize] = lS[j];
			temp[aSize] = lT[j];
			aSize++;
		}
		aSize--;
			
		// take hrg pressure or calculate by integral method
		for (int j=0;j<=aSize;j++) 
			if (j<= hlowTIndex)
				pr[j] = hP[j];
		else 
			pr[j] = pr[j-1] + 0.5*(temp[j]-temp[j-1])*(sd[j]+sd[j-1]);
			
		// energy density from therm identitiy
		for (int j=0;j<=aSize;j++) 
			ed[j] = temp[j]*sd[j] - pr[j];
			
	}
				
	// add bump to eos for statistical runs
	if ( parameter::getB(*pM,"EQOFST_ADD_BUMP",false) ) {
			
		// temperature to add bump to s/T^3
		//			double T0     = parameter::getD(*pM,"EQOFST_BUMP_START",-1.);
		double T0 = parameter::getD(*pM,"HYDRO_FOTEMP",0.15);
		// width (in temperature [GeV]) of bump
		double width  = parameter::getD(*pM,"EQOFST_BUMP_WIDTH",-1.);
		// height of bump to s/T^3 
		// given as fraction of maximum possible change (for cs2 to stay in bounds)
		double height = parameter::getD(*pM,"EQOFST_BUMP_HEIGHT",-1.); 
			
		//
		bExpTail   = false;
		bGaussTail = false;
		bPowerTail = true; 
			
		// do nothing if any value is not set
		if (height == -1. || T0 == -1. || width == -1.) {
			// print an error message if they have one but not all
			if (height != -1. || T0 != -1. || width != -1.){
				printf("\n *** ParameterMap contains insufficient information to modify the EoS ***\n");
				printf(" *** Required variables: EQOFST_BUMP_START EQOFST_BUMP_WIDTH EQOFST_BUMP_HEIGHT ***\n");
				printf("Continuing without modification....\n\n");
			}
			// else just ignore the ADD_BUMP flag
			return;
		}
		else if (height > 1. || height < 0.) {
			printf("\n *** EQOFST_BUMP_HEIGHT out of bounds (must be between 0 and 1) ***\n");
			printf("Continuing without modification....\n\n");
		}
			
		// modify height
		else addBump(T0,width,height);
	}
	/*
	for(int iT=0;iT<=1000;iT++){
		printf("%6.4f %10.6f %10.6f %10.6f\n",temp[iT],ed[iT],pr[iT],sd[iT]);
	}
	*/
	
	// so we don't have to find these each time
	eAtT = getEGivenT(svSwitchTemp);
	sAtT = getSGivenT(svSwitchTemp);
	
	//	printEos("mEos.txt");
}
