#include <iostream>
#include <fstream>
#include <string>
#include "CEos.h"

int main (int argc, char * const argv[]) {

	parameterMap* pMap = new parameterMap();
	for (int i=1;i<argc;i++) 
		parameter::ReadParsFromFile(*pMap, argv[i]);
	
	CEos mEos(pMap);
	mEos.printEos("mEos.txt");
	
}