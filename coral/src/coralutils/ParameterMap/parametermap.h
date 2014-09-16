#ifndef parameterMap_H
#define parameterMap_H

#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>

using namespace std;

//---------------------------------------------------------------
//This code helps one parse the map that contains configuration 
// parameters.
//
//Example:
// 
// MyCoolClass(parameterMap m){  //m has the parameters and is passed in
// int importantParameter = parameter::getI(m,"nameOfParameter",-1);
//  ...
//
//  If "nameOfParameter is not a key in the map, -1 is returned since that
// is the 3rd argment of the getI function.
//
//MH 22 jun04
//---------------------------------------------------------------

//This code only works with a map of the type below.  The type def is
//to make it easy to remember.
typedef  map<string,string> parameterMap;

//These functions are all in the namespace parameter.
namespace parameter {
  bool   getB(parameterMap ,string ,bool);
  int    getI(parameterMap ,string ,int);
  string getS(parameterMap ,string ,string);
  double getD(parameterMap ,string ,double);
  vector< double > getV(parameterMap, string, string);
  vector< string > getVS(parameterMap, string, string);
  vector< vector< double > > getM(parameterMap, string, double);
  void set(parameterMap&, string, double);
  void set(parameterMap&, string, int);
  void set(parameterMap&, string, bool);
  void set(parameterMap&, string, string);
  void set(parameterMap&, string, char*);
  void set(parameterMap&, string, vector< double >);
  void set(parameterMap&, string, vector< string >);
  void set(parameterMap&, string, vector< vector< double > >);
  void ReadParsFromFile(parameterMap&, const char *filename);
  void ReadParsFromFile(parameterMap&, string filename);
  void PrintPars(parameterMap&);
};

#endif
