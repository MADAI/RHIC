#ifndef XGRAPH_INCLUDE_H
#define XGRAPH_INCLUDE_H
#include <cstdio>
#include <iostream>
#include <X11/Xlib.h>
#include <cmath>
#include <cstring>
#include <cstdlib>
//#include <dlfcn.h>
#include <unistd.h>
//#include <ltdl.h>

using namespace std;

class CAxesInfo{
 public:
  double xmin,ymin,xmax,ymax;
  double xscale,yscale;
  double xintercept; //where y-axis intercepts x axis
  double yintercept; //where x-axis intercepts y axis
  int nxtics;
  int nytics;
  CAxesInfo();
};

class CXGraph{
 public:
  Display *theDisplay;
  Window theWindow;
  Screen *theScreen;
  Font theFont;
  GC theGC;
	int window_horizsize,window_vertsize,ixposition,iyposition;
  Colormap screen_colormap;
  XColor red, brown, blue, yellow, green, orange, black, cyan, violet, lightblue, navy, pink,white;
  Status rc;
  CAxesInfo axesinfo;
  void setcolor(string color);
  void setaxes(double xmin,double ymin,double xmax,double ymax);
  void drawtext(char *string,double x,double y);
  void drawline(double x1,double y1,double x2,double y2);
  void drawcircle(double x,double y,double size);
  void drawsquare(double x,double y,double size);
  void drawdiamond(double x,double y,double size);
  void drawuptriangle(double x,double y,double size);
  void drawdowntriangle(double x,double y,double size);
  void drawpoint(double x,double y);
  void drawarrow(double x1,double y1,double x2,double y2,double headsize);
  void drawaxes();
  void closedisplay();
  void getij(double x,double y,int *i,int *j);
  void plotline(double *x,double *y,int npts);
  void plotpoints(double *x,double *y,int npts);
	void clear();
	void flush();
  CXGraph(int horizsize,int vertsize,int ixposition,int iyposition);
	~CXGraph();
};
#endif
