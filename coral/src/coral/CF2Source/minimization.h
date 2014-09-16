#ifndef  MINIMIZATION_H 
#define  MINIMIZATION_H 
using namespace std;

// ABSTRACT class, can't be instantiated to be an object itself. Only good for a base(parent) class. 
class  CMinimization{ 
	public : 
	int  n;  //dimension 
	double  TOL_CG, TOL_Brent; 
	double  *  vec_x;  // variables 
	double  *  vec_dx;  // derivatives of variables 
	double  *  *  H;  // Hessian matrix as of d2f/dxidxj 
	virtual   double  fn( double  * vec_x)=0;  // get function's value with the current *x value. 
	virtual   bool  dfn( double  * vec_x)=0; // calculate the gradient of the function and the value is stored in vec_dx; 

	double  fn1( double  x);  //get function value with the current x value, this is one dimensional function. 
	// get function value at vec_x+x*vec_dir. vec_x is the starting point and x*vec_dir is the step x along vec_dir direction. 
	// This one-dimensional function will use fn(double *x) to evaluate function and used in one-dimensional minimization search. 
	// As an abstract class, this will not be used. Only for demonstration purpose and test. 
	CMinimization( int  dimension); 
	// default constructor. Only to define the object. Dimension is 1.
	//Use SetDimension() to set multi-dimensional search. 
	CMinimization(); 
	~CMinimization(); 

	bool  SetDimension( int  dimension); 

	// One dimensional search 
	// Given distinct points ax and bx, searches for new points ax, bx and cx so that they bracket a minimum of the function. 
	bool  bracket( double  &ax,  double  &bx,  double  &cx,  double  &fa,  double  &fb,  double  &fc); 
	double  brent( double  ax,  double  bx,  double  cx,  double  &xmin); 

	bool  conjugate_gradient( double  * initial_x,  int  &iteration,  double  &fmin);  // input the initial guess of x 
	// return the minimum point by reset initial_x and return minimum fmin by reference. Times of iteration is returned by iteration. 
	// the minimum point can also be found from vec_x; 
	protected : 
	bool  AllocVectors(); 
	bool  TagMultiD; 
// "linear_" means only used in one-dimensional minimization search 
	bool  linear_Min( double  &fmin);  // fmin returns the minimum value of function. vec_x is reset to the local minimum. 
	double  *  linear_dir;  // direction used in one-dimensional minimization search 
	double  *  linear_trial_x;  // To avoid using new operator each time when try new x 

// book keeping utilities 
	inline  void  shift2( double  &a,  double  &b,  double  c) 
	{ 
		a=b; 
		b=c; 
	} 
	inline  void  shift3( double  &a,  double  &b,  double  &c,  double  d) 
	{  
		a=b; 
		b=c; 
		c=d; 
	} 
	inline  void  swap2( double  &a,  double  &b) 
	{ 
		double  c=b; 
		b=a; 
		a=c; 
	} 

}; 

#endif

