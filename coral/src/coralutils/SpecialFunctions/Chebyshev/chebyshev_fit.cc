#include "chebyshev.h"
#include "tnt_array1d.h"
#include "tnt_array2d.h"
#include "chebyshev_fit.h"
#include "lsqrinvert.h"

using namespace TNT;

double CChebyshevDataFit::Error(double x) const{ 
    return 0.0;
}


void CChebyshevDataFit::ChiSquareFitCoeffs(const Array1D<double> &x, 
         const Array1D<double> &y,const Array1D<double> &dy){
         
   // Setup matrix we will invert
   int M=x.dim(), N=nc;
   Array2D<double> T(M,N);  
   
   // Compute center and range of the fit region
   double midpt=0.5*(uplim+lolim);
   double halfrange=0.5*(uplim-lolim);
   
   // Compute the matrix we want to invert
   {
      double xx;int i,j;
      for (i=0;i<M;i++){
         xx=(x[i]-midpt)/halfrange;
         T[i][0]=ChebyshevPolynomial(0,xx)-0.5;
         for (j=1;j<N;j++){
            T[i][j]=ChebyshevPolynomial(j,xx);
         }
      }
   }

   Array1D<double> new_c(N);
   Array2D<double> new_cov(N,N);

   LeastSquaresInvert(y,dy,T,new_c,new_cov);
   
   c=new_c;
   d2c=new_cov;
}         
         
void CChebyshevDataFit::ChiSquareFitCoeffs(const Array1D<double> &x, 
         const Array1D<double> &y,const Array2D<double> &covy){

   // Setup matrix we will invert
   int M=x.dim(), N=nc;
   Array2D<double> T(M,N);  
   
   // Compute center and range of the fit region
   double midpt=0.5*(uplim+lolim);
   double halfrange=0.5*(uplim-lolim);
   
   // Compute the matrix we want to invert
   {
      double xx;int i,j;
      for (i=0;i<M;i++){
         xx=(x[i]-midpt)/halfrange;
         T[i][0]=ChebyshevPolynomial(0,xx)-0.5;
         for (j=1;j<N;j++){
            T[i][j]=ChebyshevPolynomial(j,xx);
         }
      }
   }

   Array1D<double> new_c(N);
   Array2D<double> new_cov(N,N);

   LeastSquaresInvert(y,covy,T,new_c,new_cov);
   
   c=new_c;
   d2c=new_cov;
}
