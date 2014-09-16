#ifndef CHEBYSHEV_FIT_H
#define CHEBYSHEV_FIT_H

#include "chebyshev.h"
#include "tnt_array1d.h"
#include "tnt_array2d.h"

using namespace TNT;

class CChebyshevDataFit: public CChebyshevApprox1D {

public:
    CChebyshevDataFit(void){RunQuiet();}
    ~CChebyshevDataFit(void){}
    
    /* Error on fit */
    double Error(double x) const;
    /* Set up chi2 fit*/
    void ChiSquareFitCoeffs(
        const Array1D<double> &x, 
        const Array1D<double> &y,
        const Array1D<double> &dy
        );
    void ChiSquareFitCoeffs(
        const Array1D<double> &x, 
        const Array1D<double> &y,
        const Array2D<double> &covy
        );

private:
    Array2D<double> d2c;

};


#endif
