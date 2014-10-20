#include <cstdio>
#include "hist.h"

Histogram::Histogram(size_t nbins, double x_min, double x_max) {
    this->x_min = x_min;
    this->x_max = x_max;
    this->nbins = nbins;
    this->histogram.resize(nbins);
}

void Histogram::add(double x) {
    const double bin_float = (x - x_min) / ((x_max - x_min) / nbins);
    const size_t bin_int = (size_t) bin_float;
    //printf("x=%10g, bin=%10i\n", x, bin);
    if ( (bin_float >= 0) && (bin_int < nbins) )
        histogram[bin_int]++;
    else
        printf("warning: %g out bounds of histogram!\n", x);
}

