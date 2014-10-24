#pragma once

#include <vector>
#include <cstdint>
#include <cassert>

template <typename T>
class Histogram {
  public:
    T x_min;
    T x_max;
    size_t nbins;
    std::vector<int64_t> histogram;

    /// Initialize histogram.
    Histogram(size_t nbins, T x_min, T x_max);
    /// Add value x to given histogram.
    void add(T x, int64_t count=1);
    /// Return the count of the bin corresponding to x.
    int64_t get_count(T x);
  private:
    double get_bin(T x);
};

template <typename T>
Histogram<T>::Histogram(size_t nbins, T x_min, T x_max) {
    this->x_min = x_min;
    this->x_max = x_max;
    this->nbins = nbins;
    this->histogram.resize(nbins); // zero by default
}

template <typename T>
double Histogram<T>::get_bin(T x) {
    return (x - x_min) / ((x_max - x_min) / nbins);
}

template <typename T>
void Histogram<T>::add(T x, int64_t count) {
    const double bin_float = get_bin(x);
    const size_t bin_int = (size_t) bin_float;
    //printf("x=%10g, bin=%10i\n", x, bin_float);
    if ( (bin_float >= 0) && (bin_int < nbins) )
        histogram[bin_int] += count;
    else
        printf("warning: %g out bounds of histogram!\n", x);
}

template <typename T>
int64_t Histogram<T>::get_count(T x) {
    const double bin_float = get_bin(x);
    assert(bin_float >= 0);
    const size_t bin_int = (size_t) bin_float;
    assert(bin_int < nbins);
    return histogram[bin_int];
}
