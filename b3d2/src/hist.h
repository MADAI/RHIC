#include <vector>

class Histogram {
  public:
    double x_min;
    double x_max;
    size_t nbins;
    std::vector<unsigned long> histogram;

    /// Initialize histogram.
    Histogram(size_t nbins, double x_min, double x_max);
    /// Add value x to given histogram.
    void add(double x);
};
