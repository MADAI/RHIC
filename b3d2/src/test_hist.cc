#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "hist.h"

double gaussrand() {
    static double V1, V2, S;
    static int phase = 0;
    double X;

    if (phase == 0) {
        do {
            double U1 = (double) rand() / RAND_MAX;
            double U2 = (double) rand() / RAND_MAX;

            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while (S >= 1 || S == 0);

        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;

    return X;
}

int main() {
    const size_t nbins = 100;
    auto h = Histogram<double>(nbins, -5, 5);
    for (size_t i = 0; i < 100000000; i++)
        h.add(gaussrand());
    for (size_t i = 0; i < nbins; i++)
        std::cout << h.histogram[i] << "  ";
    std::cout << std::endl;
    std::cout << "\n" << h.get_count(0) << std::endl;
}
