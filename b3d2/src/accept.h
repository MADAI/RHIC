#pragma once

#include "b3d.h"

/// Reprecent whether a particle has been accepted, and with which efficiency.
struct Acceptance {
    bool accept;
    double efficiency;

    Acceptance(bool accept=false, double efficiency=0) {
        this->accept = accept;
        this->efficiency = efficiency;
    }
};

/// Calculate the acceptance at STAR.
/// The possible values for `centrality` are the following:
///  0 =>  0- 5%
///  1 =>  5-10%
///  2 => 10-20%
///  3 => 20-30%
///  4 => 30-40%
///  5 => 40-50%
///  6 => 50-60%
///  7 => 60-70%
///  8 => 70-80%
// TODO: figure out WTH dca is.
Acceptance calc_acceptance_STAR(const CPart& part, int centrality, const double* dca);
