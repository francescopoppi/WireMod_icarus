#include "Corrections.h"

#include <cmath>
#include <iostream>

Corrections::Corrections(std::map<int,float>& EE,
                         std::map<int,float>& EW,
                         std::map<int,float>& WE,
                         std::map<int,float>& WW)
    : lifetimeEE(EE), lifetimeEW(EW), lifetimeWE(WE), lifetimeWW(WW)
{}

float Corrections::getLifetime(int run, int cryo, int rawTPC) const
{
    if (cryo == 0 && (rawTPC == 0 || rawTPC == 1)) {
        auto it = lifetimeEE.find(run);
        if (it != lifetimeEE.end()) return it->second;
    }
    if (cryo == 0 && (rawTPC == 2 || rawTPC == 3)) {
        auto it = lifetimeEW.find(run);
        if (it != lifetimeEW.end()) return it->second;
    }
    if (cryo == 1 && (rawTPC == 0 || rawTPC == 1)) {
        auto it = lifetimeWE.find(run);
        if (it != lifetimeWE.end()) return it->second;
    }
    if (cryo == 1 && (rawTPC == 2 || rawTPC == 3)) {
        auto it = lifetimeWW.find(run);
        if (it != lifetimeWW.end()) return it->second;
    }

    std::cout << "Lifetime not found for run " << run << "\n";
    return 1e9; 
}

double Corrections::computeDriftTime(double tic, double t0) const
{
    // 0.4 µs per tick, 340 us offset, t0 in from ns to us
    return (tic * 0.4 - 340.0) - t0 / 1.e3;
}

double Corrections::lifetimeCorr(double driftTime, double lifetime) const
{
    return std::exp(driftTime / lifetime);
}

double Corrections::correctCharge(double q, double tic, double t0,
                                  int run, int cryo, int rawTPC) const
{
    double driftT   = computeDriftTime(tic, t0);
    double lifetime = getLifetime(run, cryo, rawTPC);
    double corr     = lifetimeCorr(driftT, lifetime);

    return q * corr;
}
