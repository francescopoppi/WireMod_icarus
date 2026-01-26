#pragma once

#include <map>

class Corrections {
public:
    Corrections(std::map<int,float>& EE,
                std::map<int,float>& EW,
                std::map<int,float>& WE,
                std::map<int,float>& WW);

    float getLifetime(int run, int cryo, int rawTPC) const;

    // drift time in ms
    double computeDriftTime(double tic, double t0) const;

    // lifetime correction factor
    double lifetimeCorr(double driftTime, double lifetime) const;

    double correctCharge(double q, double tic, double t0,
                         int run, int cryo, int rawTPC) const;

private:
    std::map<int,float>& lifetimeEE;
    std::map<int,float>& lifetimeEW;
    std::map<int,float>& lifetimeWE;
    std::map<int,float>& lifetimeWW;
};
