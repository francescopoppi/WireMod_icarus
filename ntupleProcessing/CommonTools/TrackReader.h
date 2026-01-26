#pragma once
#include "DataStructures.h"
#include <TFile.h>
#include <vector>
#include <map>

class TrackReader {
public:
    TrackReader(std::map<int,float>& EE,
                std::map<int,float>& EW,
                std::map<int,float>& WE,
                std::map<int,float>& WW);

    std::vector<Track> readFile(TFile* f, int cryo, bool applyTruthCuts);

private:
    float getLifetime(int run, int cryo, int tpc);
    std::map<int,float> lifetimeEE, lifetimeEW, lifetimeWE, lifetimeWW;
};
