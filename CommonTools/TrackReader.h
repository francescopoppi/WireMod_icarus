#pragma once
#include "DataStructures.h"
#include "Corrections.h"
#include <TFile.h>
#include <vector>
#include <map>

class TrackReader {
public:
    TrackReader(Corrections& corr);

    std::vector<Track> readFile(TFile* f, int cryo, bool applyTruthCuts);

private:
    Corrections& corrections;
};
