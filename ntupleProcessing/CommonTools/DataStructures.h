#pragma once
#include <vector>

struct Hit {
    float x, y, z;
    float integral, integral_c;
    float dqdx, dqdx_c;
    float width, goodness, pitch, rr;
    float dirx, diry, dirz;
    float tic;
    int   tpc, mult, wire, channel;
    bool  ontraj;
};

struct Track {
    int run, event, cryo;
    int selected, whicht0, G4ID;
    float t0PFP, t0CRTHit;
    float length;
    float startx, starty, startz;
    float endx, endy, endz;
    float dirx, diry, dirz;
    float PCAdirx, PCAdiry, PCAdirz;

    std::vector<Hit> hits[3];  // 0,1,2 planes
};
