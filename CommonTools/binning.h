#ifndef BINNING_H
#define BINNING_H

const int nBinsWidth      = 500;  
const int nBinsGoodness   = 400;  
const int nBinsDQDX       = 100;  
const int nBinsMult       = 5;    
const int nBinsIntegral   = 200;  
const int nBinsPitch      = 100;  
const int nBinsX = 40;
const int nBinsY = 31;
const int nBinsZ = 180;
const int nThetaBinsXW = 45; 

double thetaBinsXW[nThetaBinsXW+1];
double binsWidth[nBinsWidth+1];
double binsGoodness[nBinsGoodness+1];
double binsDQDX[nBinsDQDX+1];
double binsMult[nBinsMult+1];
double binsIntegral[nBinsIntegral+1];
double binsPitch[nBinsPitch+1];

double binsY[nBinsY+1];
double binsZ[nBinsZ+1];

inline void initBins() {

    for (int i = 0; i <= nThetaBinsXW; i++)
        thetaBinsXW[i] = 0.0 + i * (90.0 / nThetaBinsXW);

    for (int i = 0; i <= nBinsWidth; ++i)
        binsWidth[i] = 1.0 + i * (30.0 - 0.0) / nBinsWidth;

    for (int i = 0; i <= nBinsGoodness; ++i)
        binsGoodness[i] = 0.0 + i * (4.0 - 0.0) / nBinsGoodness;

    for (int i = 0; i <= nBinsDQDX; ++i)
        binsDQDX[i] = 0.0 + i * (2000.0 - 0.0) / nBinsDQDX;

    for (int i = 0; i <= nBinsMult; ++i)
        binsMult[i] = 0.0 + i * (5.0 - 0.0) / nBinsMult;

    for (int i = 0; i <= nBinsIntegral; ++i)
        binsIntegral[i] = 0.0 + i * (3000.0 - 0.0) / nBinsIntegral;

    for (int i = 0; i <= nBinsPitch; ++i)
        binsPitch[i] = 0.0 + i * (5.0 - 0.0) / nBinsPitch;
    
    for (int i = 0; i <= nBinsY; ++i)
        binsY[i] = -180 + i * (180+120) / nBinsY;

    for (int i = 0; i <= nBinsZ; ++i)
        binsZ[i] = -900 + i * 1800 / nBinsZ;
}

#endif
