#include "binning.h"

double thetaBinsXW[nThetaBinsXW+1];
double binsWidth[nBinsWidth+1];
double binsGoodness[nBinsGoodness+1];
double binsDQDX[nBinsDQDX+1];
double binsMult[nBinsMult+1];
double binsIntegral[nBinsIntegral+1];
double binsPitch[nBinsPitch+1];

double binsY[nBinsY+1];
double binsZ[nBinsZ+1];

void initBins() {
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
