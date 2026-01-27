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

extern double thetaBinsXW[nThetaBinsXW+1];
extern double binsWidth[nBinsWidth+1];
extern double binsGoodness[nBinsGoodness+1];
extern double binsDQDX[nBinsDQDX+1];
extern double binsMult[nBinsMult+1];
extern double binsIntegral[nBinsIntegral+1];
extern double binsPitch[nBinsPitch+1];
extern double binsY[nBinsY+1];
extern double binsZ[nBinsZ+1];

void initBins();

#endif
