#ifndef BINNING_H
#define BINNING_H

#include <vector>
#include <TH3.h> 

const int nBinsWidth      = 500;  
const int nBinsGoodness   = 400;  
const int nBinsDQDX       = 100;  
const int nBinsMult       = 5;    
const int nBinsIntegral   = 200;  
const int nBinsPitch      = 100;  
const int nBinsX = 40;
const int nBinsY = 31;
const int nBinsZ = 180;
const int nBinsY_2 = 15;
const int nBinsZ_2 = 45;
const int nThetaBinsXW = 45; 

extern double binsX[4][nBinsX+1]; // tpc 0-1-2-3 = EE-EW-WE-WW
extern double thetaBinsXW[nThetaBinsXW+1];
extern double binsWidth[nBinsWidth+1];
extern double binsGoodness[nBinsGoodness+1];
extern double binsDQDX[nBinsDQDX+1];
extern double binsMult[nBinsMult+1];
extern double binsIntegral[nBinsIntegral+1];
extern double binsPitch[nBinsPitch+1];
extern double binsY[nBinsY+1];
extern double binsZ[nBinsZ+1];
extern double binsY_2[nBinsY_2+1];
extern double binsZ_2[nBinsZ_2+1];


void initBins();

void getAnalysisEdges_XTheta(int tpc, std::vector<double>& xedges,
                           std::vector<double>& thetaEdges, std::vector<double>& zedges,
                           const TH3D* h3_orig);

void getAnalysisEdges_XThetaUnBinned(int tpc,
                                    std::vector<double>& xedges,
                                    std::vector<double>& thetaEdges);

#endif
