#include "binning.h"

double thetaBinsXW[nThetaBinsXW+1];
double binsWidth[nBinsWidth+1];
double binsGoodness[nBinsGoodness+1];
double binsDQDX[nBinsDQDX+1];
double binsMult[nBinsMult+1];
double binsIntegral[nBinsIntegral+1];
double binsPitch[nBinsPitch+1];
double binsX[4][nBinsX+1]; 
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

    double xMin[4] = {200, 50, 50, 200}; // Histo X min tpc 0-1-2-3 = EE-EW-WE-WW
    double xMax[4] = {370, 220, 220, 370}; // Histo X max tpc 0-1-2-3 = EE-EW-WE-WW

    for (int tpc = 0; tpc < 4; tpc++) {
        for (int i = 0; i <= nBinsX; i++)
            binsX[tpc][i] = xMin[tpc] + i * (xMax[tpc] - xMin[tpc]) / nBinsX;
    }
}

void getAnalysisEdges_XTheta(int tpc, std::vector<double>& xedges,
                           std::vector<double>& thetaEdges, std::vector<double>& zedges,
                           const TH3D* h3_orig)
{
    if (!h3_orig) return;

    xedges.clear();

    int nX = h3_orig->GetNbinsX();
    std::vector<double> allXedges;
    for (int i=1; i<=nX+1; ++i) allXedges.push_back(h3_orig->GetXaxis()->GetBinLowEdge(i));


    // Note: when I make the TH3 i have a conservative binning in X, 40 bins from xMin to xMax
    // this corresponds to (370-200)/40 = 4.25 cm bins.
    // Now, the "true" edges should be around 210 and 58.5 / 360.
    // here I am removing the first and last two bins, so that the start/end edges are 208.5 and 361.5.
    // I then merge toghether the first two and the last two bins.
    // ----> basically in this way I have 8.5 cm bins around anode and around cathode.

    int iStart=0, iEnd=allXedges.size();
    xedges.push_back(allXedges[iStart+2]);
    xedges.push_back(allXedges[iStart+4]);
    for (int i=iStart+5;i<=iEnd-4;i++) xedges.push_back(allXedges[i]);
    xedges.push_back(allXedges[iEnd-2]);

    // --- Theta edges ---
    // Here I am rebinning high angle bins
    double thetaArr[] = {0,4,8,12,16,20,24,28,32,36,40,44,48,54,65,90};
    thetaEdges.assign(thetaArr, thetaArr + sizeof(thetaArr)/sizeof(double));

    // --- Z edges: no changes
    int nZ = h3_orig->GetNbinsZ();
    zedges.clear();
    for (int i=1;i<=nZ+1;i++)
        zedges.push_back(h3_orig->GetZaxis()->GetBinLowEdge(i));
}

void getAnalysisEdges_XThetaUnBinned(int tpc,
                                    std::vector<double>& xedges,
                                    std::vector<double>& thetaEdges)
{
    xedges.clear();
    thetaEdges.clear();

    std::vector<double> allXedges;
    for (int i = 0; i <= nBinsX; ++i)
        allXedges.push_back(binsX[tpc][i]);

    int iStart = 0;
    int iEnd   = allXedges.size();

    // Remove first 2 and last 2 bins, merge first two and last two
    xedges.push_back(allXedges[iStart+2]);
    xedges.push_back(allXedges[iStart+4]);

    for (int i = iStart+5; i <= iEnd-4; ++i)
        xedges.push_back(allXedges[i]);

    xedges.push_back(allXedges[iEnd-2]);

    double thetaArr[] = {0,4,8,12,16,20,24,28,32,36,40,44,48,54,65,90};
    thetaEdges.assign(thetaArr, thetaArr + sizeof(thetaArr)/sizeof(double));
}