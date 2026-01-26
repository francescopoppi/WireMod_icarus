#pragma once

#include "DataStructures.h"

#include <TH3D.h>
#include <TFile.h>
#include <TString.h>

class HistogramManager {
public:
    HistogramManager();

    void fill(const Hit& hit, int plane, double thetaXW);
    void write(TFile& fout);

private:
    static const int nTPC   = 4;
    static const int nPlane = 3;
    static const int nVar   = 4; // width, goodness, integral, pitch

    TH3D* h3D_XTheta[nTPC][nPlane][nVar];
    TH3D* h3D_XTheta_c[nTPC][nPlane][nVar];
    TH3D* h3D_YZ[nTPC][nPlane][nVar];
    TH3D* h3D_YZ_c[nTPC][nPlane][nVar];

    void bookHistograms();
};
