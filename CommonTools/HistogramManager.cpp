#include "HistogramManager.h"
#include "binning.h"

#include <TMath.h>
#include <cmath>

HistogramManager::HistogramManager()
{
    initBins();
    bookHistograms();
}

void HistogramManager::bookHistograms()
{
    const char* vars[nVar] = {"width","goodness","hit_integral","hit_pitch"};

    double xMin[4] = {200, 50, 50, 200};
    double xMax[4] = {370, 220, 220, 370};

    for (int tpc=0; tpc<nTPC; tpc++) {

        double xBins[nBinsX+1];
        for (int i=0;i<=nBinsX;i++)
            xBins[i] = xMin[tpc] + i*(xMax[tpc]-xMin[tpc])/nBinsX;

        for (int plane=0; plane<nPlane; plane++) {
            for (int v=0; v<nVar; v++) {

                int nBinsVar = (v==0? nBinsWidth :
                                v==1? nBinsGoodness :
                                v==2? nBinsIntegral : nBinsPitch);

                double* binsVar = (v==0? binsWidth :
                                   v==1? binsGoodness :
                                   v==2? binsIntegral : binsPitch);

                TString name1  = Form("h3D_%s_TPC%d_plane%d_XTheta", vars[v], tpc, plane);
                TString title1 = Form("%s (TPC%d plane%d);x;thetaXW;%s", vars[v], tpc, plane, vars[v]);

                h3D_XTheta[tpc][plane][v] = new TH3D(name1,title1,
                                                    nBinsX,xBins,
                                                    nThetaBinsXW,thetaBinsXW,
                                                    nBinsVar,binsVar);

                TString name2  = Form("h3D_%s_TPC%d_plane%d_XTheta_c", vars[v], tpc, plane);
                h3D_XTheta_c[tpc][plane][v] = new TH3D(name2,title1,
                                                      nBinsX,xBins,
                                                      nThetaBinsXW,thetaBinsXW,
                                                      nBinsVar,binsVar);

                TString name3  = Form("h3D_%s_TPC%d_plane%d_YZ", vars[v], tpc, plane);
                TString title3 = Form("%s (TPC%d plane%d);Y;Z;%s", vars[v], tpc, plane, vars[v]);

                h3D_YZ[tpc][plane][v] = new TH3D(name3,title3,
                                                nBinsY,binsY,
                                                nBinsZ,binsZ,
                                                nBinsVar,binsVar);

                TString name4  = Form("h3D_%s_TPC%d_plane%d_YZ_c", vars[v], tpc, plane);
                h3D_YZ_c[tpc][plane][v] = new TH3D(name4,title3,
                                                  nBinsY,binsY,
                                                  nBinsZ,binsZ,
                                                  nBinsVar,binsVar);
            }
        }
    }
}


void HistogramManager::fill(const Hit& hit, int plane, double thetaXW)
{
    int tpc = hit.tpc;
    if (tpc < 0 || tpc >= nTPC) return;
    if (plane < 0 || plane >= nPlane) return;

    double x = std::fabs(hit.x);

    h3D_XTheta[tpc][plane][0]->Fill(x, thetaXW, hit.width);
    h3D_XTheta[tpc][plane][1]->Fill(x, thetaXW, hit.goodness);
    h3D_XTheta[tpc][plane][2]->Fill(x, thetaXW, hit.integral);
    h3D_XTheta[tpc][plane][3]->Fill(x, thetaXW, hit.pitch);

    h3D_XTheta_c[tpc][plane][0]->Fill(x, thetaXW, hit.width);
    h3D_XTheta_c[tpc][plane][1]->Fill(x, thetaXW, hit.goodness);
    h3D_XTheta_c[tpc][plane][2]->Fill(x, thetaXW, hit.integral_c);
    h3D_XTheta_c[tpc][plane][3]->Fill(x, thetaXW, hit.pitch);

    // YZ
    h3D_YZ[tpc][plane][0]->Fill(hit.y, hit.z, hit.width);
    h3D_YZ[tpc][plane][1]->Fill(hit.y, hit.z, hit.goodness);
    h3D_YZ[tpc][plane][2]->Fill(hit.y, hit.z, hit.integral);
    h3D_YZ[tpc][plane][3]->Fill(hit.y, hit.z, hit.pitch);

    h3D_YZ_c[tpc][plane][0]->Fill(hit.y, hit.z, hit.width);
    h3D_YZ_c[tpc][plane][1]->Fill(hit.y, hit.z, hit.goodness);
    h3D_YZ_c[tpc][plane][2]->Fill(hit.y, hit.z, hit.integral_c);
    h3D_YZ_c[tpc][plane][3]->Fill(hit.y, hit.z, hit.pitch);
}


void HistogramManager::write(TFile& fout)
{
    fout.cd();

    for (int t=0;t<nTPC;t++)
        for (int p=0;p<nPlane;p++)
            for (int v=0;v<nVar;v++) {
                h3D_XTheta[t][p][v]->Write();
                h3D_XTheta_c[t][p][v]->Write();
                h3D_YZ[t][p][v]->Write();
                h3D_YZ_c[t][p][v]->Write();
            }
}
