
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TGraph2D.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>
#include <string>


#include "AnalysisTools.h"

void makeAll_ITM_TGraphs(
    const char* infile,
    TFile* fout,
    const char* subdirName,
    const char* plotOccDir,
    const char* mode,
    Double_t sig_down = -2, Double_t sig_up = 1.75,
    Double_t tol = 1e-4,
    Int_t minEntries = 50,
    bool savePlots = false)
{
    const char* vars[4] = {"width", "goodness", "hit_integral", "hit_pitch"};

    gStyle->SetOptStat(0);

    // --- Cuts
    TF1* cut_hit = makeHitIntegralCut(mode);
    TF1* cut_width[3] = { makeWidthCut(0), makeWidthCut(1), makeWidthCut(2) };

    TFile* fin = TFile::Open(infile);
    if (!fin || fin->IsZombie()) {
        std::cerr << "Input file not found" << infile << std::endl;
        return;
    }

    fout->cd();
    TDirectory* baseDir = fout->GetDirectory(subdirName);
    if (!baseDir) baseDir = fout->mkdir(subdirName);

    TDirectory* occDir = nullptr;
    if (savePlots) {
        occDir = baseDir->GetDirectory("OccupancyMaps");
        if (!occDir) occDir = baseDir->mkdir("OccupancyMaps");
        gSystem->Exec(Form("mkdir -p %s/%s", plotOccDir, subdirName));
    }

    // --- Loop TPC / Plane / Variable
    for (int tpc = 0; tpc < 4; ++tpc) {
        for (int plane = 0; plane < 3; ++plane) {
            for (int v = 0; v < 4; ++v) {

                TString histname = Form("h3D_%s_TPC%d_plane%d_%s", vars[v], tpc, plane, mode);
                TH3D* h3_orig = (TH3D*)fin->Get(histname);
                if (!h3_orig) {
                    std::cerr << "Missing histogram: " << histname << std::endl;
                    continue;
                }

                TH3D* h3_cut = applyExpCutToTH3(h3_orig, vars[v], plane, minEntries, cut_hit, cut_width);

                // Custom additional rebinning
                TH3D* h3 = rebinTH3_Custom(h3_cut, mode, tpc);
                if (!h3) continue;

                // --- Axis labels
                TString xtitle = "X [cm]";
                TString ytitle = (TString(mode) == "XTheta" || TString(mode) == "XTheta_c") ? "Theta_{XW} [deg]" : "dQ/dx";

                TString graphname = histname; graphname.Replace(0,4,"");

                TGraph2D* gr = new TGraph2D();
                gr->SetName(Form("%s_ITMgraph", graphname.Data()));
                gr->SetTitle(Form("ITM %s; %s; %s; ITM Value", vars[v], xtitle.Data(), ytitle.Data()));

                TGraph2D* gr_err = new TGraph2D();
                gr_err->SetName(Form("%s_ITMerrgraph", graphname.Data()));
                gr_err->SetTitle(Form("ITM err %s; %s; %s; ITM Error", vars[v], xtitle.Data(), ytitle.Data()));

                // --- Occupancy histogram
                int nbx = h3->GetNbinsX();
                int nby = h3->GetNbinsY();
                
                const TAxis* xax = h3->GetXaxis();
                const TAxis* yax = h3->GetYaxis();

                std::vector<double> xedges(nbx+1), yedges(nby+1);
                for (int i = 1; i <= nbx+1; ++i) xedges[i-1] = xax->GetBinLowEdge(i);
                for (int i = 1; i <= nby+1; ++i) yedges[i-1] = yax->GetBinLowEdge(i);

                TH2D* h2_entries = new TH2D(
                    Form("h2D_entries_%s_TPC%d_plane%d_%s", vars[v], tpc, plane, mode),
                    Form("Entries map %s (%s);%s;%s;Entries", vars[v], mode, xtitle.Data(), ytitle.Data()),
                    nbx, xedges.data(),
                    nby, yedges.data()
                );


                int point = 0;

                for (int ix = 1; ix <= nbx; ++ix) {
                    double xcenter = h3->GetXaxis()->GetBinCenter(ix);
                    for (int iy = 1; iy <= nby; ++iy) {
                        double ycenter = h3->GetYaxis()->GetBinCenter(iy);

                        TH1D* projZ = h3->ProjectionZ("_pz", ix, ix, iy, iy);
                        if (!projZ) continue;

                        // These are additional rebin factors which help the binned ITM in low statistcs bins... these are very empirical, sorry!
                        int rebin_factor = 1;
                        if (strcmp(vars[v], "hit_integral") == 0 && ycenter > 60) rebin_factor = 4;
                        if (strcmp(vars[v], "width") == 0) {
                            if (plane==0 && ycenter>39 && ycenter<50) rebin_factor=2;
                            else if (plane==0 && ycenter>=50 && ycenter<70) rebin_factor=2;
                            else if ((plane==0 && ycenter>=70) || (plane==1 && ycenter>=70)) rebin_factor=5;
                        }
                        projZ->Rebin(rebin_factor);

                        double nEntries = projZ->GetEntries();
                        h2_entries->SetBinContent(ix, iy, nEntries);

                        if (nEntries >= minEntries) {
                            Double_t result[2] = {0.,0.};
                            iterative_truncated_mean(projZ, sig_down, sig_up, tol, result);
                            gr->SetPoint(point++, xcenter, ycenter, result[0]);
                            gr_err->SetPoint(point++, xcenter, ycenter, result[1]);
                        }

                        delete projZ;
                    }
                }

                // --- Write outputs
                baseDir->cd(); gr->Write(); gr_err->Write();

                if (savePlots) {
                    occDir->cd();
                    h2_entries->Write();
                
                    TCanvas cocc("cocc","",800,600);
                    cocc.SetRightMargin(0.15); cocc.SetLeftMargin(0.12); cocc.SetBottomMargin(0.12);
                    cocc.SetLogz();
                    h2_entries->Draw("COLZ");
                    cocc.SaveAs(Form("%s/%s/occ_%s_TPC%d_plane%d_%s.png",
                                     plotOccDir, subdirName, vars[v], tpc, plane, mode));
                }

                delete gr; delete gr_err; delete h2_entries; delete h3;
            }
        }
    }

    fin->Close();
    delete cut_hit;
    for (int i=0;i<3;i++) delete cut_width[i];

    std::cout << "ITM TGraphs and (eventually) occupancy plots saved in: " << plotOccDir << std::endl;
}

int main(int argc, char** argv) {

    if (argc != 5) {
        std::cout << "\nUsage:\n"
                  << "./produceTGraphFromTH3  [Offbeam|Overlay]  input.root  mode  [true|false]\n\n";
        return 1;
    }

    std::string sampleType = argv[1];
    std::string infile     = argv[2];
    std::string mode       = argv[3];
    std::string saveStr    = argv[4];

    bool savePlots = (saveStr == "true" || saveStr == "1");

    std::string subdirName = (sampleType == "Offbeam") ? "Data" : "MC";

    std::string outfile = "ITM_" + sampleType + "_" + mode + ".root";
    std::string plotDir = "occupancy_" + sampleType + "_" + mode;

    std::cout << "=== Running ITM Analysis ===\n";
    std::cout << "Sample: " << sampleType << "\n";
    std::cout << "Input : " << infile << "\n";
    std::cout << "Mode  : " << mode << "\n";
    std::cout << "Save occupancy plots: " << (savePlots ? "YES" : "NO") << "\n\n";

    TFile* fout = new TFile(outfile.c_str(), "RECREATE");

    makeAll_ITM_TGraphs(
        infile.c_str(),
        fout,
        subdirName.c_str(),
        plotDir.c_str(),
        mode.c_str(),
        -2., 1.75,     // sig_down, sig_up
        1e-4,          // tol
        50,            // minEntries
        savePlots      // save png/pdf
    );

    fout->Close();
    delete fout;

    std::cout << "\nCompleted. Output file: " << outfile << "\n";
    return 0;
}
