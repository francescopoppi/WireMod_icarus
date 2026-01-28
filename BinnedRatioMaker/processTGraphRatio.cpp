#include <TH2D.h>
#include <TGraph2D.h>
#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <string>

TH2D* graphToTH2D_usingTemplate(TGraph2D* g, TH2D* hTemplate) {
    if (!g || !hTemplate) return nullptr;

    TH2D* h = (TH2D*)hTemplate->Clone(Form("%s_TH2D", g->GetName()));
    h->Reset();
    h->SetDirectory(nullptr); 
    int n = g->GetN();
    double* x = g->GetX();
    double* y = g->GetY();
    double* z = g->GetZ();

    for (int i = 0; i < n; ++i) {
        int bin = h->FindBin(x[i], y[i]);
        h->SetBinContent(bin, z[i]);
    }
    return h;
}

void compute_ITM_ratio(
    const char* dataFile,
    const char* mcFile,
    const char* mode,
    const char* outFile
) {
    gSystem->Exec("mkdir -p ratioResults");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);

    TFile* fData = TFile::Open(dataFile, "READ");
    TFile* fMC   = TFile::Open(mcFile, "READ");
    if (!fData) { std::cerr << "Data file not found." << std::endl; return; }
    if (!fMC)   { std::cerr << "MC file not found." << std::endl; return; }

    TDirectory* dData = fData->GetDirectory("Data");
    TDirectory* dMC   = fMC->GetDirectory("MC");
    if (!dData) { std::cerr << "Data directory not found!" << std::endl; return; }
    if (!dMC)   { std::cerr << "MC directory not found!" << std::endl; return; }

    TFile* fOut = TFile::Open(outFile, "RECREATE");

    TH2D* hTemplate[4] = {nullptr, nullptr, nullptr, nullptr};
    for (int tpc = 0; tpc < 4; ++tpc) {
        TString name = Form("hEdges_XTheta_TPC%d", tpc);
        hTemplate[tpc] = (TH2D*)fData->Get(name);
        if (!hTemplate[tpc]) std::cerr << "Edges Template for TPC " << tpc << " not found!" << std::endl;
    }

    TIter next(dData->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        if (strcmp(key->GetClassName(), "TGraph2D") != 0) continue;

        TString name = key->GetName();
        if (!name.Contains("_ITMgraph")) continue;

        TString nameErr = name;
        nameErr.ReplaceAll("_ITMgraph", "_ITMerrgraph");

        TGraph2D* gData    = (TGraph2D*)dData->Get(name);
        TGraph2D* gDataErr = (TGraph2D*)dData->Get(nameErr);
        TGraph2D* gMC      = (TGraph2D*)dMC->Get(name);
        TGraph2D* gMCErr   = (TGraph2D*)dMC->Get(nameErr);

        if (!gData || !gDataErr || !gMC || !gMCErr) {
            std::cerr << "TGraph2D missing: " << name << std::endl;
            continue;
        }

        int tpc = -1;
        if      (name.Contains("TPC0")) tpc = 0;
        else if (name.Contains("TPC1")) tpc = 1;
        else if (name.Contains("TPC2")) tpc = 2;
        else if (name.Contains("TPC3")) tpc = 3;

        if (tpc < 0 || !hTemplate[tpc]) continue;

        TH2D* hData    = graphToTH2D_usingTemplate(gData, hTemplate[tpc]);
        TH2D* hDataErr = graphToTH2D_usingTemplate(gDataErr, hTemplate[tpc]);
        TH2D* hMC      = graphToTH2D_usingTemplate(gMC, hTemplate[tpc]);
        TH2D* hMCErr   = graphToTH2D_usingTemplate(gMCErr, hTemplate[tpc]);

        if (!hData || !hDataErr || !hMC || !hMCErr) continue;

        TH2D* hRatio    = (TH2D*)hData->Clone(Form("%s_ratio_TH2D", name.Data()));
        TH2D* hRatioErr = (TH2D*)hData->Clone(Form("%s_ratioErr_TH2D", name.Data()));
        hRatio->SetDirectory(nullptr);
        hRatioErr->SetDirectory(nullptr);
        hRatio->Reset();
        hRatioErr->Reset();

        for (int ix = 1; ix <= hRatio->GetNbinsX(); ++ix) {
            for (int iy = 1; iy <= hRatio->GetNbinsY(); ++iy) {
                double d  = hData->GetBinContent(ix, iy);
                double de = hDataErr->GetBinContent(ix, iy);
                double m  = hMC->GetBinContent(ix, iy);
                double me = hMCErr->GetBinContent(ix, iy);
                if (d <= 0 || m <= 0) continue;
                double r  = d / m;
                double re = r * std::sqrt( (de/d)*(de/d) + (me/m)*(me/m) );
                hRatio->SetBinContent(ix, iy, r);
                hRatioErr->SetBinContent(ix, iy, re);
            }
        }

        TGraph2D* gRatio = new TGraph2D();
        TGraph2D* gRatioErr = new TGraph2D();
        int p = 0;
        for (int ix = 1; ix <= hRatio->GetNbinsX(); ++ix) {
            double xc = hRatio->GetXaxis()->GetBinCenter(ix);
            for (int iy = 1; iy <= hRatio->GetNbinsY(); ++iy) {
                double yc = hRatio->GetYaxis()->GetBinCenter(iy);
                double z  = hRatio->GetBinContent(ix, iy);
                double ze = hRatioErr->GetBinContent(ix, iy);
                if (z == 0) continue;
                gRatio->SetPoint(p, xc, yc, z);
                gRatioErr->SetPoint(p, xc, yc, ze);
                ++p;
            }
        }

        TString simpleName = name;
        simpleName.ReplaceAll("h3D_", "");
        Ssiz_t pos = simpleName.Index("_ITMgraph");
        if (pos != kNPOS) simpleName.Remove(pos);
        if (!simpleName.EndsWith("_ratio")) simpleName += "_ratio";

        gRatio->SetName(simpleName);
        gRatio->SetTitle(simpleName);

        TString simpleNameErr = simpleName + "_err";
        gRatioErr->SetName(simpleNameErr);
        gRatioErr->SetTitle(simpleNameErr);

        hRatio->SetName(simpleName + "_TH2D");
        hRatioErr->SetName(simpleNameErr + "_TH2D");
        hRatio->SetTitle(simpleName + "_TH2D");
        hRatioErr->SetTitle(simpleNameErr + "_TH2D");

        hRatio->GetXaxis()->SetTitle("X [cm]");
        hRatio->GetYaxis()->SetTitle("#theta_{XW} [deg]");
        hRatio->GetZaxis()->SetTitle("Data / MC");
        
        hRatioErr->GetXaxis()->SetTitle("X [cm]");
        hRatioErr->GetYaxis()->SetTitle("#theta_{XW} [deg]");
        hRatioErr->GetZaxis()->SetTitle("Data / MC");

        gRatio->GetXaxis()->SetTitle("X [cm]");
        gRatio->GetYaxis()->SetTitle("#theta_{XW} [deg]");
        gRatio->GetZaxis()->SetTitle("Data / MC");

        gRatioErr->GetXaxis()->SetTitle("X [cm]");
        gRatioErr->GetYaxis()->SetTitle("#theta_{XW} [deg]");
        gRatioErr->GetZaxis()->SetTitle("Data / MC");

        TCanvas* c = new TCanvas("c","c",1000,700);
        c->SetRightMargin(0.18); 
        hRatio->Draw("COLZ");
        TString outPng = Form("ratioResults/ratio2D_%s.png", simpleName.Data());
        c->SaveAs(outPng);
        delete c; 

        fOut->WriteTObject(hRatio);
        fOut->WriteTObject(hRatioErr);
        fOut->WriteTObject(gRatio);
        fOut->WriteTObject(gRatioErr);

        delete hRatio;
        delete hRatioErr;
        delete gRatio;
        delete gRatioErr;
    }

    fOut->Close();
    fData->Close();
    fMC->Close();

    std::cout << "Data/MC ratio saved: " << outFile << std::endl;
}

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cout << "Usage: ./produceRatio inputData_TGraph2D.root inputMC_TGraph2D.root mode output.root\n";
        return 1;
    }

    const char* dataFile = argv[1];
    const char* mcFile   = argv[2];
    const char* mode     = argv[3];
    const char* outFile  = argv[4];

    compute_ITM_ratio(dataFile, mcFile, mode, outFile);
    return 0;
}
