#include "AnalysisTools.h"
#include <TH1.h>

void hist_mean_std_error(const TH1* h, Double_t* result) {

    if (h->GetEntries() == 0) {
        result[0] = 0.;
        result[1] = 0.;
        return; 
    }

    Double_t mean  = h->GetMean();
    Double_t sigma = h->GetRMS();
    Double_t N     = h->Integral();
    Double_t err   = sigma / sqrt(N);

    result[0] = mean;
    result[1] = err;
}

void iterative_truncated_mean(TH1* h, Double_t sig_down, Double_t sig_up,
                              Double_t tol, Double_t* result,
                              Int_t iter, Double_t N0) {

    if (sig_down > sig_up) std::swap(sig_down, sig_up);
    if (N0 < 0) N0 = h->Integral();

    Double_t mean = h->GetMean();
    Double_t sd   = h->GetRMS();

    if (TMath::Abs(result[0] - mean) < tol) {
        hist_mean_std_error(h, result);
        return;
    }

    result[0] = mean;
    result[1] = sd;

    Double_t probs[1] = {0.5};
    Double_t median[1];
    h->GetQuantiles(1, median, probs);

    TH1* hnew = (TH1*)h->Clone("hclone");
    hnew->Reset();

    for (int i = 1; i <= h->GetNbinsX(); i++) {
        Double_t xlow    = h->GetBinLowEdge(i);
        Double_t xupp    = xlow + h->GetBinWidth(i);

        if (xupp < median[0] + sig_down * sd) continue;
        if (xlow > median[0] + sig_up   * sd) continue;

        hnew->SetBinContent(i, h->GetBinContent(i));
        hnew->SetBinError(i,   h->GetBinError(i));
    }

    iterative_truncated_mean(hnew, sig_down, sig_up, tol,
                             result, iter + 1, N0);
    delete hnew;
}

TH3D* applyExpCutToTH3(
    TH3D* h3,
    const char* var,
    int plane,
    int minEntries,
    TF1* cut_hit,
    TF1* cut_width[3]
) {
    if (!h3) return nullptr;

    TH3D* h3_cut = (TH3D*)h3->Clone(Form("%s_expCut", h3->GetName()));
    h3_cut->SetDirectory(0);

    int nbx = h3->GetNbinsX();
    int nby = h3->GetNbinsY();
    int nbz = h3->GetNbinsZ();

    for (int ix = 1; ix <= nbx; ix++) {
        for (int iy = 1; iy <= nby; iy++) {

            TH1D* projZ = h3->ProjectionZ("_pz_tmp", ix, ix, iy, iy);
            if (!projZ) continue;

            double nEntries = projZ->GetEntries();
            if (nEntries < minEntries) {
                for (int iz = 1; iz <= nbz; iz++) {
                    h3_cut->SetBinContent(ix, iy, iz, 0.);
                    h3_cut->SetBinError(ix, iy, iz, 0.);
                }
                delete projZ;
                continue;
            }

            double ycenter = h3->GetYaxis()->GetBinCenter(iy);

            double cut = 0.;
            if (strcmp(var, "hit_integral") == 0) {
                cut = cut_hit->Eval(ycenter);
            }
            else if (strcmp(var, "width") == 0) {
                cut = cut_width[plane]->Eval(ycenter);
            }

            for (int iz = 1; iz <= nbz; iz++) {
                double zcenter = h3->GetZaxis()->GetBinCenter(iz);
                if (zcenter < cut) {
                    h3_cut->SetBinContent(ix, iy, iz, 0.);
                    h3_cut->SetBinError(ix, iy, iz, 0.);
                }
            }

            delete projZ;
        }
    }

    return h3_cut;
}

TH3D* rebinTH3_Custom(TH3D* h3_orig, const char* mode, int tpc) {
    if (!h3_orig) return nullptr;

    std::vector<double> xedges, yedges, zedges;
    if (TString(mode)=="XTheta" || TString(mode)=="XTheta_c")
        getAnalysisEdges_XTheta(tpc, xedges, yedges, zedges, h3_orig);
    else {
        std::cerr << "Unsupported mode: " << mode << std::endl;
        return nullptr;
    }

    TH3D* h3_new = new TH3D(Form("%s_custom", h3_orig->GetName()),
                            h3_orig->GetTitle(),
                            xedges.size()-1, xedges.data(),
                            yedges.size()-1, yedges.data(),
                            zedges.size()-1, zedges.data());

    // copia contenuti
    int nXorig = h3_orig->GetNbinsX();
    int nYorig = h3_orig->GetNbinsY();
    int nZorig = h3_orig->GetNbinsZ();

    for (int ix=1; ix<=nXorig; ix++)
        for (int iy=1; iy<=nYorig; iy++)
            for (int iz=1; iz<=nZorig; iz++) {
                double content = h3_orig->GetBinContent(ix,iy,iz);
                double error   = h3_orig->GetBinError(ix,iy,iz);
                if (content==0 && error==0) continue;

                int newbin = h3_new->FindBin(h3_orig->GetXaxis()->GetBinCenter(ix),
                                             h3_orig->GetYaxis()->GetBinCenter(iy),
                                             h3_orig->GetZaxis()->GetBinCenter(iz));
                double oldcontent = h3_new->GetBinContent(newbin);
                double olderror   = h3_new->GetBinError(newbin);
                h3_new->SetBinContent(newbin, oldcontent + content);
                h3_new->SetBinError(newbin, std::sqrt(olderror*olderror + error*error));
            }

    return h3_new;
}


// ---------------------------

TF1* makeHitIntegralCut(const TString& mode) {
    TF1* f = new TF1("cut_hit", "[0]*exp([1]*x)+[2]", 0, 100);
    f->SetParameters(0.631142, 0.0828856, -8.8824);

    if (mode == "XTheta_c")
        f->SetParameters(1.35696, 0.0786976, -24.1874);

    return f;
}

TF1* makeWidthCut(int plane) {

    if (plane == 0) {
        TF1* f = new TF1("cut_width0", "[0]*exp([1]*x)+[2]", 0, 100);
        f->SetParameters(0.108951, 0.0574775, 0.850119);
        return f;
    }
    if (plane == 1) {
        TF1* f = new TF1("cut_width1", "[0]*exp([1]*x)+[2]", 0, 100);
        f->SetParameters(0.0444325, 0.0709458, 0.455733);
        return f;
    }

    TF1* f = new TF1("cut_width2", "[0]*exp([1]*x)+[2]", 0, 100);
    f->SetParameters(0.0833033, 0.0544162, 0.732181);
    return f;
}
