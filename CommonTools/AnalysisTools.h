#ifndef ANALYSIS_TOOLS_H
#define ANALYSIS_TOOLS_H

#include <TH1.h>
#include <TH3.h>
#include <TF1.h>

void hist_mean_std_error(const TH1* h, Double_t* result);

void iterative_truncated_mean(TH1* h, Double_t sig_down, Double_t sig_up,
                              Double_t tol, Double_t* result,
                              Int_t iter = 0,
                              Double_t N0 = -1);

TF1* makeHitIntegralCut(const TString& mode);
TF1* makeWidthCut(int plane);

TH3D* applyExpCutToTH3(TH3D* h3,
                       const char* var,
                       int plane,
                       int minEntries,
                       TF1* cut_hit,
                       TF1* cut_width[3]);

TH3D* rebinTH3_Custom(TH3D* h3_orig, const char* mode, int tpc);
                              
#endif