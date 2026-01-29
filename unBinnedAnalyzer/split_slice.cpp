#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TString.h>

#include <iostream>
#include <vector>
#include <cstdlib>

#include "binning.h"

void splitSlice(const char* input_file,
                int tpc,
                const char* plane,
                int x_bin_index)
{
    initBins();

    std::vector<double> x_edges;
    std::vector<double> theta_edges;
    getAnalysisEdges_XThetaUnBinned(tpc, x_edges, theta_edges);

    double x_min = x_edges[x_bin_index];
    double x_max = x_edges[x_bin_index+1];

    TFile* inFile = TFile::Open(input_file, "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "File not found!" << std::endl;
        return;
    }

    TTree* tree = (TTree*) inFile->Get(plane);
    if (!tree) {
        std::cerr << "Plane not found!" << std::endl;
        return;
    }

    TString output_name = TString::Format("slice_x%d.root", x_bin_index);
    TFile* outFile = TFile::Open(output_name, "RECREATE");

    float f_dirX, f_pitch, f_x, f_integral, f_width;

    TTree* t_nominal = new TTree("nominal", "nominal bin");
    t_nominal->Branch("dirX", &f_dirX, "dirX/F");
    t_nominal->Branch("pitch", &f_pitch, "pitch/F");
    t_nominal->Branch("x", &f_x, "x/F");
    t_nominal->Branch("integral", &f_integral, "integral/F");
    t_nominal->Branch("width", &f_width, "width/F");

    float dirX, pitch, x, integral, width;
    tree->SetBranchAddress("dirX", &dirX);
    tree->SetBranchAddress("pitch", &pitch);
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("integral", &integral);
    tree->SetBranchAddress("width", &width);

    Long64_t nentries = tree->GetEntries();
    Long64_t print_every = nentries / 10;

    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        double thetaX = TMath::ATan(dirX * pitch / 0.3) * 180.0 / TMath::Pi();
        (void)thetaX;

        f_dirX = dirX;
        f_pitch = pitch;
        f_x = x;
        f_integral = integral;
        f_width = width;

        if (x >= x_min && x < x_max)
            t_nominal->Fill();

        if (print_every > 0 && i % print_every == 0) {
            int perc = int(100.0 * i / nentries);
            std::cout << "\rProcessing: " << perc << "% completed" << std::flush;
        }
    }

    std::cout << "\rProcessing: 100% completed" << std::endl;

    outFile->Write();
    outFile->Close();
    inFile->Close();

    std::cout << "File created: " << output_name.Data() << std::endl;
}

int main(int argc, char** argv)
{
    if (argc < 4) {
        std::cout << "Usage:\n"
                  << "./produceSlice input.root TPC plane [x_bin_index]\n";
        return 1;
    }

    const char* input_file = argv[1];
    int tpc = std::atoi(argv[2]);
    const char* plane = argv[3];
    int x_bin_index = (argc > 4) ? std::atoi(argv[4]) : 0;

    splitSlice(input_file, tpc, plane, x_bin_index);
    return 0;
}