#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TString.h>

#include <iostream>
#include <vector>
#include <cstdlib>

#include "binning.h"

void splitSliceYZ_4(const char* input_file,
                  int tpc,
                  const char* plane,
                  int z_bin_index)
{
    initBins();

    double z_min = binsZ_2[z_bin_index];
    double z_max = binsZ_2[z_bin_index + 1];

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

    TString output_name = TString::Format("slice_z_rebin4_%d.root", z_bin_index);
    TFile* outFile = TFile::Open(output_name, "RECREATE");

    float f_dirX, f_pitch, f_y, f_z, f_integral, f_width;

    TTree* t_nominal = new TTree("nominal", "Z slice");

    t_nominal->Branch("dirX",     &f_dirX,     "dirX/F");
    t_nominal->Branch("pitch",    &f_pitch,   "pitch/F");
    t_nominal->Branch("y",        &f_y,        "y/F");
    t_nominal->Branch("z",        &f_z,        "z/F");
    t_nominal->Branch("integral", &f_integral, "integral/F");
    t_nominal->Branch("width",    &f_width,    "width/F");

    float dirX, pitch, x, y, z, integral, width;

    tree->SetBranchAddress("dirX", &dirX);
    tree->SetBranchAddress("pitch", &pitch);
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("z", &z);
    tree->SetBranchAddress("integral", &integral);
    tree->SetBranchAddress("width", &width);

    Long64_t nentries = tree->GetEntries();
    Long64_t print_every = nentries / 10;

    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        f_dirX     = dirX;
        f_pitch    = pitch;
        f_y        = y;
        f_z        = z;
        f_integral = integral;
        f_width    = width;

        if (z >= z_min && z < z_max)
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
                  << "./produceSliceYZ_rebin4 input.root TPC plane [z_bin_index]\n";
        return 1;
    }

    const char* input_file = argv[1];
    int tpc = std::atoi(argv[2]);
    const char* plane = argv[3];
    int z_bin_index = (argc > 4) ? std::atoi(argv[4]) : 0;

    splitSliceYZ_4(input_file, tpc, plane, z_bin_index);
    return 0;
}