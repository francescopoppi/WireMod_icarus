#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <iostream>
#include <vector>
#include "../CommonTools/TrackReader.h"
#include "../CommonTools/Corrections.h"
#include "../CommonTools/Utils.h"
#include "../CommonTools/binning.h"

void ntupleAnalyzerOffbeam() {

    std::map<int,float> lifetimeEE, lifetimeEW, lifetimeWE, lifetimeWW;
    loadLifetimeFile("ntupleProcessing/lifetime_run2_tpc.csv", lifetimeEE, lifetimeEW, lifetimeWE, lifetimeWW);

    Corrections corrections(lifetimeEE, lifetimeEW, lifetimeWE, lifetimeWW);

    std::vector<std::string> files = getFileListFromRunList("running.txt");

    struct HitVars {
        float x, y, z, integral, dqdx, width, pitch;
        float dirX, dirY, dirZ;
    };

    TFile* fout[4];
    TTree* tree[4][3];
    HitVars vars[4][3];

    for (int tpc=0; tpc<4; ++tpc) {
        fout[tpc] = new TFile(Form("hits_TPC%d.root", tpc), "RECREATE");
        for (int plane=0; plane<3; ++plane) {
            tree[tpc][plane] = new TTree(Form("Plane%d", plane),
                                         Form("TPC%d Plane%d hits", tpc, plane));
            tree[tpc][plane]->Branch("x",&vars[tpc][plane].x);
            tree[tpc][plane]->Branch("y",&vars[tpc][plane].y);
            tree[tpc][plane]->Branch("z",&vars[tpc][plane].z);
            tree[tpc][plane]->Branch("integral",&vars[tpc][plane].integral);
            tree[tpc][plane]->Branch("dqdx",&vars[tpc][plane].dqdx);
            tree[tpc][plane]->Branch("width",&vars[tpc][plane].width);
            tree[tpc][plane]->Branch("pitch",&vars[tpc][plane].pitch);
            tree[tpc][plane]->Branch("dirX",&vars[tpc][plane].dirX);
            tree[tpc][plane]->Branch("dirY",&vars[tpc][plane].dirY);
            tree[tpc][plane]->Branch("dirZ",&vars[tpc][plane].dirZ);
        }
    }

    for (const auto& f : files) {

        TFile* myFile = TFile::Open(f.c_str());
        if (!myFile || myFile->IsZombie()) {
            std::cerr << "Cannot open file " << f << "\n";
            continue;
        }

        TrackReader reader(lifetimeEE, lifetimeEW, lifetimeWE, lifetimeWW);
        auto WestTracks = reader.readFile(myFile, 1, false);
        auto EastTracks = reader.readFile(myFile, 0, false);

        myFile->Close();

        auto fillHit = [&](int tpc, int plane, const auto& hit){
            vars[tpc][plane].x        = hit.x;
            vars[tpc][plane].y        = hit.y;
            vars[tpc][plane].z        = hit.z;
            vars[tpc][plane].integral = hit.integral_c;
            vars[tpc][plane].dqdx     = hit.dqdx_c;
            vars[tpc][plane].width    = hit.width;
            vars[tpc][plane].pitch    = hit.pitch;
            vars[tpc][plane].dirX     = hit.dirx;
            vars[tpc][plane].dirY     = hit.diry;
            vars[tpc][plane].dirZ     = hit.dirz;
            tree[tpc][plane]->Fill();
        };

        for (const auto& trk : WestTracks) {
            for (int plane=0; plane<3; ++plane) {
                for (const auto& hit : trk.hits[plane]) {
                    fillHit(hit.tpc, plane, hit);
                }
            }
        }
        for (const auto& trk : EastTracks) {
            for (int plane=0; plane<3; ++plane) {
                for (const auto& hit : trk.hits[plane]) {
                    fillHit(hit.tpc, plane, hit);
                }
            }
        }
    }

    for (int tpc=0; tpc<4; ++tpc) {
        fout[tpc]->cd();
        for (int plane=0; plane<3; ++plane) tree[tpc][plane]->Write();
        fout[tpc]->Close();
    }

    std::cout << "Histo Filled and Root tree created." << std::endl;
}

int main() {
    ntupleAnalyzerOffbeam();
    return 0;
}
