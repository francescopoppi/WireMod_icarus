#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <iostream>
#include <vector>

#include "../CommonTools/TrackReader.h"
#include "../CommonTools/Corrections.h"
#include "../CommonTools/Utils.h"
#include "../CommonTools/TrackCuts.h"
#include "../CommonTools/binning.h"
#include "../CommonTools/HistogramManager.h"

void runUnbinned(const std::vector<std::string>& files,
                 Corrections& corrections)
{
    struct HitVars {
        float x,y,z,integral,dqdx,width,pitch;
        float dirX,dirY,dirZ;
    };

    TFile* fout[4];
    TTree* tree[4][3];
    HitVars vars[4][3];

    for (int tpc=0;tpc<4;++tpc) {
        fout[tpc] = new TFile(Form("hits_TPC%d.root",tpc),"RECREATE");
        for (int plane=0;plane<3;++plane) {
            tree[tpc][plane] = new TTree(Form("Plane%d",plane),
                                         Form("TPC%d Plane%d hits",tpc,plane));
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

    TrackReader reader(corrections);

    for (const auto& f : files) {

        TFile* file = TFile::Open(f.c_str());
        if (!file || file->IsZombie()) continue;

        auto West = reader.readFile(file,1,false);
        auto East = reader.readFile(file,0,false);
        file->Close();

        auto fillHit = [&](int tpc,int plane,const Hit& hit){
            vars[tpc][plane].x = hit.x;
            vars[tpc][plane].y = hit.y;
            vars[tpc][plane].z = hit.z;
            vars[tpc][plane].integral = hit.integral_c;
            vars[tpc][plane].dqdx = hit.dqdx_c;
            vars[tpc][plane].width = hit.width;
            vars[tpc][plane].pitch = hit.pitch;
            vars[tpc][plane].dirX = hit.dirx;
            vars[tpc][plane].dirY = hit.diry;
            vars[tpc][plane].dirZ = hit.dirz;
            tree[tpc][plane]->Fill();
        };

        for (const auto& trk : West)
            if (Cuts::trackLength(trk))
                for (int p=0;p<3;++p)
                    for (const auto& hit : trk.hits[p])
                        if (Cuts::multiplicity1(hit) && !Cuts::danglingCable(hit))
                            fillHit(hit.tpc,p,hit);

        for (const auto& trk : East)
            if (Cuts::trackLength(trk))
                for (int p=0;p<3;++p)
                    for (const auto& hit : trk.hits[p])
                        if (Cuts::multiplicity1(hit) && !Cuts::danglingCable(hit))
                            fillHit(hit.tpc,p,hit);
    }

    for (int t=0;t<4;++t) {
        fout[t]->cd();
        for (int p=0;p<3;++p) tree[t][p]->Write();
        fout[t]->Close();
    }
}

void runBinned(const std::vector<std::string>& files,
               Corrections& corrections)
{
    HistogramManager hm;
    TrackReader reader(corrections);

    for (const auto& f : files) {

        TFile* file = TFile::Open(f.c_str());
        if (!file || file->IsZombie()) continue;

        auto West = reader.readFile(file,1,false);
        auto East = reader.readFile(file,0,false);
        file->Close();

        auto processTrack = [&](const Track& trk){
            if (!Cuts::trackLength(trk)) return;

            for (int plane=0;plane<3;++plane) {
                for (const auto& hit : trk.hits[plane]) {

                    if (!Cuts::multiplicity1(hit)) continue;
                    if (Cuts::danglingCable(hit)) continue;

                    double thetaXW =
                        std::atan(hit.dirx * hit.pitch / 0.3) * 180.0 / M_PI;

                    hm.fill(hit, plane, thetaXW);
                }
            }
        };

        for (const auto& trk : West) processTrack(trk);
        for (const auto& trk : East) processTrack(trk);
    }

    TFile fout("binnedHistograms.root","RECREATE");
    hm.write(fout);
    fout.Close();
}


int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cout << "Usage: ./ntupleAnalyzerOffbeam [binned|unbinned]\n";
        return 1;
    }

    std::string mode = argv[1];

    std::map<int,float> EE,EW,WE,WW;
    loadLifetimeFile("ntupleProcessing/lifetime_run2_tpc.csv",EE,EW,WE,WW);
    Corrections corrections(EE,EW,WE,WW);

    auto files = getFileListFromRunList("running.txt");

    if (mode == "unbinned") {
        std::cout << "Running UNBINNED mode\n";
        runUnbinned(files, corrections);
    }
    else if (mode == "binned") {
        std::cout << "Running BINNED mode\n";
        runBinned(files, corrections);
    }
    else {
        std::cout << "Unknown mode: " << mode << "\n";
        return 1;
    }

    std::cout << "Done.\n";
    return 0;
}
