#include "TrackReader.h"
#include "DataStructures.h"

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TMath.h>

#include <cmath>
#include <iostream>

TrackReader::TrackReader(std::map<int,float>& EE,
                         std::map<int,float>& EW,
                         std::map<int,float>& WE,
                         std::map<int,float>& WW)
    : lifetimeEE(EE), lifetimeEW(EW), lifetimeWE(WE), lifetimeWW(WW)
{}

float TrackReader::getLifetime(int run, int cryo, int tpc)
{
    if (cryo == 0 && (tpc == 0 || tpc == 1)) return lifetimeEE[run];
    if (cryo == 0 && (tpc == 2 || tpc == 3)) return lifetimeEW[run];
    if (cryo == 1 && (tpc == 0 || tpc == 1)) return lifetimeWE[run];
    if (cryo == 1 && (tpc == 2 || tpc == 3)) return lifetimeWW[run];

    return 1e9; 
}


std::vector<Track> TrackReader::readFile(TFile* file, int cryo, bool applyTruthCuts)
{
    std::vector<Track> tracks;

    std::string treename = (cryo==0 ? "caloskimE/TrackCaloSkim"
                                   : "caloskimW/TrackCaloSkim");

    TTreeReader reader(treename.c_str(), file);

    TTreeReaderValue<int>   runTPC(reader,   "meta.run");
    TTreeReaderValue<int>   eventTPC(reader, "meta.evt");
    TTreeReaderValue<int>   selected(reader, "selected");
    TTreeReaderValue<int>   whicht0(reader,  "whicht0");
    TTreeReaderValue<int>   truthG4(reader,  "truth.p.G4ID");

    TTreeReaderValue<float> t0PFP(reader,    "t0PFP");
    TTreeReaderValue<float> t0CRTHit(reader, "t0CRTHit");
    TTreeReaderValue<float> length(reader,   "length");

    TTreeReaderValue<float> startx(reader,"start.x");
    TTreeReaderValue<float> starty(reader,"start.y");
    TTreeReaderValue<float> startz(reader,"start.z");
    TTreeReaderValue<float> endx(reader,"end.x");
    TTreeReaderValue<float> endy(reader,"end.y");
    TTreeReaderValue<float> endz(reader,"end.z");

    TTreeReaderValue<float> dirx(reader,"dir.x");
    TTreeReaderValue<float> diry(reader,"dir.y");
    TTreeReaderValue<float> dirz(reader,"dir.z");

    TTreeReaderValue<float> PCAdirx(reader,"PCAdir.x");
    TTreeReaderValue<float> PCAdiry(reader,"PCAdir.y");
    TTreeReaderValue<float> PCAdirz(reader,"PCAdir.z");

    TTreeReaderArray<UShort_t> h_tpc[3] = {
        {reader,"hits0.h.tpc"}, {reader,"hits1.h.tpc"}, {reader,"hits2.h.tpc"} };

    TTreeReaderArray<UShort_t> h_mult[3] = {
        {reader,"hits0.h.mult"}, {reader,"hits1.h.mult"}, {reader,"hits2.h.mult"} };

    TTreeReaderArray<float> h_width[3] = {
        {reader,"hits0.h.width"}, {reader,"hits1.h.width"}, {reader,"hits2.h.width"} };

    TTreeReaderArray<float> h_goodness[3] = {
        {reader,"hits0.h.goodness"}, {reader,"hits1.h.goodness"}, {reader,"hits2.h.goodness"} };

    TTreeReaderArray<float> h_tic[3] = {
        {reader,"hits0.h.time"}, {reader,"hits1.h.time"}, {reader,"hits2.h.time"} };

    TTreeReaderArray<float> h_pitch[3] = {
        {reader,"hits0.pitch"}, {reader,"hits1.pitch"}, {reader,"hits2.pitch"} };

    TTreeReaderArray<float> h_dqdx[3] = {
        {reader,"hits0.dqdx"}, {reader,"hits1.dqdx"}, {reader,"hits2.dqdx"} };

    TTreeReaderArray<float> h_rr[3] = {
        {reader,"hits0.rr"}, {reader,"hits1.rr"}, {reader,"hits2.rr"} };

    TTreeReaderArray<float> h_x[3] = {
        {reader,"hits0.tp.x"}, {reader,"hits1.tp.x"}, {reader,"hits2.tp.x"} };

    TTreeReaderArray<float> h_y[3] = {
        {reader,"hits0.tp.y"}, {reader,"hits1.tp.y"}, {reader,"hits2.tp.y"} };

    TTreeReaderArray<float> h_z[3] = {
        {reader,"hits0.tp.z"}, {reader,"hits1.tp.z"}, {reader,"hits2.tp.z"} };

    TTreeReaderArray<float> h_dirx[3] = {
        {reader,"hits0.dir.x"}, {reader,"hits1.dir.x"}, {reader,"hits2.dir.x"} };

    TTreeReaderArray<float> h_diry[3] = {
        {reader,"hits0.dir.y"}, {reader,"hits1.dir.y"}, {reader,"hits2.dir.y"} };

    TTreeReaderArray<float> h_dirz[3] = {
        {reader,"hits0.dir.z"}, {reader,"hits1.dir.z"}, {reader,"hits2.dir.z"} };

    TTreeReaderArray<float> h_integral[3] = {
        {reader,"hits0.h.integral"}, {reader,"hits1.h.integral"}, {reader,"hits2.h.integral"} };

    TTreeReaderArray<float> h_sumadc[3] = {
        {reader,"hits0.h.sumadc"}, {reader,"hits1.h.sumadc"}, {reader,"hits2.h.sumadc"} };

    TTreeReaderArray<bool> h_ontraj[3] = {
        {reader,"hits0.ontraj"}, {reader,"hits1.ontraj"}, {reader,"hits2.ontraj"} };

    TTreeReaderArray<float> h_truth_ne[3] = {
        {reader,"hits0.h.truth.nelec"}, {reader,"hits1.h.truth.nelec"}, {reader,"hits2.h.truth.nelec"} };


    while (reader.Next()) {

        if (lifetimeEE.find(*runTPC) == lifetimeEE.end()) {
            std::cout << "Run " << *runTPC << " not in lifetime map\n";
            continue;
        }

        double t0 = *t0PFP;
        if (std::isnan(t0)) t0 = *t0CRTHit;
        if (std::isnan(t0)) continue;

        Track trk;
        trk.run = *runTPC;
        trk.event = *eventTPC;
        trk.cryo = cryo;
        trk.selected = *selected;
        trk.whicht0 = *whicht0;
        trk.G4ID = *truthG4;
        trk.t0PFP = *t0PFP;
        trk.t0CRTHit = *t0CRTHit;
        trk.length = *length;

        trk.startx=*startx; trk.starty=*starty; trk.startz=*startz;
        trk.endx=*endx; trk.endy=*endy; trk.endz=*endz;
        trk.dirx=*dirx; trk.diry=*diry; trk.dirz=*dirz;
        trk.PCAdirx=*PCAdirx; trk.PCAdiry=*PCAdiry; trk.PCAdirz=*PCAdirz;

        if (applyTruthCuts && trk.G4ID == -1) continue;

        for (int plane = 0; plane < 3; plane++) {

            for (unsigned i = 0; i < h_tpc[plane].GetSize(); i++) {

                if (!h_ontraj[plane][i]) continue; // remove hits not on trajectory (delta rays)
                if (std::isnan(h_x[plane][i]) || std::isnan(h_y[plane][i]) || std::isnan(h_z[plane][i])) continue; // remove bad hits

                if (applyTruthCuts && h_truth_ne[plane][i] == 0) continue; // clean mc hits, hopefully

                Hit hit;
                hit.x = h_x[plane][i];
                hit.y = h_y[plane][i];
                hit.z = h_z[plane][i];

                hit.width = h_width[plane][i];
                hit.goodness = h_goodness[plane][i];
                hit.pitch = h_pitch[plane][i];
                hit.rr = h_rr[plane][i];
                hit.dirx = h_dirx[plane][i];
                hit.diry = h_diry[plane][i];
                hit.dirz = h_dirz[plane][i];
                hit.tic  = h_tic[plane][i];
                hit.mult = h_mult[plane][i];
                hit.ontraj = h_ontraj[plane][i];

                int rawTPC = h_tpc[plane][i];
                if (cryo==0 && (rawTPC==0||rawTPC==1)) hit.tpc=0;
                else if (cryo==0 && (rawTPC==2||rawTPC==3)) hit.tpc=1;
                else if (cryo==1 && (rawTPC==0||rawTPC==1)) hit.tpc=2;
                else if (cryo==1 && (rawTPC==2||rawTPC==3)) hit.tpc=3;

                float lifetime = getLifetime(trk.run, cryo, rawTPC);

                double driftT = (hit.tic * 0.4 - 340) - t0 / 1.e3;
                double corr = std::exp(driftT / lifetime);

                hit.integral   = h_integral[plane][i];
                hit.integral_c = h_integral[plane][i] * corr;
                hit.dqdx       = h_dqdx[plane][i];
                hit.dqdx_c     = h_dqdx[plane][i] * corr;

                trk.hits[plane].push_back(hit);
            }
        }

        tracks.push_back(trk);
    }

    return tracks;
}
