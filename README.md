# Getting WireMod Data/MC ratios in ICARUS

## General Comments
This repository was made for:
- reproducibility (everyone with access to Data and MC can test it);
- sharing the code within the ICARUS/SBN community;
- have common executable for Data and MC which can be cloned and used on the grid.

The requirements to run these utilities are simple:
1. a Root version (lower limit to be defined, but was tested with root 6.36.04)
2. eventually python (python requirements tbd)

These utilities run both in SL7 and AL9, no need to setup icaruscode/sbncode/sbndcode.

**If you just want to obtain the ratios and you want to skip the tedious process of producing the binned and unbinned dataset, you can use those that I have already processed. They can be found on FNAL persistent storage here:
_/pnfs/icarus/persistent/users/fpoppi/WireMod/BinnedFiles/_
You will find two files:
_OffbeamRun2_XTheta_ltc_fv.root_  _OverlayRun2_XTheta_ltc_fv.root_
These two files contain the TH3D of offbeam and overlay, they have both lifetime corrected and uncorrected integrals.
Some notation I used: ltc means that they included the lifetime corrections, fv means that they included the fiducialization around the dangling cable.
These two files can be fed to **produceTGraphFromTH3** skipping **ntupleAnalyzer**.

For what concern the unbinned datasets they can be found here:
_/pnfs/icarus/persistent/users/fpoppi/WireMod/UnBinnedFiles/_


This repo is not frozen! Expect changes and updates.

## Make
All the executables can be obtained by simply running `make`. Please report if the compiler report errors.
Three executables will be produced:
1. ntupleAnalyzer;
2. produceTGraphFromTH3;
3. produceRatio;

## ntupleAnalyzer

This executable analyzes a list of Calibration Ntuples and produces reduced trees or TH3s.

**Usage: ./ntupleAnalyzer [Offbeam|Overlay] [binned|unbinned|both] [fileListPath]**

The first argument is needed to specify the input dataset, currently two options are available: Offbeam (for Offbeam Data) and Overlay (for Overlay + simulated cosmics); in the future full-MC would be additional options. This argument is needed to select only the "simulated" hits and tracks on top of the data overlay. Regular MC can still be processed in the sameway by treating it as regular Offbeam (_I know, counterintuitive, but since it is full MC there is no need of looking for truth level matches_). LifeTime corrections are applied from TrackReader to hit integral and hit dqdx. Common selection cuts are applied and defined in TrackCuts.h (fiducial volume: hanging cable; hit multiplicity = 1; Track Length > 50 cm).
Binned/Unbinned/both these are the main options.
1. Binned: this analyzes a list of files and fills TH3 ntuples with fixed bins in X/Y/Z. What is X/Y/Z can be implemented. The HistogramManager (CommonTools/HistogramManager) manages and books the relevant histograms. Currently Two "modes" are implemented: XTheta (hitX vs hit theta vs (hit width, hit integral, hit goodness of fit, hit pitch)) and YZ (hitY vs hit Z vs (hit width, hit integral, hit goodness of fit, hit pitch)). Binnings can be tuned in the corresponding header file. The TH3s are stored in the output file: binnedHistograms_Offbeam/Overlay.root .
2. Unbinned: this option does not produce a root file with TH3s but reduced trees for all the different tpcs and wire planes. The reduced trees contain the following variables for each hit: x, y, z, dqdx, width, pitch, dirX, dirY and dirZ. Three root files are produced one per TPC in Hits_TPC0/1/2/3_OffBeam/Overlay.root, each of these root files contain 3 different Trees, one per each Plane.
3. Both: does the 2 things at the same time in a currently not efficient way.
The last argument wants an input.list of files to be analyzed, I added examples in the PathFolders.

Important note: this script is meant to be run on a small list of files. My implementation is done in a way that I split the full list of files in subset of 30-40 files, each job brings on the node a list of files which is then fed to this executable. All the jobs output are then supposed to me merged toghether with hadd. It can be run on a big list nonetheless with some caveat: it takes more time, it takes more RAM, this is expecially true if ran in unbinned mode where for the full Run2 Offbeam the Hits_TPC.root files are 150 GBs, while the Overlay ones are some units of GBs each. Luckily for you, I already have those available, I just need to determine a final path for those. 

## produceTGraphFromTH3

This executable takes as input the TH3s produced by the ntuple analyzers, applies some selection cuts (removes "bad hits" based on a angle dependent cut), runs ITM calculation on the pruned histograms and obtain the TGraphs to be used for the ratio evaluation.

**Usage: ./produceTGraphFromTH3  [Offbeam|Overlay]  input.root  mode  output.root  [true|false]**

The first argument is needed to specify the input dataset, as in the previous case two options are currently available: Offbeam (for Offbeam Data) and Overlay (for Overlay + simulated cosmics).
The second argument is supposed to be the merge of the ntupleAnalyzer in binned mode production.
The third argument currently offers the possibility to process only XTheta and XTheta_c histograms, YZ mode is soon to be implemented. XTheta uses hit integral not crrected for lifetime, while XTheta_c uses hit integral corrected for lifetime.
The fourth argument allows to specifiy the output file name.
The fifth argument is a boolean to store or not canvas with the occupancy (number of hits) of X Theta bins.

The output rootfiles contains a folder called Data/MC for Offbeam/Overlay mode, the TGraphs are stored inside this folder. Two different TGraph2D are stored, one for the actual ITM, and one for the error of each bin. Out of the Data/MC folder there are 4 empty TH2D which are used as binning template from the ratio maker.

## produceRatio

This executable is the final step in the ratio evaluation. It simply generates the TGraph2D with Data/MC ratio and corresponding error and also saves TH2D of the ratios inside the output file and (as pngs) inside ratioResults.

**Usage: ./produceRatio inputData_TGraph2D.root inputMC_TGraph2D.root mode output.root**

The first argument is output of _produceTGraphFromTh3_ for Data.
The second argument is output of _produceTGraphFromTh3_ for MC.
The third argument allows to specify which kind of TGraph2D to be used, currently woring with XTheta and XTheta_c (no lifetime correction and lifetime corrected).
The fourth argument is the output name.


TBD:
Add common analysis notebook and python scripts for final ratio plots.  --> Very Soon
Add binned analysis python script. --> Very Soon
Add other modalities. --> Not so Soon

For questions ask Francesco Poppi on Slack or by mail (poppi@bo.infn.it)