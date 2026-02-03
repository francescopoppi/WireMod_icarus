ROOTCFLAGS    := $(shell root-config --cflags)
ROOTLIBS      := $(shell root-config --libs)
ROOTGLIBS     := $(shell root-config --glibs)

CXX           := g++
CXXFLAGS      := -O2 -Wall -fPIC -std=c++17 -I./CommonTools

SRC_ANALYZER  := ntupleProcessing/ntupleAnalyzer.cpp \
                 CommonTools/TrackReader.cpp \
                 CommonTools/Corrections.cpp \
                 CommonTools/Utils.cpp \
                 CommonTools/HistogramManager.cpp \
                 CommonTools/binning.cpp

SRC_TH3       := BinnedRatioMaker/processTH3.cpp \
                 CommonTools/AnalysisTools.cpp \
                 CommonTools/binning.cpp

SRC_RATIO     := BinnedRatioMaker/processTGraphRatio.cpp \
                 CommonTools/binning.cpp

SRC_SLICE := unBinnedAnalyzer/split_slice.cpp \
             CommonTools/binning.cpp

SRC_SLICEYZ := unBinnedAnalyzer/split_slice_z.cpp \
             CommonTools/binning.cpp


SRC_BINEDGE := unBinnedAnalyzer/returnEdges.cpp \
             CommonTools/binning.cpp


TARGET_ANALYZER := ntupleAnalyzer
TARGET_TH3      := produceTGraphFromTH3
TARGET_RATIO    := produceRatio
TARGET_SLICE    := produceSlice
TARGET_SLICEYZ  := produceSliceYZ

TARGET_BINEDGE  := returnEdge

all: $(TARGET_ANALYZER) $(TARGET_TH3) $(TARGET_RATIO) $(TARGET_SLICE) $(TARGET_SLICEYZ) $(TARGET_BINEDGE)

$(TARGET_ANALYZER): $(SRC_ANALYZER)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(SRC_ANALYZER) -o $(TARGET_ANALYZER) $(ROOTLIBS)

$(TARGET_TH3): $(SRC_TH3)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(SRC_TH3) -o $(TARGET_TH3) $(ROOTLIBS)

$(TARGET_RATIO): $(SRC_RATIO)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(SRC_RATIO) -o $(TARGET_RATIO) $(ROOTLIBS)

$(TARGET_SLICE): $(SRC_SLICE)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(SRC_SLICE) -o $(TARGET_SLICE) $(ROOTLIBS)

$(TARGET_SLICEYZ): $(SRC_SLICEYZ)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(SRC_SLICEYZ) -o $(TARGET_SLICEYZ) $(ROOTLIBS)

$(TARGET_BINEDGE): $(SRC_BINEDGE)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(SRC_BINEDGE) -o $(TARGET_BINEDGE) $(ROOTLIBS)

clean:
	rm -f $(TARGET_ANALYZER) $(TARGET_TH3) $(TARGET_RATIO) $(TARGET_SLICE) $(TARGET_SLICEYZ) $(TARGET_BINEDGE) *.o
