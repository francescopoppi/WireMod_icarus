ROOTCFLAGS    := $(shell root-config --cflags)
ROOTLIBS      := $(shell root-config --libs)
ROOTGLIBS     := $(shell root-config --glibs)

CXX           := g++
CXXFLAGS      := -O2 -Wall -fPIC -std=c++17

SRC           := ntupleProcessing/ntupleAnalyzer.cpp \
                 CommonTools/TrackReader.cpp \
                 CommonTools/Corrections.cpp \
                 CommonTools/Utils.cpp \
                 CommonTools/HistogramManager.cpp \
                 CommonTools/binning.cpp

TARGET        := ntupleAnalyzer
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(SRC) -o $(TARGET) $(ROOTLIBS)

clean:
	rm -f $(TARGET) *.o
