#include <iostream>
#include <vector>
#include <string>
#include "binning.h"

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <tpc> <bin_index> <min|max>\n";
        return 1;
    }
    initBins();
    int tpc = std::stoi(argv[1]);
    int bin_index = std::stoi(argv[2]);
    std::string opt = argv[3];

    std::vector<double> xedges;
    std::vector<double> thetaEdges;

    getAnalysisEdges_XThetaUnBinned(tpc, xedges, thetaEdges);

    if (bin_index < 0 || bin_index >= static_cast<int>(xedges.size()) - 1) {
        std::cerr << "Error: bin_index out of range (0-" << xedges.size()-2 << ")\n";
        return 1;
    }

    if (opt == "min") {
        std::cout << xedges[bin_index] << "\n";
    } else if (opt == "max") {
        std::cout << xedges[bin_index + 1] << "\n";
    } else {
        std::cerr << "Error: option must be 'min' or 'max'\n";
        return 1;
    }

    return 0;
}
