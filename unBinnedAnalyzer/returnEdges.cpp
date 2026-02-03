#include <iostream>
#include <vector>
#include <string>
#include "binning.h"
#include <iostream>
#include <vector>
#include <string>
#include "binning.h"

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <X|Z|Z2> <tpc> <bin_index> <min|max>\n";
        return 1;
    }

    initBins();

    std::string mode = argv[1];
    int tpc = std::stoi(argv[2]);
    int bin_index = std::stoi(argv[3]);
    std::string opt = argv[4];

    if (mode == "X") {
        std::vector<double> xedges;
        std::vector<double> thetaEdges;

        getAnalysisEdges_XThetaUnBinned(tpc, xedges, thetaEdges);

        if (bin_index < 0 || bin_index >= static_cast<int>(xedges.size()) - 1) {
            std::cerr << "Error: X bin_index out of range (0-" << xedges.size()-2 << ")\n";
            return 1;
        }

        if (opt == "min")
            std::cout << xedges[bin_index] << "\n";
        else if (opt == "max")
            std::cout << xedges[bin_index + 1] << "\n";
        else {
            std::cerr << "Error: option must be 'min' or 'max'\n";
            return 1;
        }
    }

    else if (mode == "Z") {

        if (bin_index < 0 || bin_index >= nBinsZ) {
            std::cerr << "Error: Z bin_index out of range (0-" << nBinsZ-1 << ")\n";
            return 1;
        }

        if (opt == "min")
            std::cout << binsZ[bin_index] << "\n";
        else if (opt == "max")
            std::cout << binsZ[bin_index + 1] << "\n";
        else {
            std::cerr << "Error: option must be 'min' or 'max'\n";
            return 1;
        }
    }

    else if (mode == "Z2") {

        if (bin_index < 0 || bin_index >= nBinsZ_2) {
            std::cerr << "Error: Z bin_index out of range (0-" << nBinsZ_2-1 << ")\n";
            return 1;
        }

        if (opt == "min")
            std::cout << binsZ_2[bin_index] << "\n";
        else if (opt == "max")
            std::cout << binsZ_2[bin_index + 1] << "\n";
        else {
            std::cerr << "Error: option must be 'min' or 'max'\n";
            return 1;
        }
    }

    else {
        std::cerr << "Error: first argument must be 'X' or 'Z' or 'Z2'\n";
        return 1;
    }

    return 0;
}
