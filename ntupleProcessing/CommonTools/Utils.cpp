#include "Utils.h"
#include <fstream>
#include <iostream>

std::vector<std::string> getFileListFromRunList(const std::string& filename) {
    std::vector<std::string> files;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Cannot open run list file: " << filename << std::endl;
        return files;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (!line.empty()) files.push_back(line);
    }

    infile.close();
    return files;
}


void loadLifetimeFile(const std::string& filename,
                      std::map<int,float>& tpc0,
                      std::map<int,float>& tpc1,
                      std::map<int,float>& tpc2,
                      std::map<int,float>& tpc3)
{
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error opening lifetime file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        std::istringstream ss(line);
        std::string token;
        int run;
        float v0, v1, v2, v3;

        try {
            std::getline(ss, token, ';'); run = std::stoi(token);
            std::getline(ss, token, ';'); v0 = std::stof(token);
            std::getline(ss, token, ';'); v1 = std::stof(token);
            std::getline(ss, token, ';'); v2 = std::stof(token);
            std::getline(ss, token, ';'); v3 = std::stof(token);
        }
        catch (const std::exception& e) {
            std::cerr << "Error parsing line in " << filename << ": " << e.what() << std::endl;
            continue;
        }

        tpc0[run] = v0;
        tpc1[run] = v1;
        tpc2[run] = v2;
        tpc3[run] = v3;
    }

    infile.close();
}