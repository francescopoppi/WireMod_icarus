#pragma once

#include <map>
#include <string>
#include <sstream> 
#include <fstream> 
#include <iostream>
#include <vector>

std::vector<std::string> getFileListFromRunList(const std::string& filename);

void loadLifetimeFile(const std::string& filename,
                      std::map<int,float>& tpc0,
                      std::map<int,float>& tpc1,
                      std::map<int,float>& tpc2,
                      std::map<int,float>& tpc3);