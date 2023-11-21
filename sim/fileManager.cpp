//
#include <iostream>
#include "fileManager.hpp"
#include <fstream>

std::ifstream openFile(const std::string &inputFile) {
    std::ifstream input_file(inputFile);
    return input_file;
}