//
#include <iostream>
#include <string>

#include <cstdlib>  // For exit()

void throwException(const std::string &errorMessage, int exitCode) {
    std::cerr << errorMessage << "\n";
    exit(exitCode);
}