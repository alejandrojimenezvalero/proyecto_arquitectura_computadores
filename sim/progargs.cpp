//
#include "exceptionHandler.hpp"

#include <fstream>
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace std;

int validateParameters(const std::vector<std::string> &arguments) {
  size_t n_args = arguments.size();
  exceptionHandler(n_args != 3,"Error: Invalid number of arguments: " + std::to_string(n_args) + ".", -1);
  ifstream inputFile(arguments[1]); ofstream outputFile(arguments[2]);
  try {
    exceptionHandler(stoi(arguments[0]) < 0, "Error: Invalid number of time steps.", -2);
  } catch (const invalid_argument &) {throwException("Error: time steps must be numeric.", -1);}
  exceptionHandler(!inputFile, "Cannot open " + arguments[1] + " for reading", -3);
  exceptionHandler(!outputFile, "Cannot open " + arguments[2] + " for writing", -4);
  inputFile.close();
  outputFile.close();
  return 0;
}
