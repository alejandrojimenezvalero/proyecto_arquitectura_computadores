//
#include "exceptionHandler.hpp"

#include <fstream>
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace std;

int validateParameters(const std::vector<std::string> &arguments) {
  size_t n_args = arguments.size();
  if(n_args != 3){throwException( "Error: Invalid number of arguments: " + std::to_string(n_args) + ".", -1);};
  try {
    if(stoi(arguments[0]) < 0){
      throwException("Error: Invalid number of time steps.",-2 );
    }
  } catch (const invalid_argument &) {
    throwException("Error: time steps must be numeric.", -1);
  }
  try {
    ifstream inputFile(arguments[1]);
    inputFile.close();
  } catch (const ifstream::failure &e) {
    throwException( "Cannot open " + arguments[1] + " for reading", -3);
  }
  try {
    ofstream outputFile(arguments[2]);
    outputFile.close();
  } catch (const ofstream::failure &e) {
    throwException( "Cannot open " + arguments[2] + " for writing", -4);
  }
  return 0;
  }

