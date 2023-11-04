//
#include <iostream>
#include "fileManager.hpp"
#include <fstream>
#include "exceptionHandler.hpp"
std::ifstream openFile(const std::string& inputFile){
  std::ifstream input_file(inputFile);
  exceptionHandler(!input_file, "Cannot open " + inputFile + " for reading", -3);
  return input_file;
}