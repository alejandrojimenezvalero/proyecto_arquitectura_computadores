#include "constants.hpp"
#include "grid.hpp"
using namespace simulationConstants;

gridSize calculateGridSize(double smoothing_length){
  gridSize grid {0,0,0};
  grid.nx = (UPPER_LIMIT[0] - LOWER_LIMIT[0])/smoothing_length;
  grid.ny = (UPPER_LIMIT[1] - LOWER_LIMIT[1])/smoothing_length;
  grid.nz = (UPPER_LIMIT[2] - LOWER_LIMIT[2])/smoothing_length;
  return grid;
};