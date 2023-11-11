
#include "constants.hpp"
#include "grid.hpp"

#include "cmath"
using namespace simulationConstants;

gridSize calculateGridSize(double smoothing_length){
  gridSize grid {0,0,0, {}};
  //Usamos una función techo para redondear hacia arriba
  grid.nx = std::floor((UPPER_LIMIT[0] - LOWER_LIMIT[0])/smoothing_length);
  grid.ny = std::floor((UPPER_LIMIT[1] - LOWER_LIMIT[1])/smoothing_length);
  grid.nz = std::floor((UPPER_LIMIT[2] - LOWER_LIMIT[2])/smoothing_length);
  return grid;
}