
#include "constants.hpp"
#include "block.hpp"
#include "grid.hpp"
#include "cmath"
using namespace simulationConstants;

blockSize calculateBlockSize(gridSize grid){
  blockSize block {0,0,0};
  block.sx = (UPPER_LIMIT[0] - LOWER_LIMIT[0])/grid.nx;
  block.sy = (UPPER_LIMIT[1] - LOWER_LIMIT[1])/grid.ny;
  block.sz = (UPPER_LIMIT[2] - LOWER_LIMIT[2])/grid.nz;
  return block;
};