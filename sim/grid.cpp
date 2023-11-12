
#include "constants.hpp"
#include "grid.hpp"

#include "cmath"
using namespace simulationConstants;

void initGrid(Grid& grid, double smoothing_length) {
    // Actualiza los miembros de grid según smoothing_length
    double nx, ny, nz;
    nx = std::floor((UPPER_LIMIT[0] - LOWER_LIMIT[0]) / smoothing_length);
    ny = std::floor((UPPER_LIMIT[1] - LOWER_LIMIT[1]) / smoothing_length);
    nz = std::floor((UPPER_LIMIT[2] - LOWER_LIMIT[2]) / smoothing_length);
    grid.ngrid = {nx, ny, nz};
    grid.block_dimensions={};
    grid.particleMap={};
}

void calculateBlockSize(Grid& grid){
    double sx = (UPPER_LIMIT[0] - LOWER_LIMIT[0])/grid.ngrid[0];
    double sy = (UPPER_LIMIT[1] - LOWER_LIMIT[1])/grid.ngrid[1];
    double sz = (UPPER_LIMIT[2] - LOWER_LIMIT[2])/grid.ngrid[2];
    grid.block_dimensions = {sx, sy, sz};
}