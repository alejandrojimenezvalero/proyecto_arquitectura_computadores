
#include "constants.hpp"
#include "grid.hpp"

#include "cmath"
using namespace simulationConstants;

void initGrid(Grid& grid, double smoothing_length) {
    // Actualiza los miembros de grid según smoothing_length
    double nx, ny, nz;
    if(smoothing_length > 0) {
      nx = std::floor((UPPER_LIMIT[0] - LOWER_LIMIT[0]) / smoothing_length);
      ny = std::floor((UPPER_LIMIT[1] - LOWER_LIMIT[1]) / smoothing_length);
      nz = std::floor((UPPER_LIMIT[2] - LOWER_LIMIT[2]) / smoothing_length);
    }else{nx =0, ny = 0, nz= 0;}
    grid.grid_dimensions = {nx, ny, nz};
    grid.block_dimensions={};
    grid.grid_blocks={};
}

void calculateBlockSize(Grid& grid){
    double sx = 0, sy = 0, sz = 0;
    if (grid.grid_dimensions[0] > 0){sx = (UPPER_LIMIT[0] - LOWER_LIMIT[0])/grid.grid_dimensions[0];}
    if (grid.grid_dimensions[1] > 0){sy = (UPPER_LIMIT[1] - LOWER_LIMIT[1])/grid.grid_dimensions[1];}
    if (grid.grid_dimensions[2] > 0){sz = (UPPER_LIMIT[2] - LOWER_LIMIT[2])/grid.grid_dimensions[2];}

    grid.block_dimensions = {sx, sy, sz};
}
std::vector<int> calcParticleIndex(Particle& particle, Grid& grid){
    std::vector<int> particle_block_index{};
    double i =0, j=0, k=0;

    if (grid.grid_dimensions[0]>0){i = std::floor((particle.pos[0] - LOWER_LIMIT[0]) / grid.grid_dimensions[0]);}
    if (grid.grid_dimensions[1]>0){j = std::floor((particle.pos[1] - LOWER_LIMIT[1]) / grid.grid_dimensions[1]);}
    if (grid.grid_dimensions[2]>0){k = std::floor((particle.pos[2] - LOWER_LIMIT[2]) / grid.grid_dimensions[2]);}

    if (i < 0) {i = 0;} else if (i > grid.grid_dimensions[0] - 1) { i = grid.grid_dimensions[0] - 1;}

    if (j < 0) {j = 0;} else if (j > grid.grid_dimensions[1] - 1) { j = grid.grid_dimensions[1] - 1;}

    if (k < 0) {k = 0;} else if (k > grid.grid_dimensions[2] - 1) { k = grid.grid_dimensions[2] - 1;}

    particle_block_index = {static_cast<int>(i), static_cast<int>(j), static_cast<int>(k)};
    return particle_block_index;
}