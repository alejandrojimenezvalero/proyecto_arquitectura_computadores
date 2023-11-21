#include "constants.hpp"
#include "grid.hpp"

#include "cmath"
using namespace simulationConstants;

void initGrid(Grid& grid, double smoothing_length) {
    double n_x=0;double n_y=0;double n_z=0;double s_x=0;double s_y=0;double s_z=0;
    if(smoothing_length > 0) {
        n_x = std::floor((UPPER_LIMIT[0] - LOWER_LIMIT[0]) / smoothing_length);
        n_y = std::floor((UPPER_LIMIT[1] - LOWER_LIMIT[1]) / smoothing_length);
        n_z = std::floor((UPPER_LIMIT[2] - LOWER_LIMIT[2]) / smoothing_length);
    }
    grid.grid_dimensions = {n_x, n_y, n_z};
    if (grid.grid_dimensions[0] > 0){s_x = (UPPER_LIMIT[0] - LOWER_LIMIT[0])/grid.grid_dimensions[0];}
    if (grid.grid_dimensions[1] > 0){s_y = (UPPER_LIMIT[1] - LOWER_LIMIT[1])/grid.grid_dimensions[1];}
    if (grid.grid_dimensions[2] > 0){s_z = (UPPER_LIMIT[2] - LOWER_LIMIT[2])/grid.grid_dimensions[2];}
    grid.block_dimensions = {s_x, s_y, s_z};
    grid.grid_blocks={};
}

std::vector<int> calcParticleIndex(Particle& particle, Grid& grid){
    std::vector<int> particle_block_index{};

    double index_i=0; double index_j=0; double index_k=0;

    if (grid.grid_dimensions[0]>0){index_i = std::floor((particle.pos[0] - LOWER_LIMIT[0]) / grid.block_dimensions[0]);}
    if (grid.grid_dimensions[1]>0){index_j = std::floor((particle.pos[1] - LOWER_LIMIT[1]) / grid.block_dimensions[1]);}
    if (grid.grid_dimensions[2]>0){index_k = std::floor((particle.pos[2] - LOWER_LIMIT[2]) / grid.block_dimensions[2]);}

    if (index_i < 0) {index_i = 0;} else if (index_i > grid.grid_dimensions[0] - 1) { index_i = grid.grid_dimensions[0] - 1;}
    if (index_j < 0) {index_j = 0;} else if (index_j > grid.grid_dimensions[1] - 1) {index_j = grid.grid_dimensions[1] - 1;}
    if (index_k < 0) {index_k = 0;} else if (index_k > grid.grid_dimensions[2] - 1) {index_k = grid.grid_dimensions[2] - 1;}

    particle_block_index = {static_cast<int>(index_i), static_cast<int>(index_j), static_cast<int>(index_k)};
    return particle_block_index;
}

int calcParticleIndexVector(Grid& grid, const std::vector<int>& block_cords){
    const int index_in_vector = grid.adjacent_index_map[block_cords];
    return index_in_vector;
}

bool blockExists(int i, int j, int k, Grid& grid){
    bool block_exist =true;
    if ((i < 0) or (i > grid.grid_dimensions[0] - 1) or (j < 0) or (i > grid.grid_dimensions[1] - 1) or (k < 0) or (k > grid.grid_dimensions[2] - 1)){block_exist = false;}
    return block_exist;
}