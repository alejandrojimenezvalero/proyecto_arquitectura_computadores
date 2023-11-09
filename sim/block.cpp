
#include "constants.hpp"
#include "block.hpp"
#include "grid.hpp"
#include "particle.hpp"

using namespace simulationConstants;




std::vector<double> calculateBlockSize(gridSize grid){
    double sx = (UPPER_LIMIT[0] - LOWER_LIMIT[0])/grid.nx;
    double sy = (UPPER_LIMIT[1] - LOWER_LIMIT[1])/grid.ny;
    double sz = (UPPER_LIMIT[2] - LOWER_LIMIT[2])/grid.nz;
    std::vector<double> block_dimensions = {sx, sy, sz};
    return block_dimensions;
};


Block createBlock(int i, int j, int k){
    Block block {{0,0,0}, {}};
    block.block_index[0] = i; block.block_index[1] = j, block.block_index[2] = k;
    return block;
}

std::vector<int> calcParticleIndex(Particle particle, std::vector<double> block_dimensions){
    std::vector<int> particle_block_index{};
    int i, j, k;
    i = std::floor((particle.px - LOWER_LIMIT[0])/block_dimensions[0]);
    j = std::floor((particle.py - LOWER_LIMIT[1])/block_dimensions[1]);
    k = std::floor((particle.pz - LOWER_LIMIT[2])/block_dimensions[2]);
    particle_block_index = {i, j, k};
    return particle_block_index;
}