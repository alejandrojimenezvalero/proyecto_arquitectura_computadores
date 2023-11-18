#ifndef GRID_HPP
#define GRID_HPP

#include "particle.hpp"
#include "block.hpp"

#include <vector>
#include <map>

struct Grid{
    std::vector<double> grid_dimensions;
    std::vector<double> block_dimensions;
    std::map<std::vector<int>,Block> grid_blocks;
};


void calculateBlockSize(Grid& grid);
void initGrid(Grid& grid, double smoothing_length);
std::vector<int> calcParticleIndex(Particle& particle, Grid& grid);
#endif