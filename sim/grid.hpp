#ifndef GRID_HPP
#define GRID_HPP

#include "particle.hpp"
#include "block.hpp"

#include <vector>
#include <map>

struct Grid{
    std::vector<double> grid_dimensions;
    std::vector<double> block_dimensions;
    std::vector<Block> grid_blocks;
    std::map<std::vector<int>,int> adjacent_index_map;
};


void calculateBlockSize(Grid& grid);
void initGrid(Grid& grid, double smoothing_length);
std::vector<int> calcParticleIndex(Particle& particle, Grid& grid);
int calcParticleIndexVector(Particle& particle, Grid& grid);
bool blockExists(int i, int j, int k, Grid& grid);
#endif