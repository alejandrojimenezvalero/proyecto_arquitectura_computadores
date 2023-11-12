#ifndef GRID_HPP
#define GRID_HPP

#include "particle.hpp"
#include "block.hpp"

#include <vector>
#include <map>

struct Grid{
    std::vector<double> ngrid;
    std::vector<double> block_dimensions;
    std::vector<Block> particleMap;
};

struct SubGrid {
    Block& current_block;
    std::vector<std::reference_wrapper<Block>> particleSubMap;

    SubGrid(Block& initialBlock)
            : current_block(initialBlock), particleSubMap() {}
};



void calculateBlockSize(Grid& grid);
void initGrid(Grid& grid, double smoothing_length);
#endif