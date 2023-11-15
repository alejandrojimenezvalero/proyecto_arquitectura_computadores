#ifndef BLOCK_HPP
#define BLOCK_HPP

#include "particle.hpp"
#include <memory>

struct Block{
    std::vector<int> block_index{};
    std::vector<Particle> block_particles;
    std::vector<std::vector<int>> adj_blocks;
    bool updated = false;
};

Block createBlock(int i, int j, int k);

#endif