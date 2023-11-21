#ifndef BLOCK_HPP
#define BLOCK_HPP

#include "particle.hpp"
#include <memory>

struct Block {
    std::vector<int> block_index{};
    std::vector<Particle> block_particles_v1;
    std::vector<Particle> block_particles_v2;
    std::vector<int> adj_index_vector_blocks;
    std::vector<std::vector<int>> adj_blocks_cords;
    bool updated_density = false;
    bool updated_acceleration = false;
};

Block createBlock(int i, int j, int k);

#endif