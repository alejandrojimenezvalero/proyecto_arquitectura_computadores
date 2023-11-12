#ifndef BLOCK_HPP
#define BLOCK_HPP

#include "particle.hpp"
#include <memory>

struct Block{
    std::vector<int> block_index{};
    std::vector<Particle> block_particles;
};

Block createBlock(int i, int j, int k);
std::vector<int> calcParticleIndex(Particle& particle, std::vector<double> block_dimensions);

#endif