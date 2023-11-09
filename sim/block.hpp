#ifndef BLOCK_HPP
#define BLOCK_HPP

#include "grid.hpp"
#include "particle.hpp"

struct Block{
    std::vector<int> block_index{};
    std::vector<Particle> block_particles;
};

std::vector<double> calculateBlockSize(gridSize grid);
Block createBlock(int i, int j, int k);
std::vector<int> calcParticleIndex(Particle particle, std::vector<double> block_dimensions);

#endif