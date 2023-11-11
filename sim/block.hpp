#ifndef BLOCK_HPP
#define BLOCK_HPP

#include "grid.hpp"
#include "particle.hpp"
#include <memory>

struct Block{
    std::vector<int> block_index{};
    std::shared_ptr<std::vector<Particle>> block_particles;
};

std::vector<double> calculateBlockSize(gridSize grid);
std::shared_ptr<std::vector<Particle>> initializeBlockParticles();
void insertParticle(Block& block, const Particle& particle);
Block createBlock(int i, int j, int k);
std::vector<int> calcParticleIndex(Particle& particle, std::vector<double> block_dimensions);

#endif