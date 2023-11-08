#ifndef BLOCK_HPP
#define BLOCK_HPP

#include "grid.hpp"
#include "particle.hpp"

struct blockSize{
    double sx;
    double sy;
    double sz;
};

blockSize calculateBlockSize(gridSize grid);
void calcParticleIndex(Particle particle, blockSize block);

#endif