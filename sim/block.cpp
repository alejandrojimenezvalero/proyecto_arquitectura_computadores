
#include "constants.hpp"
#include "block.hpp"
#include "particle.hpp"

#include <memory>

using namespace simulationConstants;

Block createBlock(int i, int j, int k) {
    Block block {{i, j, k},{}};  // Inicializar con un puntero nulo
    return block;
}

std::vector<int> calcParticleIndex(Particle& particle, std::vector<double> block_dimensions){
    std::vector<int> particle_block_index{};
    int i, j, k;

    i = std::floor((particle.pos[0] - LOWER_LIMIT[0])/block_dimensions[0]);
    j = std::floor((particle.pos[1] - LOWER_LIMIT[1])/block_dimensions[1]);
    k = std::floor((particle.pos[2] - LOWER_LIMIT[2])/block_dimensions[2]);

    particle_block_index = {i, j, k};
    return particle_block_index;
}