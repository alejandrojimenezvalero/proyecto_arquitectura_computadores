
#include "constants.hpp"
#include "block.hpp"
#include "grid.hpp"
#include "particle.hpp"

#include <memory>

using namespace simulationConstants;

std::vector<double> calculateBlockSize(gridSize grid){
  double sx = 0,sy = 0,sz = 0;
    if(grid.nx >0){sx = (UPPER_LIMIT[0] - LOWER_LIMIT[0])/grid.nx;}
    if(grid.ny >0){sy = (UPPER_LIMIT[1] - LOWER_LIMIT[1])/grid.ny;}
    if(grid.nz >0){sz = (UPPER_LIMIT[2] - LOWER_LIMIT[2])/grid.nz;}
    std::vector<double> block_dimensions = {sx, sy, sz};
    return block_dimensions;
};

std::shared_ptr<std::vector<Particle>> initializeBlockParticles() {
    return std::make_shared<std::vector<Particle>>();
}

void insertParticle(Block& block, const Particle& particle) {
    if (!block.block_particles) {
        // Si el vector de partículas no ha sido inicializado o está vacío, inicializarlo ahora.
        block.block_particles = initializeBlockParticles();
    }
    // Agregar la partícula al vector de partículas
    block.block_particles->push_back(particle);
}

Block createBlock(int i, int j, int k) {
    Block block {{i, j, k}, nullptr};  // Inicializar con un puntero nulo
    return block;
}

std::vector<int> calcParticleIndex(Particle& particle, std::vector<double> block_dimensions){
    std::vector<int> particle_block_index{};
    int i, j, k;
    i = std::floor((particle.px - LOWER_LIMIT[0])/block_dimensions[0]);
    j = std::floor((particle.py - LOWER_LIMIT[1])/block_dimensions[1]);
    k = std::floor((particle.pz - LOWER_LIMIT[2])/block_dimensions[2]);
    particle_block_index = {i, j, k};
    return particle_block_index;
}