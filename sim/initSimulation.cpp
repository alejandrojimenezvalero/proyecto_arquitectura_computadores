//
#include "initSimulation.hpp"

#include "block.hpp"
#include "constants.hpp"
#include "exceptionHandler.hpp"
#include "fileManager.hpp"
#include "grid.hpp"
#include "particle.hpp"
#include "processSimulation.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <chrono>

using namespace simulationConstants;

/*
void createAdjacentBlocks(Block& block, Grid grid) {
    int i = block.block_index[0];
    int j = block.block_index[1];
    int k = block.block_index[2];
    for (int di = -1; di <= 1; di++) {
        for (int dj = -1; dj <= 1; dj++) {
            for (int dk = -1; dk <= 1; dk++) {
                if (blockExists(i + di, j + dj, k + dk, grid)) {
                    // Agregar las coordenadas del bloque adyacente al vector
                    block.adj_blocks.push_back({i + di, j + dj, k + dk});
                }
            }
        }
    }
}
void createGridBlocks(Grid& grid) {
    for (int i = 0; i < grid.grid_dimensions[0]; ++i) {
        for (int j = 0; j < grid.grid_dimensions[1]; ++j) {
            for (int k = 0; k < grid.grid_dimensions[2]; ++k) {
                Block block = createBlock(i, j, k);
                createAdjacentBlocks(block, grid);

                grid.grid_blocks.push_back(block);
            }
        }
    }
}*/
void createGridBlocks(Grid& grid) {
    int cont = 0;
    for (int i =0; i < grid.grid_dimensions[0]; ++i) {
        for (int j=0; j < grid.grid_dimensions[1]; ++j) {
            for (int k=0; k < grid.grid_dimensions[2]; ++k){
                Block block = createBlock(i, j, k);
                grid.adjacent_index_map[{i,j,k}]= cont;
                cont++;
                for (int di = -1; di <= 1; di++) {
                    for (int dj = -1; dj <= 1; dj++) {
                        for (int dk = -1; dk <= 1; dk++) {
                            const int target_i = i + di;
                            const int target_j = j + dj;
                            const int target_k = k + dk;

                            if (blockExists(target_i, target_j, target_k, grid)) {
                                // Agregar las coordenadas del bloque adyacente al vector
                                block.adj_blocks_cords.push_back({target_i, target_j, target_k});
                            }
                        }
                    }
                }
                grid.grid_blocks.push_back(block);
            }
        }
    }
}

void initAdjIndexVectorBlocks(Grid& grid) {
    for (Block& block : grid.grid_blocks) {
        for(const std::vector<int>& adj_cords : block.adj_blocks_cords) {
            block.adj_index_vector_blocks.push_back(grid.adjacent_index_map[adj_cords]);
        }
    }
}


void checkBlockIndex(int &index_i, int & index_j, int & index_k, Grid& grid){
    index_i = (index_i < 0)? 0:index_i;
    index_i = (index_i > grid.grid_dimensions[0] - 1) ? static_cast<int>(grid.grid_dimensions[0]) - 1 : index_i;
    index_j = (index_j < 0)? 0: index_j;
    index_j = (index_j > grid.grid_dimensions[1] - 1) ? static_cast<int>(grid.grid_dimensions[1]) - 1 : index_j;
    index_k = (index_k < 0)? 0: index_k;
    index_k = (index_k > grid.grid_dimensions[2] - 1) ? static_cast<int>(grid.grid_dimensions[2]) - 1 : index_k;
}
void calculateParameters(double ppm, int np, SimulationData& data) {
    const double smoothing_length = RADIO_MULTIPLICATOR / ppm;
    const double particle_mass = FLUID_DENSITY / pow(ppm, 3);
    initGrid(data.grid, smoothing_length);
    calculateBlockSize(data.grid);
    data.smoothing_length= smoothing_length; data.particle_mass= particle_mass;
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    data.smoothing_length_2 = pow(data.smoothing_length, 2);
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    data.smoothing_length_6 = pow(smoothing_length, 6);
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    data.smoothing_length_9 = pow(smoothing_length, 9);
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    data.escalar_pos  = (15 / (simulationConstants::PI * data.smoothing_length_6)) * ((3 * particle_mass * simulationConstants::STIFFNESS_PRESSURE) / 2);
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    data.escalar_vel = (45 / (simulationConstants::PI * data.smoothing_length_6));
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    data.escalar_density = (315/(64*PI*data.smoothing_length_9)) * data.particle_mass;

    std::cout << "Number of particles: " << np << '\n';
    std::cout << "Particles per meter: " << ppm << '\n';
    std::cout << "Smoothing length: " << smoothing_length << '\n';
    std::cout << "Particle mass: " << particle_mass << '\n';
    std::cout << "Grid size: " << data.grid.grid_dimensions[0] << " x " << data.grid.grid_dimensions[1] << " x " << data.grid.grid_dimensions[2] << '\n';
    std::cout << "Number of blocks: " << data.grid.grid_dimensions[0] * data.grid.grid_dimensions[1] * data.grid.grid_dimensions[2] << '\n';
    std::cout << "Block size: " << data.grid.block_dimensions[0] << " x " << data.grid.block_dimensions[1] << " x " << data.grid.block_dimensions[2] << '\n';
}

int setParticleData(const std::string& inputFile, SimulationData& data){
    std::ifstream input_file = openFile(inputFile);
    createGridBlocks(data.grid);
    initAdjIndexVectorBlocks(data.grid);
    std::cout << "check0" << '\n';
    int real_particles = 0;
    //NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    input_file.seekg(8,std::ios::beg);
    while(!input_file.eof()){
        Particle particle;
        readParticleFields(input_file, particle);
        //Calculamos los indices
        std::vector<int> particle_block_index = calcParticleIndex(particle, data.grid);

        checkBlockIndex(particle_block_index[0], particle_block_index[1], particle_block_index[2], data.grid);

        if (input_file.eof()) {
          break;
        }
        particle.id = real_particles;
        addParticleToBlock(particle, data.grid, particle_block_index);

        ++real_particles;
        }
    return real_particles;
}

void readParticleFields(std::ifstream& input_file, Particle& particle) {
    float floatField = 0;
    for (const auto& field : {particle.pos.data(), &particle.pos[1], &particle.pos[2], particle.hv.data(), &particle.hv[1], &particle.hv[2], particle.vel.data(), &particle.vel[1], &particle.vel[2]}) {
        //NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
        input_file.read(reinterpret_cast<char*>(&floatField), sizeof(float));
        *field = static_cast<double>(floatField);
    }
}

void addParticleToBlock(const Particle& particle, Grid& grid, const std::vector<int>& particle_block_index) {
    for (Block& block : grid.grid_blocks) {
        if (block.block_index == particle_block_index) {
            block.block_particles.push_back(particle);
        }
    }
}

int initiateSimulation(const std::string& n_iterations, const std::string& inputFile) {
    auto start = std::chrono::high_resolution_clock::now();
    std::ifstream input_file = openFile(inputFile);
    if(!input_file){throwException("Cannot open " + inputFile + " for reading", -3);}

    const int n_iterations_int = std::stoi(n_iterations);

    float ppmFloat = 0.0;
    input_file.read(reinterpret_cast<char *>(&ppmFloat), sizeof(float)); //NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
    const auto ppm = static_cast<double>(ppmFloat);

    int num_particles = 0;
    input_file.read(reinterpret_cast<char *>(&num_particles), sizeof(int)); //NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
    if (num_particles < 0) {throwException("Error: Invalid number of particles: " + std::to_string(num_particles), -5);} //NOLINT(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)

    Grid grid;
    SimulationData data(grid);
    calculateParameters(ppm, num_particles, data);

    const int real_particles = setParticleData(inputFile, data);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Init: " << duration.count() << '\n';
    if (num_particles != real_particles) {throwException("Error: Number of particles mismatch. Header:  " + std::to_string(num_particles) +", Found: " + std::to_string(real_particles),-5); } //NOLINT(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)

    for (int i = 0; i < n_iterations_int; ++i) {processSimulation(data);}

    return 0;
}
