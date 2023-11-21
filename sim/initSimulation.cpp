//
#include "initSimulation.hpp"

#include "block.hpp"
#include "constants.hpp"
#include "exceptionHandler.hpp"
#include "fileManager.hpp"
#include "grid.hpp"
#include "particle.hpp"
#include "processSimulation.hpp"
#include "simulationResults.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <chrono>

using namespace simulationConstants;


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


void checkBlockIndex(std::vector<int>& particle_block_index, Grid& grid){
    particle_block_index[0] = (particle_block_index[0] < 0)? 0:particle_block_index[0];
    particle_block_index[0] = (particle_block_index[0] > grid.grid_dimensions[0] - 1) ? static_cast<int>(grid.grid_dimensions[0]) - 1 : particle_block_index[0];
    particle_block_index[1] = (particle_block_index[1] < 0)? 0: particle_block_index[1];
    particle_block_index[1] = (particle_block_index[1] > grid.grid_dimensions[1] - 1) ? static_cast<int>(grid.grid_dimensions[1]) - 1 : particle_block_index[1];
    particle_block_index[2] = (particle_block_index[2] < 0)? 0: particle_block_index[2];
    particle_block_index[2] = (particle_block_index[2] > grid.grid_dimensions[2] - 1) ? static_cast<int>(grid.grid_dimensions[2]) - 1 : particle_block_index[2];
}
void initializeData(SimulationData& data, const double smoothing_length, const double particle_mass){
    data.smoothing_length= smoothing_length;data.particle_mass= particle_mass;
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    data.smoothing_length_2 = pow(data.smoothing_length, 2);
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    data.smoothing_length_6 = pow(data.smoothing_length_2, 3);
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    data.smoothing_length_9 = pow(smoothing_length, 9);
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    data.escalar_pos  = (15 / (simulationConstants::PI * data.smoothing_length_6)) * ((3 * particle_mass * simulationConstants::STIFFNESS_PRESSURE) / 2);
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    data.escalar_vel = (45 / (simulationConstants::PI * data.smoothing_length_6));
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    data.escalar_density = (315/(64*PI*data.smoothing_length_9)) * data.particle_mass;
}
void calculateParameters(double ppm, int np, SimulationData& data) {
    const double smoothing_length = RADIO_MULTIPLICATOR / ppm;
    const double particle_mass = FLUID_DENSITY / pow(ppm, 3);
    initGrid(data.grid, smoothing_length);
    data.np = np; data.ppm = ppm;

    initializeData(data, smoothing_length, particle_mass);

    std::cout << "Number of particles: " << np << '\n';
    std::cout << "Particles per meter: " << ppm << '\n';
    std::cout << "Smoothing length: " << smoothing_length << '\n';
    std::cout << "Particle mass: " << particle_mass << '\n';
    std::cout << "Grid size: " << data.grid.grid_dimensions[0] << " x " << data.grid.grid_dimensions[1] << " x " << data.grid.grid_dimensions[2] << '\n';
    std::cout << "Number of blocks: " << data.grid.grid_dimensions[0] * data.grid.grid_dimensions[1] * data.grid.grid_dimensions[2] << '\n';
    std::cout << "Block size: " << data.grid.block_dimensions[0] << " x " << data.grid.block_dimensions[1] << " x " << data.grid.block_dimensions[2] << '\n';
}

void readParticleFields(std::ifstream& input_file, Particle& particle) {
    float floatField = 0;
    for (const auto& field : {particle.pos.data(), &particle.pos[1], &particle.pos[2], particle.hv.data(), &particle.hv[1], &particle.hv[2], particle.vel.data(), &particle.vel[1], &particle.vel[2]}) {
        //NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
        input_file.read(reinterpret_cast<char*>(&floatField), sizeof(float));
        *field = static_cast<double>(floatField);
    }
}


int setParticleData(std::ifstream& input_file, SimulationData& data){
    //std::ifstream input_file = openFile(inputFile);
    createGridBlocks(data.grid);
    initAdjIndexVectorBlocks(data.grid);
    //std::cout << "check0" << '\n';
    int real_particles = 0;
    //NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    while(!input_file.eof()){
        Particle particle;
        readParticleFields(input_file, particle);
        if (input_file.eof()) {
            break;
        }
        //Calculamos los indices
        std::vector<int> particle_block_index = calcParticleIndex(particle, data.grid);
        checkBlockIndex(particle_block_index, data.grid);
        const int particle_block_vector_index = calcParticleIndexVector(data.grid, particle_block_index);
        particle.id = real_particles;
        data.grid.grid_blocks[particle_block_vector_index].block_particles_v1.push_back(particle);
        ++real_particles;
        }
    if (data.np != real_particles) {throwException("Error: Number of particles mismatch. Header:  " + std::to_string(data.np) +", Found: " + std::to_string(real_particles),-5); } //NOLINT(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)

    return 0;
}

int initiateSimulation(const std::vector<std::string>& stringInitVector) {
    auto start = std::chrono::high_resolution_clock::now();
    std::ifstream input_file = openFile(stringInitVector[1]);
    const int n_iterations_int = std::stoi(stringInitVector[0]);
    float ppmFloat = 0.0;
    input_file.read(reinterpret_cast<char *>(&ppmFloat), sizeof(float)); //NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
    const auto ppm = static_cast<double>(ppmFloat);
    int num_particles = 0;
    input_file.read(reinterpret_cast<char *>(&num_particles), sizeof(int)); //NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
    if (num_particles < 0) {throwException("Error: Invalid number of particles: " + std::to_string(num_particles), -5);} //NOLINT(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
    SimulationData data;
    calculateParameters(ppm, num_particles, data);

    setParticleData(input_file, data);
    input_file.close();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Init: " << duration.count() << '\n';

    for (int i = 0; i < n_iterations_int; ++i) {processSimulation(data);}
    writeParticleData(stringInitVector[2], data);
    return 0;
}
