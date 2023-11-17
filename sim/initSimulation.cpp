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
#include <map>
#include <tuple>
#include <memory>

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
    for (int i =0; i < grid.grid_dimensions[0]; ++i) {
        for (int j=0; j < grid.grid_dimensions[1]; ++j) {
            for (int k=0; k < grid.grid_dimensions[2]; ++k){
                Block block = createBlock(i, j, k);
                for (int di = -1; di <= 1; di++) {
                    for (int dj = -1; dj <= 1; dj++) {
                        for (int dk = -1; dk <= 1; dk++) {
                            int target_i = i + di;
                            int target_j = j + dj;
                            int target_k = k + dk;

                            if (blockExists(target_i, target_j, target_k, grid)) {
                                // Agregar las coordenadas del bloque adyacente al vector
                                block.adj_blocks.push_back({target_i, target_j, target_k});
                            }
                        }
                    }
                }
                grid.grid_blocks[block.block_index]  = block;
            }
        }
    }
}


int checkBlockIndex(int &i, int &j, int &k, Grid& grid){
  i = (i < 0)? 0:i;
  i = (i > grid.grid_dimensions[0] - 1) ? grid.grid_dimensions[0] - 1 : i;
  j = (j < 0)? 0:j;
  j = (j > grid.grid_dimensions[1] - 1) ? grid.grid_dimensions[1] - 1 : j;
  k = (k < 0)? 0:k;
  k = (k > grid.grid_dimensions[2] - 1) ? grid.grid_dimensions[2] - 1 : k;
  return 0;
}
void calculateParameters(double ppm, int np, SimulationData& data) {
  double smoothing_length = RADIO_MULTIPLICATOR / ppm;
  double particle_mass = FLUID_DENSITY / pow(ppm, 3);
  initGrid(data.grid, smoothing_length);
  calculateBlockSize(data.grid);
  data.smoothing_length= smoothing_length; data.particle_mass= particle_mass;data.all_particles_density_updated = false;
  data.smoothing_length_2 = pow(data.smoothing_length,2);data.smoothing_length_6 = pow(smoothing_length, 6); data.smoothing_length_9 = pow(smoothing_length, 9);
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
    std::cout << "check0" << '\n';
    int real_particles = 0;
    input_file.seekg(8,std::ios::beg);
    while(!input_file.eof()){
        Particle particle;
        float floatField;
        for (auto& field :{&particle.pos[0], &particle.pos[1], &particle.pos[2], &particle.hv[0], &particle.hv[1], &particle.hv[2],&particle.vel[0], &particle.vel[1], &particle.vel[2]}){
          input_file.read(reinterpret_cast<char*>(&floatField), sizeof(float));
          *field = static_cast<double>(floatField);
        }
        //Calculamos los indices
        std::vector<int> particle_block_index = calcParticleIndex(particle, data.grid);

        checkBlockIndex(particle_block_index[0], particle_block_index[1], particle_block_index[2], data.grid);

        if (input_file.eof()) {
          break;
        }
        particle.id = real_particles;
        for (auto& block : data.grid.grid_blocks){
            if (block.second.block_index == particle_block_index){
                block.second.block_particles.push_back(particle);
            }
        }

        ++real_particles;
        }
    return real_particles;
}

int initiateSimulation(const std::string& n_iterations, const std::string& inputFile) {
    auto start = std::chrono::high_resolution_clock::now();
    std::ifstream input_file = openFile(inputFile);
    if(!input_file){throwException("Cannot open " + inputFile + " for reading", -3);}
    int n_iterations_int = std::stoi(n_iterations);
    // double doubleNumber = static_cast<double>(floatNumber);
    float ppmFloat;
    input_file.read(reinterpret_cast<char *>(&ppmFloat), sizeof(float));
    double ppm = static_cast<double>(ppmFloat);
    int np;
    input_file.read(reinterpret_cast<char *>(&np), sizeof(int));

    if (np < 0) {
        throwException("Error: Invalid number of particles: " + std::to_string(np), -5);
    }
    Grid grid;
    SimulationData data(grid);
    calculateParameters(ppm, np, data);
    // Debemos llamar a la función que guarda los parámetros de las partículas
    int real_particles;
    real_particles = setParticleData(inputFile, data);

    if (np != real_particles) {
        throwException("Error: Number of particles mismatch. Header:  " + std::to_string(np) +
                           ", Found: " + std::to_string(real_particles),
                       -5);
    }
    //int c = 0;
    for (int i = 0; i < n_iterations_int; ++i) {
        processSimulation(data);
        //std::cout << "----" << '\n';
        //std::cout << c << '\n';
        //std::cout << "----" << '\n';
        //c += 1;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << duration.count()/1000000 << '\n';
    return 0;
  }
