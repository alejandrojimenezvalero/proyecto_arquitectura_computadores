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
#include <map>
#include <tuple>
#include <memory>

using namespace simulationConstants;

void createMap(Grid& grid) {
    for (int i = 0; i < grid.ngrid[0]; ++i) {
        for (int j = 0; j < grid.ngrid[1]; ++j) {
            for (int k = 0; k < grid.ngrid[2]; ++k) {
                Block block = createBlock(i, j, k);
                grid.particleMap.push_back(block);
            }
        }
    }
}

int checkBlockIndex(int &i, int &j, int &k, Grid& grid){
  i = (i < 0)? 0:i;
  i = (i > grid.ngrid[0]-1)? grid.ngrid[0]-1:i;
  j = (j < 0)? 0:j;
  j = (j > grid.ngrid[1]-1)? grid.ngrid[1]-1:j;
  k = (k < 0)? 0:k;
  k = (k > grid.ngrid[2]-1)? grid.ngrid[2]-1:k;
  return 0;
}
void calculateParameters(double ppm, int np, SimulationData& data) {
  double smoothing_length = RADIO_MULTIPLICATOR / ppm;
  double particle_mass = FLUID_DENSITY / pow(ppm, 3);
  initGrid(data.grid, smoothing_length);
  calculateBlockSize(data.grid);
  data.smoothing_length= smoothing_length; data.particle_mass= particle_mass;data.all_particles_density_updated = false;
  std::cout << "Number of particles: " << np << '\n';
  std::cout << "Particles per meter: " << ppm << '\n';
  std::cout << "Smoothing length: " << smoothing_length << '\n';
  std::cout << "Particle mass: " << particle_mass << '\n';
  std::cout << "Grid size: " << data.grid.ngrid[0] << " x " << data.grid.ngrid[1] << " x " << data.grid.ngrid[2] << '\n';
  std::cout << "Number of blocks: " << data.grid.ngrid[0]  * data.grid.ngrid[1]  * data.grid.ngrid[2]  << '\n';
  std::cout << "Block size: " << data.grid.block_dimensions[0] << " x " << data.grid.block_dimensions[1] << " x " << data.grid.block_dimensions[2] << '\n';
}
int setParticleData(const std::string& inputFile, SimulationData& data){
    std::ifstream input_file = openFile(inputFile);
    createMap(data.grid);
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
        std::vector<int> particle_block_index = calcParticleIndex(particle, data.grid.block_dimensions);

        checkBlockIndex(particle_block_index[0], particle_block_index[1], particle_block_index[2], data.grid);

        if (input_file.eof()) {
          break;
        }
        particle.id = real_particles;
        for (Block& block : data.grid.particleMap){
            if (block.block_index == particle_block_index){
                block.block_particles.push_back(particle);
            }
        }

        ++real_particles;
        }
    return real_particles;
}

int initiateSimulation(const std::string& n_iterations, const std::string& inputFile) {
    std::ifstream input_file = openFile(inputFile);
    int n_iterations_int = std::stoi(n_iterations);
    //double doubleNumber = static_cast<double>(floatNumber);
    float ppmFloat;
    input_file.read(reinterpret_cast<char *>(&ppmFloat), sizeof(float));
    double ppm = static_cast<double>(ppmFloat);
    int np;
    input_file.read(reinterpret_cast<char *>(&np), sizeof(int));

    exceptionHandler(np < 0, "Error: Invalid number of particles: " + std::to_string(np));
    Grid grid;
    SimulationData data(grid);
    calculateParameters(ppm, np,data);
    //Debemos llamar a la función que guarda los parámetros de las partículas
    int real_particles;
    real_particles = setParticleData(inputFile, data);

    exceptionHandler(np != real_particles, "Error: Number of particles mismatch. Header:  "+ std::to_string(np) + ", Found: " + std::to_string(real_particles));
    for(int i=0; i< n_iterations_int; ++i){
      processSimulation(data);
    }
    return 0;
}
