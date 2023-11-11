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

std::vector<Block> createMap(SimulationData data) {
    std::vector<Block> initMap;
    for (int i = 0; i < data.grid.nx; ++i) {
        for (int j = 0; j < data.grid.ny; ++j) {
            for (int k = 0; k < data.grid.nz; ++k) {
                Block block = createBlock(i, j, k);
                block.block_particles = initializeBlockParticles();
                initMap.emplace_back(block);
            }
        }
    }
    return initMap;
}


int checkBlockIndex(int &i, int &j, int &k, SimulationData data){
  i = (i < 0)? 0:i;
  i = (i > data.grid.nx-1)? data.grid.nx-1:i;
  j = (j < 0)? 0:j;
  j = (j > data.grid.ny-1)? data.grid.ny-1:j;
  k = (k < 0)? 0:k;
  k = (k > data.grid.nz-1)? data.grid.nz-1:k;
  return 0;
}
SimulationData calculateParameters(double ppm, int np) {
  double smoothing_length = RADIO_MULTIPLICATOR / ppm;
  double particle_mass = FLUID_DENSITY / pow(ppm, 3);

  gridSize grid = calculateGridSize(smoothing_length);
  std::vector<double> block_dimensions = calculateBlockSize(grid);
  grid.block_dimensions = block_dimensions;
  SimulationData data{grid, smoothing_length, particle_mass, false};
  std::cout << "Number of particles: " << np << '\n';
  std::cout << "Particles per meter: " << ppm << '\n';
  std::cout << "Smoothing length: " << smoothing_length << '\n';
  std::cout << "Particle mass: " << particle_mass << '\n';
  std::cout << "Grid size: " << grid.nx << " x " << grid.ny << " x " << grid.nz << '\n';
  std::cout << "Number of blocks: " << grid.nx  * grid.ny  * grid.nz  << '\n';
  std::cout << "Block size: " << block_dimensions[0] << " x " << block_dimensions[1] << " x " << block_dimensions[2] << '\n';
  return data;
}
std::tuple< int, std::vector<Block> > setParticleData(const std::string& inputFile, SimulationData data){
    std::ifstream input_file = openFile(inputFile);
    std::vector<Block> particleMap = createMap(data);
    int real_particles = 0;
    input_file.seekg(8,std::ios::beg);
    while(!input_file.eof()){
        Particle particle;
        float floatField;
        for (auto& field :{&particle.px, &particle.py, &particle.pz, &particle.hvx, &particle.hvy, &particle.hvz,&particle.vx, &particle.vy, &particle.vz}){
          input_file.read(reinterpret_cast<char*>(&floatField), sizeof(float));
          *field = static_cast<double>(floatField);
        }
        //Calculamos los indices
        std::vector<int> particle_block_index = calcParticleIndex(particle, data.grid.block_dimensions);

        checkBlockIndex(particle_block_index[0], particle_block_index[1], particle_block_index[2], data);
        //std::cout << "Particle index:" << particle.i << " " << particle.j << " " << particle.k << " " << '\n';
        if (input_file.eof()) {
          break;
        }
        particle.id = real_particles;
        for (Block& block:particleMap){
            if (block.block_index == particle_block_index){
                block.block_particles->push_back(particle);
            }
        }

        ++real_particles;
        }
    return std::make_tuple(real_particles, particleMap);
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

    SimulationData data = calculateParameters(ppm, np);
    //Debemos llamar a la función que guarda los parámetros de las partículas
    int real_particles;
    std::vector<Block> particleMap;
    std::tie(real_particles, particleMap) = setParticleData(inputFile, data);

    exceptionHandler(np != real_particles, "Error: Number of particles mismatch. Header:  "+ std::to_string(np) + ", Found: " + std::to_string(real_particles));
    for(int i=0; i< n_iterations_int; ++i){
      processSimulation(particleMap,data);
    }
    return 0;
}
