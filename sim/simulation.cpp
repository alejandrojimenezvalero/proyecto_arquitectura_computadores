//
#include "exceptionHandler.hpp"
#include "constants.hpp"
#include "grid.hpp"
#include "block.hpp"
#include "fileManager.hpp"
#include "particle.hpp"
#include "simulation.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

using namespace simulationConstants;

int checkBlockIndex(int &i, int &j, int &k, gridSize grid){
  i = (i < 0)? 0:i;
  i = (i > grid.nx-1)? grid.nx-1:i;
  j = (j < 0)? 0:j;
  j = (j > grid.ny-1)? grid.ny-1:j;
  k = (k < 0)? 0:k;
  k = (k > grid.nz-1)? grid.nz-1:k;
  return 0;
}
SimulationData calculateParameters(double ppm, int np) {
  double smoothing_length = RADIO_MULTIPLICATOR / ppm;
  double particle_mass = FLUID_DENSITY / pow(ppm, 3);

  gridSize grid = calculateGridSize(smoothing_length);
  blockSize block = calculateBlockSize(grid);

  SimulationData data{grid, block};
  std::cout << "Number of particles: " << np << '\n';
  std::cout << "Particles per meter: " << ppm << '\n';
  std::cout << "Smoothing length: " << smoothing_length << '\n';
  std::cout << "Particle mass: " << particle_mass << '\n';
  std::cout << "Grid size: " << grid.nx << " x " << grid.ny << " x " << grid.nz << '\n';
  std::cout << "Number of blocks: " << grid.nx  * grid.ny  * grid.nz  << '\n';
  std::cout << "Block size: " << block.sx << " x " << block.sy << " x " << block.sz << '\n';
  return data;
}
int setParticleData(const std::string& inputFile, SimulationData data){
  std::ifstream input_file = openFile(inputFile);
  std::unordered_map<int, Particle> particleMap;
  int particle_index = 0;
  input_file.seekg(8,std::ios::beg);
  while(!input_file.eof()){
    Particle particle;
    float floatField;
    for (auto& field :{&particle.px, &particle.py, &particle.pz, &particle.hvx, &particle.hvy, &particle.hvz,&particle.vx, &particle.vy, &particle.vz}){
      input_file.read(reinterpret_cast<char*>(&floatField), sizeof(float));
      *field = static_cast<double>(floatField);
    }
    //Calculamos los indices
    blockSize block = data.block;
    particle.i = std::floor((particle.px - LOWER_LIMIT[0])/block.sx);
    particle.j = std::floor((particle.py - LOWER_LIMIT[1])/block.sy);
    particle.k = std::floor((particle.pz - LOWER_LIMIT[2])/block.sz);

    gridSize grid = data.grid;
    checkBlockIndex(particle.i, particle.j, particle.k, grid);
    if (input_file.eof()) {
      break;
    }
    particleMap[particle_index] = particle;
    ++particle_index;
  };
  /*
  Particle particula2 = particleMap[989];
  std::cout << "Partícula con índice 1675:\n"
  << "px: " << particula2.px << "\n"
  << "py: " << particula2.py << "\n"
  << "pz: " << particula2.pz << "\n"
  << "hvx: " << particula2.hvx << "\n"
  << "hvy: " << particula2.hvy << "\n"
  << "hvz: " << particula2.hvz << "\n"
  << "vx: " << particula2.vx << "\n"
  << "vy: " << particula2.vy << "\n"
  << "vz: " << particula2.vz << "\n"
  << "i: " << particula2.i << "\n"
  << "j: " << particula2.j << "\n"
  << "k: " << particula2.k << std::endl;
  */
  return particle_index;
}

int initiateSimulation(const std::string& inputFile) {
  std::ifstream input_file = openFile(inputFile);
  //double doubleNumber = static_cast<double>(floatNumber);
  float ppmFloat;
  input_file.read(reinterpret_cast<char *>(&ppmFloat), sizeof(float));
  double ppm = static_cast<double>(ppmFloat);
  int np;
  input_file.read(reinterpret_cast<char *>(&np), sizeof(int));

  exceptionHandler(np < 0, "Error: Invalid number of particles: " + std::to_string(np) , -5);

  SimulationData data = calculateParameters(ppm, np);
  //Debemos llamar a la función que guarda los parámetros de las partículas
  int real_particles_number = setParticleData(inputFile, data);
  exceptionHandler(np != real_particles_number, "Error: Number of particles mismatch. Header:  "+ std::to_string(np) + ", Found: " + std::to_string(real_particles_number)  , -5);
  return 0;
}
