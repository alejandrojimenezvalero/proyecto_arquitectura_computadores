//
#include "initSimulation.hpp"

#include "block.hpp"
#include "constants.hpp"
#include "exceptionHandler.hpp"
#include "fileManager.hpp"
#include "grid.hpp"
#include "particle.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <map>

using namespace simulationConstants;


std::map<std::vector<int>, std::vector<Particle>> createMap(gridSize grid){

  std::map<std::vector<int>, std::vector<Particle>> initMap;
  for (int i =0; i < grid.nx; ++i) {
    for (int j=0; j < grid.ny; ++j) {
      for (int k=0; k < grid.nz; ++k){
        std::vector<Particle> addVector{};
        std::vector<int> key = {i, j, k};
        initMap[key] = addVector;
      }
    }
  }
  return initMap;

}
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

  SimulationData data{grid, block, smoothing_length, particle_mass};
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
  //std::map<std::vector<int>, std::unordered_map<int, Particle>> particleMap;
  std::map<std::vector<int>, std::vector<Particle>> particleMap = createMap(data.grid);
  //std::cout << "Particle Map: " << particleMap.size() << '\n';
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
    blockSize block = data.block;
    particle.i = std::floor((particle.px - LOWER_LIMIT[0])/block.sx);
    particle.j = std::floor((particle.py - LOWER_LIMIT[1])/block.sy);
    particle.k = std::floor((particle.pz - LOWER_LIMIT[2])/block.sz);

    gridSize grid = data.grid;
    checkBlockIndex(particle.i, particle.j, particle.k, grid);
    //std::cout << "Particle index:" << particle.i << " " << particle.j << " " << particle.k << " " << '\n';
    if (input_file.eof()) {
      break;
    }
    particleMap[{particle.i,particle.j,particle.k}].push_back(particle);
    ++real_particles;
  }
  /* COMPROBAMOS QUE EL MAPA ESTÁ BIEN CREADO, VIENDO QUE EN CADA BLOQUE HAY X PARTICULAS Y VIENDO QUE EL NÚ,ERO TOTAL ES EXACTAMENTE 4800
  int contador = 0;
  for(auto &keyValue: particleMap){
    std::vector<int> key = keyValue.first;
    std::vector<Particle> particles = keyValue.second;
    //std::cout << "Bloque de indices " << "i: "<< key[0] << " j: "<< key[1] << " k: " << key[2] << '\n' << " => Particulas del bloque: "<< '\n';
    //std::cout << "Número de partículas de bloque: " << particles.size() << '\n';
    contador += particles.size();
    for (Particle particle : particles) {
      std::cout << "Partícula: " << particle.px << "," << '\n';
    }
  }
  std::cout << "Numero de particulas total: " << contador << '\n';
  */
  /*
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
  return real_particles;
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
  int real_particles= setParticleData(inputFile, data);
  exceptionHandler(np != real_particles, "Error: Number of particles mismatch. Header:  "+ std::to_string(np) + ", Found: " + std::to_string(real_particles)  , -5);
  return 0;
}