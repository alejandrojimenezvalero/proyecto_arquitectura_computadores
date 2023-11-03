#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include "exceptionHandler.hpp"
#include "constants.hpp"
#include "grid.hpp"
#include "block.hpp"

using namespace std;
using namespace simulationConstants;

void calculateParameters(double ppm, int np) {
  double smoothing_length = RADIO_MULTIPLICATOR / ppm;
  double particle_mass = FLUID_DENSITY / pow(ppm, 3);

  gridSize grid = calculateGridSize(smoothing_length);
  blockSize block = calculateBlockSize(grid);

  cout << "Number of particles: " << np << '\n';
  cout << "Particles per meter: " << ppm << '\n';
  cout << "Smoothing length: " << smoothing_length << '\n';
  cout << "Particle mass: " << particle_mass << '\n';
  cout << "Grid size: " << grid.nx << " x " << grid.ny << " x " << grid.nz << '\n';
  cout << "Number of blocks: " << grid.nx  * grid.ny  * grid.nz  << '\n';
  cout << "Block size: " << block.sx << " x " << block.sy << " x " << block.sz << '\n';
}

std::ifstream openFile(const std::string& inputFile){
  std::ifstream input_file(inputFile);
  exceptionHandler(!input_file, "Cannot open " + inputFile + " for reading", -3);
  return input_file;
}
int saveParticles(const std::string& inputFile){
  std::ifstream input_file = openFile(inputFile);
  std::vector<double> buffer;
  float value;
  while (input_file.read(reinterpret_cast<char *>(&value), sizeof(float))) {
    buffer.push_back(value);
  }

  size_t buffer_size = buffer.size();
  int bsize = static_cast<int>(buffer_size);

  int real_particles_number = bsize / 9;
  return real_particles_number;
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

  calculateParameters(ppm, np);
  //Debemos llamar a la función que guarda los parámetros de las partículas
  int real_particles_number = saveParticles(inputFile);
  exceptionHandler(np != real_particles_number, "Error: Number of particles mismatch. Header:  "+ std::to_string(np) + ", Found: " + std::to_string(real_particles_number)  , -5);
  return 0;
}

int validateParameters(const std::vector<std::string> &arguments) {
  size_t n_args = arguments.size();
  exceptionHandler(n_args != 3,"Error: Invalid number of arguments: " + std::to_string(n_args) + ".", -1);
  ifstream inputFile(arguments[1]); ofstream outputFile(arguments[2]);
  try {
    exceptionHandler(stoi(arguments[0]) < 0, "Error: Invalid number of time steps.", -2);
  } catch (const invalid_argument &) {throwException("Error: time steps must be numeric.", -1);}
  exceptionHandler(!inputFile, "Cannot open " + arguments[1] + " for reading", -3);
  exceptionHandler(!outputFile, "Cannot open " + arguments[2] + " for writing", -4);
  inputFile.close();
  outputFile.close();
  return 0;
}

//Mover las funciones que no sean de validacion de argumentos de entrada del main a otros ficheros
//Comenzar a leer los datos de la particulas: particles sera un struct
//Calcular el los indices de bloque para cada particula
//Refactorizar la lectura de particulas como struct of arrays
