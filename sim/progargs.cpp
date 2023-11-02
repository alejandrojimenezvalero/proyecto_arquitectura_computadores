#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include "../sim/exceptionHandler.cpp"
using namespace std;
namespace simulationConstants {
  const float radio_multiplicator = 1.695;
  const int fluid_density = 1000;
  /*
  const float stiffness_pressure = 3.0;
  const int stiffness_colisions = 30000;
  const float damping = 128.0;
  const float viscosity = 0.4;
  const float particle_size = 0.0002;
  const float timestamp = 0.001;
   */
  const vector<float> upper_limit{0.065, 0.1, 0.065};
  const vector<float> lower_limit{-0.065,-0.08, -0.065};
}
void calculateParameters(double ppm, int np) {
  using namespace simulationConstants;
  double smoothing_length = radio_multiplicator / ppm;
  double particle_mass = fluid_density / pow(ppm, 3);

  vector<double> enclousure {upper_limit[0] - lower_limit[0], upper_limit[1] - lower_limit[1],upper_limit[2] - lower_limit[2]};
  vector<double> mesh {round((enclousure[0]) / smoothing_length), round((enclousure[1]) / smoothing_length), round((enclousure[2]) / smoothing_length)};
  vector<double> block_size{enclousure[0]/mesh[0],enclousure[1]/mesh[1],enclousure[2]/mesh[2]};

  cout << "Number of particles: " << np << '\n';
  cout << "Particles per meter: " << ppm << '\n';
  cout << "Smoothing length: " << smoothing_length << '\n';
  cout << "Particle mass: " << particle_mass << '\n';
  cout << "Grid size: " << mesh[0] << " x " << mesh[1] << " x " << mesh[2] << '\n';
  cout << "Number of blocks: " << mesh[0]  * mesh[1]  * mesh[2]  << '\n';
  cout << "Block size: " << block_size[0] << " x " << block_size[1] << " x " << block_size[2] << '\n';
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


//Definir que variables son constantes
//Refactorizar la lectura de particulas como struct of arrays
/*Decidir que funciones van en los difernetes modulos sobre todo lo referido a grid y block*/