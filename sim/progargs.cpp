#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdexcept>

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
int initiateSimulation(std::ifstream &input_file) {
  int res = 0;
  if (!input_file) {
    res = -3;
    return res;
  }

  float ppm;
  int np;
  input_file.read(reinterpret_cast<char *>(&ppm), sizeof(float));
  input_file.read(reinterpret_cast<char *>(&np), sizeof(int));

  if (np < 0) {
    std::cerr << "Error: Invalid number of particles: " << np << "." << '\n';
    res = -5;
    return res;
  }

  std::vector<float> buffer;
  float value;
  while (input_file.read(reinterpret_cast<char *>(&value), sizeof(float))) {
    buffer.push_back(value);
  }

  size_t buffer_size = buffer.size();
  int bsize = static_cast<int>(buffer_size);

  int real_particles_number = bsize / 9;

  if (np != real_particles_number) {
    std::cerr << "Error: Number of particles mismatch. Header:  " << np
              << ", Found: " << real_particles_number << '\n';
    res = -5;
    return res;
  }

  double smoothing_length = simulationConstants::radio_multiplicator / ppm;
  double particle_mass = simulationConstants::fluid_density / pow(ppm, 3);

  std::cout << "Number of particles: " << np << std::endl;
  std::cout << "Particles per meter: " << ppm << std::endl;
  std::cout << "Smoothing length: " << smoothing_length << std::endl;
  std::cout << "Particle mass: " << particle_mass << std::endl;

  double x_length = simulationConstants::upper_limit[0] - simulationConstants::lower_limit[0];
  double y_length = simulationConstants::upper_limit[1] - simulationConstants::lower_limit[1];
  double z_length = simulationConstants::upper_limit[2] - simulationConstants::lower_limit[2];

  double nx = round((x_length) / smoothing_length);
  double ny = round((y_length) / smoothing_length);
  double nz = round((z_length) / smoothing_length);

  double sx = (x_length) / nx;
  double sy = (y_length) / ny;
  double sz = (z_length) / nz;

  std::cout << "Grid size: " << nx << " x " << ny << " x " << nz << std::endl;
  std::cout << "Number of blocks: " << nx * ny * nz << std::endl;
  std::cout << "Block size: " << sx << " x " << sy << " x " << sz << std::endl;

  // No es necesario cerrar el archivo aquí, ya que se cierra automáticamente al salir de la función
  return res;
}


int ProgArgs(const std::vector<std::string> &arguments) {
  int res = 0;
    const size_t n_args = arguments.size();
  if (n_args != 3) {
    std::cerr << "Error: Invalid number of arguments: " << n_args << "." << '\n';
    return res = -1;
  }

    const std::string inputFile = arguments[1];
    const std::string outputFile = arguments[2];

  std::ifstream input_file(inputFile, std::ios::binary);
  std::ofstream output_file(outputFile);

  if (!input_file) {
    std::cerr << "Error: Could not open the input file." << '\n';
    return -3;
  }

  if (!output_file) {
    std::cerr << "Error: Could not open the output file." << '\n';
    return -4;
  }

  try {
        const int numSteps = std::stoi(arguments[0]);
    if (numSteps < 0) {
      std::cerr << "Error: Invalid number of time steps." << '\n';
      return -2;
    }
  } catch (const std::invalid_argument &) {
    std::cerr << "Error: Time steps must be numeric." << '\n';
    return -1;
  }

  res = initiateSimulation(input_file);

  if (res < 0) {
    return res;
  }

  // Cierra los archivos al salir de la función
  input_file.close();
  output_file.close();

  return res;
}
