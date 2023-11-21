

#ifndef FLUID_SIMULATION_RESULTS_HPP
#define FLUID_SIMULATION_RESULTS_HPP

#include "particle.hpp"
#include "initSimulation.hpp"
#include <vector>


void sortParticlesById(std::vector<Particle> &particles);

void writeFloatToStream(float value, std::ofstream &stream);

void saveData(const Particle &particle, std::ofstream &outputFile);

void writeParticleData(const std::string &filename, SimulationData &data);

#endif  // FLUID_SIMULATION_RESULTS_HPP
