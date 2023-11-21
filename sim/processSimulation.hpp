//
#ifndef FLUID_PROCESSSIMULATION_HPP
#define FLUID_PROCESSSIMULATION_HPP

#include "particle.hpp"
#include "grid.hpp"
#include "initSimulation.hpp"
#include <map>
#include <unordered_map>
#include <vector>
#include <string>

double calculateNorm(const std::vector<double> &particlei, const std::vector<double> &particlej);

void transformDensity(Particle &particle, const double &smoothing_length_6, const double &escalar_density);

void transferAcceleration(Particle &particlei, Particle &particlej, const double &dist, const SimulationData &data);

void updateBlocksDensity(SimulationData &data);

void updateBlocksAcceleration(SimulationData &data);

double calcCord(Particle &particle, int index);

double calcVariation(int &index_block, const double &cordParticle, const double &grid_dimensions, int &index);

void calcAcceleration(Particle &particle, const double &var, SimulationData &data, int &index);

double updatePosition(Particle &particle, int &index);

double updateVelocity(Particle &particle, int &index);

double updateHv(Particle &particle, int &index);

void checkBorderLimits(Particle &particle, const double &grid_dimensions, int index, int actual_block_index);

void updateParticleBlockBelonging(SimulationData &data);

void updateParticle(std::vector<int> &block_index, Particle &particle, SimulationData &data);

void establishParticleFunctionality(SimulationData &data);

int processSimulation(SimulationData &data);

#endif  // FLUID_PROCESSSIMULATION_HPP
