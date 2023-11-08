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

bool blockExists(int i, int j, int k, gridSize grid);
std::map<std::vector<int>,std::vector<std::vector<int>>> createAdjacentBlocks(gridSize grid);
std::map<std::vector<int>, std::vector<Particle>> createSubMap(std::map<std::vector<int>, std::vector<Particle>> particleMap, std::vector<std::vector<int>>adjacent_blocks);
double calculateNorm(std::vector<double>particlei,std::vector<double>particlej);
double transformDensity(double density, SimulationData data);
double transferAcceleration(double acceleration, double dist, SimulationData data);
int updateDensityFalse(std::map<std::vector<int>, std::vector<Particle>> particleMap);
void updateBlock(std::vector<Particle> current_block_particles, std::map<std::vector<int>, std::vector<Particle>> particleSubMap, SimulationData data, std::string mode);
void updateBlock2(std::vector<Particle> current_block_particles, std::map<std::vector<int>, std::vector<Particle>> particleSubMap, SimulationData data, std::string mode);
int modifyBlock(std::vector<int> current_block_key, std::map<std::vector<int>, std::vector<Particle>> particleMap, std::map<std::vector<int>,std::vector<std::vector<int>>> adjacent_blocks_map, SimulationData data);
double calcCord(Particle particle, int index);
double calcVariation(Particle particle, double cordParticle, double ngrid, int index);
double calcAcceleration(Particle particle, double var, double ngrid, int index);
void calcCollision(std::map<std::vector<int>, std::vector<Particle>> particleMap, gridSize grid);
double updatePosition(Particle particle, int index);
double updateVelocity(Particle particle, int index);
double updateHv(Particle particle, int index);
void calcParticleMovement(std::map<std::vector<int>, std::vector<Particle>> particleMap);
int processSimulation(std::map<std::vector<int>, std::vector<Particle>> particleMap, SimulationData data);

#endif  // FLUID_PROCESSSIMULATION_HPP
