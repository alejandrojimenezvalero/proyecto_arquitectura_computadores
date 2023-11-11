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

bool blockExists(int i, int j, int k, gridSize& grid);
std::map<std::vector<int>,std::vector<std::vector<int>>> createAdjacentBlocks(gridSize& grid);
std::tuple<std::shared_ptr<std::vector<Particle>>, std::vector<Block> > createSubMap(std::vector<int> current_block_key, std::vector<Block>& particleMap, std::vector<std::vector<int>> adjacent_blocks);
double calculateNorm(const std::vector<double>& particlei, const std::vector<double>& particlej);
double transformDensity(double density, const SimulationData& data);
std::vector<double> transferAcceleration(Particle& particlei, Particle& particlej, double dist, const SimulationData& data);
//void updateBlock(std::vector<Particle> current_block_particles, std::map<std::vector<int>, std::vector<Particle>> particleSubMap, SimulationData data, std::string mode);
void initializeDensityAcceleration(std::map<std::vector<int>, std::shared_ptr<std::vector<Particle>>>& particleMap);
void updateBlock2(std::shared_ptr<std::vector<Particle>> current_block_particles, std::map<std::vector<int>, std::shared_ptr<std::vector<Particle>>>& particleSubMap, SimulationData& data, std::string mode);
int modifyBlock(std::vector<int> current_block_key, std::map<std::vector<int>, std::shared_ptr<std::vector<Particle>>>& particleMap, std::map<std::vector<int>,std::vector<std::vector<int>>>& adjacent_blocks_map, SimulationData& data);
double calcCord(Particle& particle, int index);
double calcVariation(int index_block, double cordParticle, double ngrid, int index);
double calcAcceleration(Particle& particle, double var, double ngrid, int index);
//void calcCollision(std::map<std::vector<int>, std::shared_ptr<std::vector<Particle>>>& particleMap, gridSize& grid);
double updatePosition(Particle& particle, int index);
double updateVelocity(Particle& particle, int index);
double updateHv(Particle& particle, int index);
void removeParticlesFromBlock(std::shared_ptr<std::vector<Particle>>& old_block_particles, std::vector<Particle>& particles_to_remove);
void updateParticleBlockBelonging(std::vector<int>old_block_index, Particle& particle, std::map<std::vector<int>, std::shared_ptr<std::vector<Particle>>>& particleMap, SimulationData& data);
void updateParticle(std::vector<int>block_index, Particle& particle, SimulationData data);
void establishParticleFunctionality(std::map<std::vector<int>, std::shared_ptr<std::vector<Particle>>>& particleMap, gridSize& grid);
int processSimulation(std::vector<Block>& particleMap, SimulationData& data);

#endif  // FLUID_PROCESSSIMULATION_HPP
