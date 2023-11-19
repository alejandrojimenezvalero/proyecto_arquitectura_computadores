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

bool blockExists(int i, int j, int k, Grid& grid);
//std::map<std::vector<int>,std::vector<std::vector<int>>> createAdjacentBlocks(Grid& grid);
double calculateNorm(const std::vector<double>& particlei, const std::vector<double>& particlej);
void transformDensity(Particle& particle, const SimulationData& data);
void transferAcceleration(Particle& particlei, Particle& particlej, double& dist, const SimulationData& data);
//void updateBlock(std::vector<Particle> current_block_particles, std::map<std::vector<int>, std::vector<Particle>> particleSubMap, SimulationData data, std::string mode);
//void initializeDensityAcceleration(SimulationData& data);
void updateBlocksDensity(SimulationData& data);
void updateBlocksAcceleration(SimulationData& data);
//void updateBlock2(SubGrid& particleSubMap , SimulationData& data, std::string mode);
//int modifyBlock(Block& current_block, std::map<std::vector<int>, std::vector<std::vector<int>>>& adjacent_blocks_map, SimulationData& data);
double calcCord(Particle& particle, int index);
double calcVariation(int& index_block, double& cordParticle, const SimulationData& data, int& index);
double calcAcceleration(Particle& particle, double& var, SimulationData& data, int& index);
//void calcCollision(std::map<std::vector<int>, std::shared_ptr<std::vector<Particle>>>& grid_blocks, Grid& grid);
double updatePosition(Particle& particle, int& index);
double updateVelocity(Particle& particle, int& index);
double updateHv(Particle& particle, int& index);
void removeParticlesFromBlock(Block& block, std::vector<Particle>& particles_to_remove);
void updateParticleBlockBelonging(SimulationData& data);
void updateParticle(std::vector<int>& block_index, Particle& particle, SimulationData& data);
void establishParticleFunctionality(SimulationData& data);
int processSimulation(SimulationData& data);

#endif  // FLUID_PROCESSSIMULATION_HPP
