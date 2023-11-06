//


#ifndef FLUID_PROCESSSIMULATION_HPP
#define FLUID_PROCESSSIMULATION_HPP
#include "particle.hpp"
#include "grid.hpp"
#include "initSimulation.hpp"
#include <map>
#include <unordered_map>
#include <vector>

std::map<std::vector<int>,std::vector<std::vector<int>>> createAdjacentBlocks(gridSize grid);
std::map<std::vector<int>, std::vector<Particle>> createSubMap(std::map<std::vector<int>, std::vector<Particle>> particleMap, std::vector<std::vector<int>>adjacent_blocks);
double transferAcceleration(double acceleration, double dist, SimulationData data);
double transformDensity(double density, SimulationData data);
int increaseDensity(std::unordered_map<int, Particle> particleMap);

#endif  // FLUID_PROCESSSIMULATION_HPP
