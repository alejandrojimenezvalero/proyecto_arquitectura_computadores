//

#ifndef FLUID_INITSIMULATION_HPP
#define FLUID_INITSIMULATION_HPP

#include "block.hpp"
#include "grid.hpp"
#include "particle.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <memory>

struct SimulationData {
    Grid grid;
    double smoothing_length = 0.0;
    double smoothing_length_2 = 0.0;
    double smoothing_length_6 = 0.0;
    double smoothing_length_9 = 0.0;
    double particle_mass = 0.0;
    double ppm = 0.0;
    double np = 0.0;
    double escalar_pos = 0.0;
    double escalar_vel = 0.0;
    double escalar_density = 0.0;
};


void createGridBlocks(Grid& grid);
void initAdjIndexVectorBlocks(Grid& grid);
void checkBlockIndex(std::vector<int>& particle_block_index, Grid& grid);
void calculateParameters(double ppm, int np, SimulationData& data);
int setParticleData(std::ifstream& input_file, SimulationData& data);
int initiateSimulation(const std::vector<std::string>& stringInitVector);
void readParticleFields(std::ifstream& input_file, Particle& particle);
void addParticleToBlock(const Particle& particle, Grid& grid, const std::vector<int>& particle_block_index);

#endif  // FLUID_INITSIMULATION_HPP
