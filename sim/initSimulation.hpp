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

void createGridBlocks(Grid &grid);

void initAdjIndexVectorBlocks(Grid &grid);

void checkBlockIndex(std::vector<int> &particle_block_index, Grid &grid);

void initializeData(SimulationData &data, double smoothing_length, double particle_mass);

void calculateParameters(double ppm, int np, SimulationData &data);

void readParticleFields(std::ifstream &input_file, Particle &particle);

int setParticleData(std::ifstream &input_file, SimulationData &data);

int initiateSimulation(const std::vector<std::string> &stringInitVector);

#endif  // FLUID_INITSIMULATION_HPP
