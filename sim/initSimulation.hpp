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

struct SimulationData{
    Grid& grid;
    double smoothing_length;
    double smoothing_length_2;
    double smoothing_length_6;
    double smoothing_length_9;
    double particle_mass;
    double ppm;
    double np;
    bool all_particles_density_updated;
    double escalar_pos;
    double escalar_vel;
    double escalar_density;

    SimulationData(Grid& initialGrid)
            :grid(initialGrid), smoothing_length(0.0), particle_mass(0.0){}
};

void createGridBlocks(Grid& grid);
void initAdjIndexVectorBlocks(Grid& grid);
void checkBlockIndex(int &i, int & index_j, int & index_k, Grid& grid);
void calculateParameters(double ppm, int np, SimulationData& data);
int setParticleData(std::ifstream& input_file, SimulationData& data);
int initiateSimulation(const std::string& n_iterations, const std::string& inputFile, const std::string& outputFile);
void readParticleFields(std::ifstream& input_file, Particle& particle);
void addParticleToBlock(const Particle& particle, Grid& grid, const std::vector<int>& particle_block_index);

#endif  // FLUID_INITSIMULATION_HPP
