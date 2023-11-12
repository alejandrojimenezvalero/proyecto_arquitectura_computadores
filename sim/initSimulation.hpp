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
    double particle_mass;
    bool all_particles_density_updated;

    SimulationData(Grid& initialGrid)
            :grid(initialGrid), smoothing_length(0.0), particle_mass(0.0), all_particles_density_updated(false) {}
};

void createMap(Grid& grid);
int checkBlockIndex(int &i, int &j, int &k, Grid& grid);
void calculateParameters(double ppm, int np, SimulationData& data);
int setParticleData(const std::string& inputFile, SimulationData& data);
int initiateSimulation(const std::string& n_iterations, const std::string& inputFile);

#endif  // FLUID_INITSIMULATION_HPP
