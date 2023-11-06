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


struct SimulationData{
    gridSize grid;
    blockSize block;
    double smoothing_length;
    double particle_mass;
};
std::map<std::vector<int>, std::vector<Particle>> createMap(gridSize grid);
int checkBlockIndex(int &i, int &j, int &k, gridSize grid);
SimulationData calculateParameters(double ppm, int np);
int setParticleData(const std::string& inputFile, SimulationData data);
int initiateSimulation(const std::string& inputFile);



#endif  // FLUID_INITSIMULATION_HPP
