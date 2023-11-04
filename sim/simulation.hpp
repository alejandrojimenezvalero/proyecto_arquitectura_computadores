//

#ifndef FLUID_SIMULATION_HPP
#define FLUID_SIMULATION_HPP

#include "block.hpp"
#include "grid.hpp"

#include <iostream>
#include <fstream>
#include <vector>


struct SimulationData{
    gridSize grid;
    blockSize block;
};
int checkBlockIndex(int &i, int &j, int &k, SimulationData data);
SimulationData calculateParameters(double ppm, int np);
int setParticleData(const std::string& inputFile, SimulationData data);
int initiateSimulation(const std::string& inputFile);




#endif  // FLUID_SIMULATION_HPP
