#ifndef GRID_HPP
#define GRID_HPP

#include "particle.hpp"

#include <vector>
#include <map>

struct gridSize{
    double nx;
    double ny;
    double nz;
};

gridSize calculateGridSize(double smoothing_length);
#endif