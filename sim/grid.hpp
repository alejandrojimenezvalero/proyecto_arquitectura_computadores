#ifndef GRID_HPP
#define GRID_HPP


struct gridSize{
    double nx;
    double ny;
    double nz;
};

gridSize calculateGridSize(double smoothing_length);
#endif