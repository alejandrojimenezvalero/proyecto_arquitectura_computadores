#ifndef BLOCK_HPP
#define BLOCK_HPP

#include "grid.hpp"

struct blockSize{
    double sx;
    double sy;
    double sz;
};

blockSize calculateBlockSize(gridSize grid);

#endif