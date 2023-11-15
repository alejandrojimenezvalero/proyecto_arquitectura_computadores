
#include "constants.hpp"
#include "block.hpp"
#include "particle.hpp"

#include <memory>

using namespace simulationConstants;

Block createBlock(int i, int j, int k) {
    Block block {{i, j, k},{}, {}};  // Inicializar con un puntero nulo
    return block;
}


#include <cmath>
#include <vector>


