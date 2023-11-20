#include "block.hpp"


using namespace simulationConstants;

Block createBlock(int i, int j, int k) {
    Block block {{i, j, k},{}, {}, {},{}};  // Inicializar con un puntero nulo
    return block;
}
