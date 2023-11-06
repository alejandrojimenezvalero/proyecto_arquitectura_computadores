//

#ifndef FLUID_PARTICLE_HPP
#define FLUID_PARTICLE_HPP
#include "constants.hpp"

using namespace simulationConstants;

struct Particle {
    double px, py, pz; // Coordenadas de la posición
    double hvx, hvy, hvz; // Coordenadas del vector hv
    double vx, vy, vz; // Coordenadas de la velocidad
    int i, j, k = 0 ; //Indices de bloque para la particula
    double density = 0;
    bool density_updated = false;
    bool acceleration_updated = false;
    std::vector<float> acceleration = GRAVITY;

};


#endif  // FLUID_PARTICLE_HPP
