//

#ifndef FLUID_PARTICLE_HPP
#define FLUID_PARTICLE_HPP
#include "constants.hpp"

using namespace simulationConstants;

struct Particle {
    int id;
    double px, py, pz; // Coordenadas de la posición
    double hvx, hvy, hvz; // Coordenadas del vector hv
    double vx, vy, vz; // Coordenadas de la velocidad
    int i, j, k = 0; //Indices de bloque para la particula
    double density = 0;
    bool density_updated = false;
    bool density_updated2 = false;
    bool acceleration_updated = false;
    std::vector<double> acceleration = GRAVITY;
    bool operator==(const Particle& other) const {
        return (id == other.id);
    }
};


#endif  // FLUID_PARTICLE_HPP
