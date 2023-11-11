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
    bool acceleration_updated = false;
    std::vector<double> acceleration = GRAVITY;
    Particle() : id(0), px(0.0), py(0.0), pz(0.0), hvx(0.0), hvy(0.0), hvz(0.0),
                 vx(0.0), vy(0.0), vz(0.0), i(0), j(0), k(0), density(0.0),
                 density_updated(false), acceleration_updated(false),
                 acceleration(GRAVITY) {}
    bool operator==(const Particle& other) const {
        return (id == other.id);
    }
};


#endif  // FLUID_PARTICLE_HPP
