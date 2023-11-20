//

#ifndef FLUID_PARTICLE_HPP
#define FLUID_PARTICLE_HPP
#include "constants.hpp"

using namespace simulationConstants;

struct Particle {
    int id;
    std::vector<double> pos;
    std::vector<double> hv;
    std::vector<double> vel;
    double density = 0;
    bool density_updated = false;
    bool acceleration_updated = false;
    bool already_checked = false;
    std::vector<double> acceleration = GRAVITY;
    Particle() : id(0), pos{0.0,0.0,0.0}, hv{0.0,0.0,0.0}, vel{0.0,0.0,0.0},
                 density(0.0), density_updated(false), acceleration_updated(false),
                 acceleration(GRAVITY) {}
    bool operator==(const Particle& other) const {
        return (id == other.id);
    }
};


#endif  // FLUID_PARTICLE_HPP
