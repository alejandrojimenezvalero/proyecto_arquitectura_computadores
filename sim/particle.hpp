//

#ifndef FLUID_PARTICLE_HPP
#define FLUID_PARTICLE_HPP

#include "constants.hpp"

using namespace simulationConstants;

struct Particle {
    int id = 0;
    std::vector<double> pos;
    std::vector<double> hv;
    std::vector<double> vel;
    double density = 0.0;
    bool density_updated = false;
    bool acceleration_updated = false;
    std::vector<double> acceleration = GRAVITY;

    Particle() : pos{0.0, 0.0, 0.0}, hv{0.0, 0.0, 0.0}, vel{0.0, 0.0, 0.0} {}

    //boool operatoor
    bool operator==(const Particle &other) const {
        return (id == other.id);
    }
};


#endif  // FLUID_PARTICLE_HPP
