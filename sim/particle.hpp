//

#ifndef FLUID_PARTICLE_HPP
#define FLUID_PARTICLE_HPP

struct Particle {
    double px, py, pz; // Coordenadas de la posición
    double hvx, hvy, hvz; // Coordenadas del vector hv
    double vx, vy, vz; // Coordenadas de la velocidad
    int i, j, k = 0 ; //Indices de bloque para la particula
};


#endif  // FLUID_PARTICLE_HPP
