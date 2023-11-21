
#include "simulationResults.hpp"
#include "initSimulation.hpp"
#include "grid.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

void sortParticlesById(std::vector<Particle> &particles) {
    std::sort(particles.begin(), particles.end(), [](const Particle &a, const Particle &b) {
        return a.id < b.id;
    });
}


void writeFloatToStream(float value, std::ofstream& stream) {
    //NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    stream.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

void saveData(const Particle &particle, std::ofstream &outputFile) {
    for (int i = 0; i < 3; ++i) {
        writeFloatToStream(static_cast<float>(particle.pos[i]), outputFile);
    }
    for (int i = 0; i < 3; ++i) {
        writeFloatToStream(static_cast<float>(particle.hv[i]), outputFile);
    }
    for (int i = 0; i < 3; ++i) {
        writeFloatToStream(static_cast<float>(particle.vel[i]), outputFile);
    }
}

void writeParticleData(const std::string &filename, SimulationData &data) {
    std::ofstream outputFile(filename, std::ios::out | std::ios::binary);
    //En el archivo se deben escribir primero los datos ppm y np
    auto particles_per_meter = static_cast<float>(data.ppm);
    outputFile.write(
            reinterpret_cast<const char *>(&particles_per_meter),//NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
            sizeof(particles_per_meter));
    auto num_particles = static_cast<int>(data.np);
    outputFile.write(
            reinterpret_cast<const char *>(&num_particles),//NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
            sizeof(num_particles));
    //Debemos para cada particle y en order de su id imprimir pos, v, hv
    //Las particles las tenemos stored por bloques, por lo que almacenamos todas las particulas en un vector
    std::vector<Particle> grid_particles;
    for (const auto &block: data.grid.grid_blocks) {
        for (const Particle &particle: block.block_particles_v1) {
            grid_particles.push_back(particle);
        }
    }
    sortParticlesById(grid_particles);
    for (const auto &particle: grid_particles) {
        saveData(particle, outputFile);
    }
    outputFile.close();
}
