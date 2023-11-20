
#include "simulationResults.hpp"
#include "exceptionHandler.hpp"
#include "initSimulation.hpp"
#include "grid.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

void sortParticlesById(std::vector<Particle>& particles) {
  std::sort(particles.begin(), particles.end(), [](const Particle& a, const Particle& b) {
    return a.id < b.id;
  });
}
void writeParticleData(const std::string& filename, SimulationData data) {
    std::ofstream outputFile(filename, std::ios::out | std::ios::binary);
    //En el archivo se deben escribir primero los datos ppm y np
    auto particles_per_meter = static_cast<float>(data.ppm);
    outputFile.write(reinterpret_cast<const char*>(&particles_per_meter), sizeof(particles_per_meter)); //NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
    auto num_particles = static_cast<int>(data.np);
    outputFile.write(reinterpret_cast<const char*>(&num_particles), sizeof(num_particles)); //NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
    //Debemos para cada particle y en order de su id imprimir pos, v, hv
    //Las particles las tenemos stored por bloques, por lo que almacenamos todas las particulas en un vector
    std::vector<Particle> grid_particles;
    for(auto block: data.grid.grid_blocks){
        for(Particle particle: block.block_particles_v1){
            grid_particles.push_back(particle);
        }
    }
    sortParticlesById(grid_particles);
    for (const auto& particle : grid_particles) {
        float posx = static_cast<float>(particle.pos[0]);float posy = static_cast<float>(particle.pos[1]);float posz = static_cast<float>(particle.pos[2]);
        float hvx = static_cast<float>(particle.hv[0]);float hvy = static_cast<float>(particle.hv[1]);float hvz = static_cast<float>(particle.hv[2]);
        float velx = static_cast<float>(particle.vel[0]);float vely = static_cast<float>(particle.vel[1]);float velz = static_cast<float>(particle.vel[2]);
        outputFile.write(reinterpret_cast<const char*>(&posx), sizeof(posx));//NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
        outputFile.write(reinterpret_cast<const char*>(&posy), sizeof(posy));//NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
        outputFile.write(reinterpret_cast<const char*>(&posz), sizeof(posz));//NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
        outputFile.write(reinterpret_cast<const char*>(&hvx), sizeof(hvx));//NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
        outputFile.write(reinterpret_cast<const char*>(&hvy), sizeof(hvy));//NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
        outputFile.write(reinterpret_cast<const char*>(&hvz), sizeof(hvz));//NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
        outputFile.write(reinterpret_cast<const char*>(&velx), sizeof(velx));//NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
        outputFile.write(reinterpret_cast<const char*>(&vely), sizeof(vely));//NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
        outputFile.write(reinterpret_cast<const char*>(&velz), sizeof(velz));//NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
    }
    outputFile.close();
}
