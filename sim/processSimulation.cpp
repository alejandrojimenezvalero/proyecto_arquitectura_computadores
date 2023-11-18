//

#include "particle.hpp"
#include "grid.hpp"
#include "block.hpp"
#include "initSimulation.hpp"
#include "constants.hpp"

#include <unordered_map>
#include <map>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <chrono>

using namespace std;
using namespace simulationConstants;

bool blockExists(int i, int j, int k, Grid& grid){
  if ((i < 0) or (i > grid.grid_dimensions[0] - 1) or (j < 0) or (i > grid.grid_dimensions[1] - 1) or (k < 0) or (k > grid.grid_dimensions[2] - 1)){return false;}
  return true;
}
/*
void createAdjacentBlocks(Grid& grid) {
    map<vector<int>, vector<vector<int>>> adjacentBlocks;
    for (int i =0; i < grid.grid_dimensions[0]; ++i) {
        for (int j=0; j < grid.grid_dimensions[1]; ++j) {
            for (int k=0; k < grid.grid_dimensions[2]; ++k){
                vector<vector<int>> addVector;
                for (int di = -1; di <= 1; di++) {
                    for (int dj = -1; dj <= 1; dj++) {
                        for (int dk = -1; dk <= 1; dk++) {
                            int target_i = i + di;
                            int target_j = j + dj;
                            int target_k = k + dk;
                            if (blockExists(target_i, target_j, target_k, grid)) {
                                // Agregar las coordenadas del bloque adyacente al vector
                                addVector.push_back({target_i, target_j, target_k});
                            }
                        }
                    }
                }
                std::vector<int> key = {i, j, k};
                adjacentBlocks[key]  = addVector;
            }
        }
    }
    return adjacentBlocks;
}*/

double calculateNorm(const std::vector<double>& particlei, const std::vector<double>& particlej) {
    double sumOfSquares = 0.0;
    for (int i = 0; i < 3; i++) {
        double diff = particlei[i] - particlej[i];
        sumOfSquares += diff * diff;
    }
    return sumOfSquares;
}
void transformDensity(Particle& particle, const SimulationData& data){
  particle.density = (particle.density + data.smoothing_length_6) * (315/(64*PI*data.smoothing_length_9)) * data.particle_mass;
}

void transferAcceleration(Particle& particlei, Particle& particlej, double& dist, const SimulationData& data) {
    const double escalar_pos = (15 / (simulationConstants::PI * data.smoothing_length_6)) * ((3 * data.particle_mass * simulationConstants::STIFFNESS_PRESSURE) / 2) *
                               ((std::pow(data.smoothing_length - dist, 2)) / dist) * (particlei.density + particlej.density - 2 * simulationConstants::FLUID_DENSITY);
    const double escalar_vel = (45 / (simulationConstants::PI * data.smoothing_length_6)) * simulationConstants::VISCOSITY * data.particle_mass;

    std::vector<double> diff_pos = {particlei.pos[0] - particlej.pos[0], particlei.pos[1] - particlej.pos[1], particlei.pos[2] - particlej.pos[2]};
    std::vector<double> diff_vel = {particlej.vel[0] - particlei.vel[0], particlej.vel[1] - particlei.vel[1], particlej.vel[2] - particlei.vel[2]};

    std::vector<double> variation_acc = {
            (diff_pos[0] * escalar_pos + diff_vel[0] * escalar_vel) / (particlei.density * particlej.density),
            (diff_pos[1] * escalar_pos + diff_vel[1] * escalar_vel) / (particlei.density * particlej.density),
            (diff_pos[2] * escalar_pos + diff_vel[2] * escalar_vel) / (particlei.density * particlej.density)
    };

    particlei.acceleration = {particlei.acceleration[0] + variation_acc[0], particlei.acceleration[1] + variation_acc[1], particlei.acceleration[2] + variation_acc[2]};
    particlej.acceleration = {particlej.acceleration[0] - variation_acc[0], particlej.acceleration[1] - variation_acc[1], particlej.acceleration[2] - variation_acc[2]};
}


/*void initializeDensityAcceleration(SimulationData& data){
    std::vector<int> block_index;
    std::shared_ptr<std::vector<Particle>> block_particles;
    for (auto& block: data.grid.grid_blocks){
        for(Particle& particle: block.second.block_particles){
            //if (particle.id == 204){cout << 204 << " initializeDensityAcceleration" << '\n';}
            particle.density = 0.0;
            particle.density_updated = false;
            particle.acceleration[0] = 0.0; particle.acceleration[1] = -9.8; particle.acceleration[2] = 0.0;
            particle.acceleration_updated = false;
        }
    }
}*/
/*void updateParticle(SimulationData& data, Particle& particle, int index_adj, string mode){
    double norm_2, dist;
    for(Particle& adj_particle: data.grid.grid_blocks[index_adj].block_particles){

        norm_2 = calculateNorm(particle.pos, adj_particle.pos);
        if (mode =="density"){
            if ((adj_particle.density_updated == false) && particle.id != adj_particle.id){
                if(norm_2<data.smoothing_length_2){
                    //if (particle.id == 204 ){cout << "Particle: " << particle.density << '\n';}
                    //if (adj_particle.id == 204){cout << "Adj_Particle: " << adj_particle.density << '\n';}
                    particle.density += pow((data.smoothing_length_2-norm_2),3), adj_particle.density += pow((data.smoothing_length_2-norm_2),3);
                }
                else{particle.density += 0;}
            }
        }
        else if (mode=="acceleration"){
            if (particle.density_updated && adj_particle.density_updated && adj_particle.acceleration_updated == false && particle != adj_particle) {
                if (norm_2 < data.smoothing_length_2) {
                    dist = sqrt(max(norm_2, MINIMUM_DISTANCE));
                    transferAcceleration(particle, adj_particle, dist, data);
                }
            }
        }
    }
}*/

void updateBlocksDensity(SimulationData& data){
    double norm_2;
    //int c=0; int k= 0;
    for(auto& block: data.grid.grid_blocks){
        for (Particle& particle: block.second.block_particles){
            for(vector<int> index_adj:block.second.adj_blocks){
                Block& adj_block = data.grid.grid_blocks[index_adj];
                //cout << "index_adj: " << index_adj[0] << ", " << index_adj[1] << ", " << index_adj[2] << ", " <<'\n';
                //cout << "index_block: " << block.second.block_index[0] << ", " << block.second.block_index[1] << ", " << block.second.block_index[2] << ", " <<'\n';
                if(!adj_block.updated_density){
                    //if (particle.id == 204 ){cout << "Particle 204: " << particle.density << '\n';}
                    for(Particle& adj_particle: adj_block.block_particles){
                        norm_2 = calculateNorm(particle.pos, adj_particle.pos);
                        //k+=1;
                        if (!adj_particle.density_updated && particle.id != adj_particle.id){
                            //c+=1;
                            if(norm_2<data.smoothing_length_2){
                                double added_norm = std::pow((data.smoothing_length_2-norm_2),3);
                                //if (particle.id == 204 ){cout << "Particle: " << particle.density << '\n';}
                                //if (adj_particle.id == 204){cout << "Adj_Particle: " << adj_particle.density << '\n';}
                                particle.density += added_norm, adj_particle.density += added_norm;
                            }
                            else{particle.density += 0;}
                        }
                    }
                }

            }
            transformDensity(particle, data);
            particle.density_updated = true;
            particle.acceleration_updated = false;
            //cout << "Id: "<< particle.id << ", Density: " << particle.density << '\n';

            /*if (particle.id == 204){
                cout << "Id: "<< particle.id << ", Density: " << particle.density << '\n';
                cout << "Pos: " << particle.pos[0] << ", " << particle.pos[1] << ", " << particle.pos[2] << '\n';
                cout << "ACC: " << particle.acceleration[0] << ", " << particle.acceleration[1] << ", " << particle.acceleration[2] << '\n';
            }*/
        }
        block.second.updated_density = true;
    }
    data.all_particles_density_updated = true;
    //cout << "c: " << c << '\n';
    //cout << "k: " << k << '\n';
}
void updateBlocksAcceleration(SimulationData& data){
    double norm_2, dist;
    for(auto& block: data.grid.grid_blocks){
        std::vector<Particle>& particles = block.second.block_particles;
        for (size_t i = 0; i < particles.size(); ++i){
            Particle& particle = particles[i];
            for(vector<int> index_adj:block.second.adj_blocks){
                Block& adj_block = data.grid.grid_blocks[index_adj];
                //cout << "index_adj: " << index_adj[0] << ", " << index_adj[1] << ", " << index_adj[2] << ", " <<'\n';
                //cout << "index_block: " << block.second.block_index[0] << ", " << block.second.block_index[1] << ", " << block.second.block_index[2] << ", " <<'\n';
                //k+=1;
                if(!adj_block.updated_acceleration){
                    //c+=1;
                    //if (particle.id == 204 ){cout << "Particle 204: " << particle.density << '\n';}
                    for(Particle& adj_particle: adj_block.block_particles){

                        norm_2 = calculateNorm(particle.pos, adj_particle.pos);

                        if (!adj_particle.acceleration_updated && particle != adj_particle) {
                            if (norm_2 < data.smoothing_length_2) {
                                dist = sqrt(max(norm_2, MINIMUM_DISTANCE));
                                transferAcceleration(particle, adj_particle, dist, data);
                            }
                        }

                    }
                }

            }
            particle.acceleration_updated = true;
            /*if (particle.id == 204){
                cout << "Id: "<< particle.id << ", Density: " << particle.density << '\n';
                cout << "Pos: " << particle.pos[0] << ", " << particle.pos[1] << ", " << particle.pos[2] << '\n';
                cout << "ACC: " << particle.acceleration[0] << ", " << particle.acceleration[1] << ", " << particle.acceleration[2] << '\n';
            }*/
        }
        block.second.updated_acceleration=true;
    }
    //cout << "c: " << c << '\n';
    //cout << "k: " << k << '\n';
}
double calcCord(Particle& particle, int index){
    double cord;
    cord = particle.pos[index] + particle.hv[index]*TIMESTAMP;
    return cord;
}
double calcVariation(int& index_block, double& cordParticle, const SimulationData& data, int& index){
    double var = 0.0;
    if(index_block == 0){
        var = PARTICLE_SIZE - (cordParticle - LOWER_LIMIT[index]);
    }
    else if(index_block == data.grid.grid_dimensions[index] - 1){
        var = PARTICLE_SIZE - (UPPER_LIMIT[index] - cordParticle);
    }
    return var;
}
void calcAcceleration(Particle& particle, double& var, SimulationData& data, int& index){
    double v;
    int cordBlock = calcParticleIndex(particle, data.grid)[index];
    v = particle.vel[index];
    if(cordBlock == 0){
        particle.acceleration[index] += (STIFFNESS_COLLISIONS*var - DAMPING*v);
    }
    else if(cordBlock == data.grid.grid_dimensions[index] - 1){
        particle.acceleration[index] -= (STIFFNESS_COLLISIONS*var + DAMPING*v);
    }
}

void updatePosition(Particle& particle, int& index){
    particle.pos[index] += particle.hv[index]*TIMESTAMP + particle.acceleration[index]*TIMESTAMP_2;
}
void updateVelocity(Particle& particle, int& index){
    particle.vel[index] = particle.hv[index] + (particle.acceleration[index]*TIMESTAMP)/2;
}
void updateHv(Particle& particle, int& index){
    particle.hv[index] += particle.acceleration[index]*TIMESTAMP;
}
void checkBorderLimits(Particle& particle, SimulationData& data, int& index, int& actual_block_index){
    double d = 0.0;
    if (actual_block_index==0){
        d = particle.pos[index] - LOWER_LIMIT[index];
    }
    else if (actual_block_index == data.grid.grid_dimensions[index] - 1){
        d = UPPER_LIMIT[index] - particle.pos[index];
    }
    /*if (particle.id == 2229){
        cout << d << '\n';
        cout << "Pos xxx: "  << particle.pos[0] << ", " << particle.pos[1] << ", " << particle.pos[2] << '\n';
        //cout << "Block Index: x:" << actual_block_index[0] << "y: " << actual_block_index[1] << "z: " << actual_block_index[2] << '\n';
    }*/
    if (d<0){
        double particleAxisPos;
        if (actual_block_index == 0){
            particleAxisPos = LOWER_LIMIT[index] - d;
        }
        else if (actual_block_index == data.grid.grid_dimensions[index] - 1){
            particleAxisPos = UPPER_LIMIT[index] + d;
        }
        else{particleAxisPos = particle.pos[index];}
        //if (particle.id == 204){cout << particleAxixPos << '\n';}
        particle.pos[index] = particleAxisPos;
        particle.vel[index] = -particle.vel[index];
        particle.hv[index] = -particle.hv[index];
    }
    /*if (particle.id == 204){
        cout << d << '\n';
        cout << "Pos yyy: "  << particle.pos[0] << ", " << particle.pos[1] << ", " << particle.pos[2] << '\n';}*/
}

void removeParticlesFromBlock(Block& block, std::vector<Particle>& particles_to_remove){
    auto new_end = std::remove_if(block.block_particles.begin(), block.block_particles.end(),
                                  [&](const Particle& particle) {
                                      // La lambda devuelve true si la partícula está en particles_to_remove
                                      return std::find(particles_to_remove.begin(), particles_to_remove.end(), particle) != particles_to_remove.end();
                                  });

    // Utilizar std::erase para eliminar las partículas trasladadas
    block.block_particles.erase(new_end, block.block_particles.end());
}

void updateParticleBlockBelonging(SimulationData& data){
    vector<int> new_block_particle_index;
    std::vector<int>old_block_index;
    std::vector<int>check_block_index{1,0,7};
    for(auto& current_block:data.grid.grid_blocks){
        old_block_index = current_block.second.block_index;
        std::vector<Particle> particles_to_remove{};
        for (Particle& particle: current_block.second.block_particles){
            particle.density = 0.0;
            particle.density_updated = false;
            particle.acceleration[0] = 0.0; particle.acceleration[1] = -9.8; particle.acceleration[2] = 0.0;
            particle.acceleration_updated = false;
            new_block_particle_index = calcParticleIndex(particle, data.grid);
            if (new_block_particle_index != old_block_index){
                //if (particle.id == 204){cout << "Id Block 204: " << new_block_particle_index[0] << ", " << new_block_particle_index[1] << ", " << new_block_particle_index[2] << '\n';}
                //cout << "here " << particle.id <<'\n';
                // añadir a new_block_particles
                data.grid.grid_blocks[new_block_particle_index].block_particles.push_back(particle);
                // eliminar de old_block_particles
                particles_to_remove.push_back(particle);
            }

        }
        current_block.second.updated_density=false;
        current_block.second.updated_acceleration=false;
        if(!particles_to_remove.empty()){removeParticlesFromBlock(current_block.second, particles_to_remove);}
    }
}
void updateParticle(std::vector<int>& block_index, Particle& particle, SimulationData& data){
    double cord, var;
    std::vector<double> ngrids{data.grid.grid_dimensions[0], data.grid.grid_dimensions[1], data.grid.grid_dimensions[2]};
    // Aqui esta la salida del fichero partcol-base
    /*if (particle.id == 2229) {
        cout << "---------------------------------" << '\n';
        cout << "Density before collision: " << particle.density << '\n';
        cout << "Pos before collision: "  << particle.pos[0] << ", " << particle.pos[1] << ", " << particle.pos[2] << '\n';
        cout << "Hv before collision: "  << particle.hv[0] << ", " << particle.hv[1] << ", " << particle.hv[2] << '\n';
        cout << "V before collision: "  << particle.vel[0] << ", " << particle.vel[1] << ", " << particle.vel[2] << '\n';
        cout << "Acc before collision: "  << particle.acceleration[0] << ", " << particle.acceleration[1] << ", " << particle.acceleration[2] << '\n';
        cout << "---------------------------------" << '\n';
    }*/
    /*for(int i=0; i < 3 ; ++i) {
        cord = calcCord(particle, i);
        var = calcVariation(block_index[i], cord, data, i);
        if (var > MINIMUM_VARIATION) { calcAcceleration(particle, var, data, i);}
        updatePosition(particle, i);
        updateVelocity(particle, i);
        updateHv(particle, i);
        checkBorderLimits(particle, data, i, block_index[i]);
    }*/
    for(int i=0; i < 3 ; ++i) {
        cord = calcCord(particle, i);
        var = calcVariation(block_index[i], cord, data, i);
        if (var > MINIMUM_VARIATION) { calcAcceleration(particle, var, data, i);}
    }
    /*if (particle.id == 2229) {
        cout << "---------------------------------" << '\n';
        cout << "Density after collision/before movement: " << particle.density << '\n';
        cout << "Pos after collision/before movement: "  << particle.pos[0] << ", " << particle.pos[1] << ", " << particle.pos[2] << '\n';
        cout << "Hv after collision/before movement: "  << particle.hv[0] << ", " << particle.hv[1] << ", " << particle.hv[2] << '\n';
        cout << "V after collision/before movement: "  << particle.vel[0] << ", " << particle.vel[1] << ", " << particle.vel[2] << '\n';
        cout << "Acc after collision/before movement: "  << particle.acceleration[0] << ", " << particle.acceleration[1] << ", " << particle.acceleration[2] << '\n';
        cout << "---------------------------------" << '\n';
    }*/

    for(int i=0; i < 3 ; ++i) {
        updatePosition(particle, i);
        updateVelocity(particle, i);
        updateHv(particle, i);
    }

    // Aqui está la salida del fichero motion-base

     /*if (particle.id == 2229) {
        cout << "Density after movement/before limit process: " << particle.density << '\n';
        cout << "Pos after movement/before limit process: "  << particle.pos[0] << ", " << particle.pos[1] << ", " << particle.pos[2] << '\n';
        cout << "Hv after movement/before limit process: "  << particle.hv[0] << ", " << particle.hv[1] << ", " << particle.hv[2] << '\n';
        cout << "V after movement/before limit process: "  << particle.vel[0] << ", " << particle.vel[1] << ", " << particle.vel[2] << '\n';
        cout << "Acc after movement/before limit process: "  << particle.acceleration[0] << ", " << particle.acceleration[1] << ", " << particle.acceleration[2] << '\n';
        cout << "---------------------------------" << '\n';
    }*/

    for(int i=0; i < 3 ; ++i){
        checkBorderLimits(particle, data, i, block_index[i]);
    }
    // Aqui está la salida del fichero boundint-base
    /*if (particle.id == 2229) {
        cout << "Density after limit process: " << particle.density << '\n';
        cout << "Pos after limit process: "  << particle.pos[0] << ", " << particle.pos[1] << ", " << particle.pos[2] << '\n';
        cout << "Hv after limit process: "  << particle.hv[0] << ", " << particle.hv[1] << ", " << particle.hv[2] << '\n';
        cout << "V after limit process: "  << particle.vel[0] << ", " << particle.vel[1] << ", " << particle.vel[2] << '\n';
        cout << "Acc after limit process: "  << particle.acceleration[0] << ", " << particle.acceleration[1] << ", " << particle.acceleration[2] << '\n';
        cout << "---------------------------------" << '\n';
    }*/
}

void establishParticleFunctionality(SimulationData& data){
    std::vector<int>block_index;
    std::shared_ptr<std::vector<Particle>> block_particles;
    for(auto& block:data.grid.grid_blocks){
        for (Particle& particle: block.second.block_particles){
            //if (particle.id == 204){cout << 204 << '\n';}
            updateParticle(block.second.block_index, particle, data);
        }
    }
    //cout << "here-4" << '\n';
}


int processSimulation(SimulationData& data){
    //map<vector<int>,vector<vector<int>>> adjacent_blocks = createAdjacentBlocks(data.grid);

    //cout << "check" << '\n';
    //auto start = std::chrono::high_resolution_clock::now();
    updateParticleBlockBelonging(data);
    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    //std::cout << "updateParticleBelonging: " << duration.count() << '\n';

    //auto start2 = std::chrono::high_resolution_clock::now();
    updateBlocksDensity(data);
    //auto end2 = std::chrono::high_resolution_clock::now();
    //auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2);
    //std::cout << "updateBlocksDensity: " << duration2.count() << '\n';

    //auto start3 = std::chrono::high_resolution_clock::now();
    updateBlocksAcceleration(data);
    //auto end3 = std::chrono::high_resolution_clock::now();
    //auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(end3 - start3);
    //std::cout << "updateBlocksAcceleration: " << duration3.count() << '\n';

    //std::cout << "-----------END ACCELERATION---------------" << '\n';
    data.all_particles_density_updated = false;


    establishParticleFunctionality(data);


    //initializeDensityAcceleration(data);
    //std::cout << "-----------END ITERATION---------------" << '\n';
    return 0;
}
