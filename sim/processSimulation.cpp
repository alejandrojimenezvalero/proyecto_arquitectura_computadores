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


double calculateNorm(const std::vector<double>& particlei, const std::vector<double>& particlej) {
    double sumOfSquares = 0.0;
    for (int i = 0; i < 3; i++) {
        sumOfSquares += (particlei[i] - particlej[i])*(particlei[i] - particlej[i]);
    }
    return sumOfSquares;
}
void transformDensity(Particle& particle, const SimulationData& data){
  particle.density = (particle.density + data.smoothing_length_6) * data.escalar_density;
}

void transferAcceleration(Particle& particlei, Particle& particlej, const double& dist, const SimulationData& data) {
    const double escalar_pos = data.escalar_pos *
                               ((std::pow(data.smoothing_length - dist, 2)) / dist) * (particlei.density + particlej.density - 2 * simulationConstants::FLUID_DENSITY);
    const double escalar_vel = data.escalar_vel * simulationConstants::VISCOSITY * data.particle_mass;
    const double density_product = particlei.density * particlej.density;

    std::vector<double> diff_pos = {particlei.pos[0] - particlej.pos[0], particlei.pos[1] - particlej.pos[1], particlei.pos[2] - particlej.pos[2]};
    std::vector<double> diff_vel = {particlej.vel[0] - particlei.vel[0], particlej.vel[1] - particlei.vel[1], particlej.vel[2] - particlei.vel[2]};

    std::vector<double> variation_acc = {
            (diff_pos[0] * escalar_pos + diff_vel[0] * escalar_vel) / (density_product),
            (diff_pos[1] * escalar_pos + diff_vel[1] * escalar_vel) / (density_product),
            (diff_pos[2] * escalar_pos + diff_vel[2] * escalar_vel) / (density_product)
    };

    particlei.acceleration = {particlei.acceleration[0] + variation_acc[0], particlei.acceleration[1] + variation_acc[1], particlei.acceleration[2] + variation_acc[2]};
    particlej.acceleration = {particlej.acceleration[0] - variation_acc[0], particlej.acceleration[1] - variation_acc[1], particlej.acceleration[2] - variation_acc[2]};
}

void updateBlocksDensity(SimulationData& data){
    //int c=0; int k= 0;
    for(Block& block: data.grid.grid_blocks){
        for (Particle& particle: block.block_particles){
            for(const int index_adj:block.adj_index_vector_blocks){
                Block& adj_block = data.grid.grid_blocks[index_adj];
                //cout << "index_adj: " << index_adj[0] << ", " << index_adj[1] << ", " << index_adj[2] << ", " <<'\n';
                //cout << "index_block: " << block.second.block_index[0] << ", " << block.second.block_index[1] << ", " << block.second.block_index[2] << ", " <<'\n';
                if(!adj_block.updated_density){
                    //if (particle.id == 204 ){cout << "Particle 204: " << particle.density << '\n';}
                    for(Particle& adj_particle: adj_block.block_particles){

                        //k+=1;
                        if (!adj_particle.density_updated && particle.id != adj_particle.id){
                            double norm_2 = calculateNorm(particle.pos, adj_particle.pos);
                            //c+=1;
                            if(norm_2<data.smoothing_length_2){
                                const double added_norm = std::pow((data.smoothing_length_2-norm_2),3);
                                //if (particle.id == 204 ){cout << "Particle: " << particle.density << '\n';}
                                //if (adj_particle.id == 204){cout << "Adj_Particle: " << adj_particle.density << '\n';}
                                particle.density += added_norm, adj_particle.density += added_norm;
                            }
                        }
                    }
                }

            }
            transformDensity(particle, data);
            particle.density_updated = true;
            particle.already_checked = false;
            //cout << "Id: "<< particle.id << ", Density: " << particle.density << '\n';

            /*if (particle.id == 204){
                cout << "Id: "<< particle.id << ", Density: " << particle.density << '\n';
                cout << "Pos: " << particle.pos[0] << ", " << particle.pos[1] << ", " << particle.pos[2] << '\n';
                cout << "ACC: " << particle.acceleration[0] << ", " << particle.acceleration[1] << ", " << particle.acceleration[2] << '\n';
            }*/
        }
        block.updated_density = true;
    }
    //cout << "c: " << c << '\n';
    //cout << "k: " << k << '\n';
}
void updateBlocksAcceleration(SimulationData& data){
    for(Block& block: data.grid.grid_blocks){
        for (Particle& particle: block.block_particles){
            for(int index_adj:block.adj_index_vector_blocks){
                Block& adj_block = data.grid.grid_blocks[index_adj];
                //cout << "index_adj: " << index_adj[0] << ", " << index_adj[1] << ", " << index_adj[2] << ", " <<'\n';
                //cout << "index_block: " << block.second.block_index[0] << ", " << block.second.block_index[1] << ", " << block.second.block_index[2] << ", " <<'\n';
                //k+=1;
                if(!adj_block.updated_acceleration){
                    //c+=1;
                    //if (particle.id == 204 ){cout << "Particle 204: " << particle.density << '\n';}
                    for(Particle& adj_particle: adj_block.block_particles){

                        if (!adj_particle.acceleration_updated && particle.id != adj_particle.id) {
                            const double norm_2 = calculateNorm(particle.pos, adj_particle.pos);
                            if (norm_2 < data.smoothing_length_2) {
                                const double dist = sqrt(max(norm_2, MINIMUM_DISTANCE));
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
        block.updated_acceleration=true;
    }
    //cout << "c: " << c << '\n';
    //cout << "k: " << k << '\n';
}
double calcCord(Particle& particle, int index){
    const double cord = particle.pos[index] + particle.hv[index]*TIMESTAMP;
    return cord;
}
double calcVariation(int& index_block, const double& cordParticle, const SimulationData& data, int& index){
    double var = 0.0;
    if(index_block == 0){
        var = PARTICLE_SIZE - (cordParticle - LOWER_LIMIT[index]);
    }
    else if(index_block == data.grid.grid_dimensions[index] - 1){
        var = PARTICLE_SIZE - (UPPER_LIMIT[index] - cordParticle);
    }
    return var;
}
void calcAcceleration(Particle& particle, const double& var, SimulationData& data, int& index){
    const double v = particle.vel[index];
    const int cordBlock = calcParticleIndex(particle, data.grid)[index];
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
void checkBorderLimits(Particle& particle, SimulationData& data, int index, int actual_block_index) {
    double d = 0.0;
    const double lower_limit = LOWER_LIMIT[index];
    const double upper_limit = UPPER_LIMIT[index];

    if (actual_block_index == 0) {
        d = particle.pos[index] - lower_limit;
    } else if (actual_block_index == data.grid.grid_dimensions[index] - 1) {
        d = upper_limit - particle.pos[index];
    }
    /*if (particle.id == 2229){
        cout << d << '\n';
        cout << "Pos xxx: "  << particle.pos[0] << ", " << particle.pos[1] << ", " << particle.pos[2] << '\n';
        //cout << "Block Index: x:" << actual_block_index[0] << "y: " << actual_block_index[1] << "z: " << actual_block_index[2] << '\n';
    }*/
    if (d < 0) {
        double particleAxisPos;
        if (actual_block_index == 0) {
            particleAxisPos = lower_limit - d;
        } else if (actual_block_index == data.grid.grid_dimensions[index] - 1) {
            particleAxisPos = upper_limit + d;
        } else {
            particleAxisPos = particle.pos[index];
        }

        particle.pos[index] = particleAxisPos;
        particle.vel[index] = -particle.vel[index];
        particle.hv[index] = -particle.hv[index];
    }
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
/*void removeParticlesFromBlock(Block& block, const std::vector<Particle>& particles_to_remove) {
    // Utilizar la función erase-remove idiom para eliminar las partículas por igualdad de id
    block.block_particles.erase(std::remove_if(block.block_particles.begin(), block.block_particles.end(),
                                               [&particles_to_remove](const Particle& particle) {
                                                   return std::find(particles_to_remove.begin(), particles_to_remove.end(), particle) != particles_to_remove.end();
                                               }),
                                block.block_particles.end());
}*/

void updateParticleBlockBelonging(SimulationData& data){
    //auto start = std::chrono::high_resolution_clock::now();
    //auto second_for = std::chrono::microseconds(0);
    //auto total_remove_time = std::chrono::microseconds(0);
    //auto total_push_back_new = std::chrono::microseconds(0);
    //auto total_calc_particle_index = std::chrono::microseconds(0);
    //auto total_reset_parameters = std::chrono::microseconds(0);
    //auto total_index_in_vector = std::chrono::microseconds(0);
    //auto total_push_back_to_remove = std::chrono::microseconds(0);
    //auto total_if = std::chrono::microseconds(0);
    for (std::size_t i = 0; i < data.grid.grid_blocks.size(); i += 1){
        // Acceder al bloque actual
        Block& current_block = data.grid.grid_blocks[i];
        const vector<int>& current_block_index = current_block.block_index;
        std::vector<Particle> particles_to_remove{};
        //auto start9 = std::chrono::high_resolution_clock::now();
        for (Particle& particle: current_block.block_particles){
            if(!particle.already_checked){
                //auto start7 = std::chrono::high_resolution_clock::now();
                const vector<int>& new_block_particle_index = calcParticleIndex(particle, data.grid);
                //auto end7 = std::chrono::high_resolution_clock::now();
                //auto duration7 = std::chrono::duration_cast<std::chrono::microseconds>(end7 - start7);
                //total_calc_particle_index += duration7;
                //auto start8 = std::chrono::high_resolution_clock::now();
                particle.density = 0.0;
                particle.density_updated = false;
                particle.acceleration[0] = 0.0; particle.acceleration[1] = -9.8; particle.acceleration[2] = 0.0;
                particle.acceleration_updated = false;
                //auto end8 = std::chrono::high_resolution_clock::now();
                //auto duration8 = std::chrono::duration_cast<std::chrono::microseconds>(end8 - start8);
                //total_reset_parameters += duration8;
                //auto start12 = std::chrono::high_resolution_clock::now();
                if (new_block_particle_index != current_block_index){
                    //auto start10 = std::chrono::high_resolution_clock::now();
                    const int index_in_vector = calcParticleIndexVector(data.grid, new_block_particle_index);
                    //auto end10 = std::chrono::high_resolution_clock::now();
                    //auto duration10 = std::chrono::duration_cast<std::chrono::microseconds>(end10 - start10);
                    //total_index_in_vector += duration10;
                    //if (particle.id == 204){cout << "Id Block 204: " << new_block_particle_index[0] << ", " << new_block_particle_index[1] << ", " << new_block_particle_index[2] << '\n';}
                    //cout << "here " << particle.id <<'\n';
                    // añadir a new_block_particles
                    //auto start6 = std::chrono::high_resolution_clock::now();
                    particle.already_checked = true;
                    data.grid.grid_blocks[index_in_vector].block_particles.push_back(particle);
                    //auto end6 = std::chrono::high_resolution_clock::now();
                    //auto duration6 = std::chrono::duration_cast<std::chrono::microseconds>(end6 - start6);
                    //total_push_back_new += duration6;
                    // eliminar de old_block_particles
                    //auto start11 = std::chrono::high_resolution_clock::now();
                    particles_to_remove.push_back(particle);
                    //auto end11 = std::chrono::high_resolution_clock::now();
                    //auto duration11 = std::chrono::duration_cast<std::chrono::microseconds>(end11 - start11);
                    //total_push_back_to_remove += duration11;
                }
                //auto end12 = std::chrono::high_resolution_clock::now();
                //auto duration12 = std::chrono::duration_cast<std::chrono::microseconds>(end12 - start12);
                //total_if += duration12;
            }

        }
        //auto end9 = std::chrono::high_resolution_clock::now();
        //auto duration9 = std::chrono::duration_cast<std::chrono::microseconds>(end9 - start9);
        //second_for += duration9;
        current_block.updated_density=false;
        current_block.updated_acceleration=false;
        //auto start5 = std::chrono::high_resolution_clock::now();
        if (!particles_to_remove.empty()) {
            removeParticlesFromBlock(current_block, particles_to_remove);
        }
        //auto end5 = std::chrono::high_resolution_clock::now();
        //auto duration5 = std::chrono::duration_cast<std::chrono::microseconds>(end5 - start5);

        //total_remove_time += duration5;

        //std::cout << "remove: " << duration5.count() << '\n';
    }
    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    //std::cout << "Total time: " << duration.count() << "\n";
    //std::cout << "Total second for: " << second_for.count() << "\n";
    //std::cout << "Total remove time: " << total_remove_time.count() << "\n";
    //std::cout << "Total pushback new: " << total_push_back_new.count() << "\n";
    //std::cout << "Total calc particle index: " << total_calc_particle_index.count() << "\n";
    //std::cout << "Total reset parameters: " << total_reset_parameters.count() << "\n";
    //std::cout << "Total index in vector: " << total_index_in_vector.count() << "\n";
    //std::cout << "Total pushback remove: " << total_push_back_to_remove.count() << "\n";
    //std::cout << "Total if: " << total_if.count() << "\n";
}
void updateParticle(std::vector<int>& block_index, Particle& particle, SimulationData& data){
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
        const double cord = calcCord(particle, i);
        const double var = calcVariation(block_index[i], cord, data, i);
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
    for(Block& block:data.grid.grid_blocks){
        for (Particle& particle: block.block_particles){
            //if (particle.id == 204){cout << 204 << '\n';}
            updateParticle(block.block_index, particle, data);
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

    //auto start4 = std::chrono::high_resolution_clock::now();
    establishParticleFunctionality(data);
    //auto end4 = std::chrono::high_resolution_clock::now();
    //auto duration4 = std::chrono::duration_cast<std::chrono::microseconds>(end4 - start4);
    //std::cout << "establishParticleFunctionality: " << duration4.count() << '\n';

    //initializeDensityAcceleration(data);
    //std::cout << "-----------END ITERATION---------------" << '\n';
    return 0;
}
