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
#include <numeric>
#include <tuple>

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
}
void createSubMap(Grid& grid, vector<vector<int>> adjacent_blocks, SubGrid& particleSubMap){
    //cout << "-------------------------" << '\n';
    //cout << current_block_key[0] << ", " << current_block_key[1] << ", " << current_block_key[2] << '\n';
    //cout << "-------------------------" << '\n';
    std::vector<std::reference_wrapper<Block>> insert_wrapper{};
    for (auto& adj_block_key: adjacent_blocks){
        for (auto& adj_block: grid.grid_blocks){
            if (adj_block_key == adj_block.block_index){
                insert_wrapper.emplace_back(std::ref(adj_block));
            }
        }

    }
    particleSubMap.particleSubMap = insert_wrapper;
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

void transferAcceleration(Particle& particlei, Particle& particlej, double dist, const SimulationData& data){
  std::vector<double> variation_acc{};
  double escalar_pos, escalar_vel;
  vector<double> diff_pos = {particlei.pos[0] - particlej.pos[0], particlei.pos[1] - particlej.pos[1], particlei.pos[2] - particlej.pos[2]};
  vector<double> diff_vel = {particlej.vel[0] - particlei.vel[0], particlej.vel[1] - particlei.vel[1], particlej.vel[2] - particlei.vel[2]};
  escalar_pos = ((15/(PI*data.smoothing_length_6))*((3*data.particle_mass*STIFFNESS_PRESSURE)/2)*((pow(data.smoothing_length-dist,2))/dist)*(particlei.density+particlej.density-2*FLUID_DENSITY));
  escalar_vel = (45/(PI*data.smoothing_length_6))*VISCOSITY*data.particle_mass;
  variation_acc = {(diff_pos[0]*(escalar_pos)+diff_vel[0]*(escalar_vel))/(particlei.density*particlej.density),
                   (diff_pos[1]*(escalar_pos)+diff_vel[1]*(escalar_vel))/(particlei.density*particlej.density),
                   (diff_pos[2]*(escalar_pos)+diff_vel[2]*(escalar_vel))/(particlei.density*particlej.density)};
    particlei.acceleration = {particlei.acceleration[0] + variation_acc[0],particlei.acceleration[1] + variation_acc[1],particlei.acceleration[2] + variation_acc[2]};
    particlej.acceleration = {particlej.acceleration[0] - variation_acc[0],particlej.acceleration[1] - variation_acc[1],particlej.acceleration[2] - variation_acc[2]};
}
/*
void updateBlock(std::vector<Particle> current_block_particles, std::map<std::vector<int>, std::vector<Particle>> particleSubMap, SimulationData data, string mode){
    // IMPORTANTE -> Esta función necesita que al inicio de cada iteracion de processSimulation se utilice la funcion updateDensityFalse
    // REFACTORIZAR CON VECTORES
    double h = data.smoothing_length;
    double norm, dist;
    for(Particle particle: current_block_particles){
        for(auto block: particleSubMap){
            vector<Particle> adj_particles = block.second;
            for(auto adj_particle: adj_particles){
                if (mode =="density"){
                    if (!adj_particle.density_updated && particle != adj_particle){
                        norm = calculateNorm({particle.px,particle.py,particle.pz}, {adj_particle.px, adj_particle.py, adj_particle.pz});
                        if(pow(norm,2)<pow(h,2)){
                            particle.density += pow((pow(h,2)-pow(norm,2)), 3), adj_particle.density += pow((pow(h,2)-pow(norm,2)), 3);
                            particle.density = transformDensity(particle.density, data), adj_particle.density = transformDensity(adj_particle.density, data);
                        }
                    }
                }
                else if (mode=="acceleration"){
                    if (particle.density_updated && adj_particle.density_updated && !adj_particle.acceleration_updated && particle != adj_particle) {
                        if (pow(norm, 2) < pow(h, 2)) {
                            dist = sqrt(max(pow(norm, 2), pow(10, -12)));
                            std::vector<double> transfered_acceleration = transferAcceleration(particle, adj_particle, dist, data);
                            particle.acceleration = {particle.acceleration[0] + transfered_acceleration[0],particle.acceleration[1] + transfered_acceleration[1],particle.acceleration[2] + transfered_acceleration[2]};
                            adj_particle.acceleration = {adj_particle.acceleration[0] + transfered_acceleration[0],adj_particle.acceleration[1] + transfered_acceleration[1],adj_particle.acceleration[2] + transfered_acceleration[2]};
                        }
                    }
                }
            }
            particle.density_updated = (mode=="density")?true:particle.density_updated;
            particle.acceleration_updated = (mode=="acceleration")?true:false;
        }
    }
}
 */
void initializeDensityAcceleration(SimulationData& data){
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
}

void updateBlocks(SimulationData& data, std::string mode){
    double norm_2, dist;
    int c = 0, k = 0;
    for(auto& block: data.grid.grid_blocks){
        std::vector<Particle>& particles = block.second.block_particles;
        for (size_t i = 0; i < particles.size(); ++i){
            Particle& particle = particles[i];
            for(vector<int> index_adj:block.second.adj_blocks){

                //cout << "index_adj: " << index_adj[0] << ", " << index_adj[1] << ", " << index_adj[2] << ", " <<'\n';
                //cout << "index_block: " << block.second.block_index[0] << ", " << block.second.block_index[1] << ", " << block.second.block_index[2] << ", " <<'\n';
                k+=1;
                if (index_adj != block.second.block_index){
                    c+=1;
                    //if (particle.id == 204 ){cout << "Particle 204: " << particle.density << '\n';}
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
                }

                else if (index_adj == block.second.block_index){
                    //if (particle.id == 204 ){cout << "Particle 204: " << particle.density << '\n';}
                    for (size_t j = i + 1; j < particles.size(); ++j){
                        Particle& adj_particle = particles[j];
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
                }
            }
            if (mode == "density"){
                transformDensity(particle, data);
                particle.density_updated = true;
                particle.acceleration_updated = false;
                //cout << "Id: "<< particle.id << ", Density: " << particle.density << '\n';
            }
            if (mode == "acceleration"){particle.acceleration_updated = true;}
            /*if (particle.id == 204){
                cout << "Id: "<< particle.id << ", Density: " << particle.density << '\n';
                cout << "Pos: " << particle.pos[0] << ", " << particle.pos[1] << ", " << particle.pos[2] << '\n';
                cout << "ACC: " << particle.acceleration[0] << ", " << particle.acceleration[1] << ", " << particle.acceleration[2] << '\n';
            }*/
        }
    }
    cout << "c: " << c << '\n';
    cout << "k: " << k << '\n';
}

/*
// ACTUALIZADA
void updateBlock2(SubGrid& particleSubMap , SimulationData& data, std::string mode){
    // REFACTORIZAR CON VECTORES
    //cout << mode << '\n';
    double norm_2, dist;
    //check_iter=0-> USAMOS density_updated, check_iter=1-> USAMOS density_updated2
    //el using_update_adj, usará en la ADJUNTA density_updated si estamos en check_iter=0, y usará density_updated2 si estamos en check_iter=1
    for(Particle& particle : particleSubMap.current_block.block_particles){
        //if (particle.id == 204){cout << 204 << " updateBlock2" << '\n';}
        for (std::reference_wrapper<Block> refBlock  : particleSubMap.particleSubMap) {
            Block& block = refBlock.get();
            for(Particle& adj_particle: block.block_particles){
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
        }
        if (mode == "density"){
            transformDensity(particle, data);
            particle.density_updated = true;
            particle.acceleration_updated = false;
            //cout << "Id: "<< particle.id << ", Density: " << particle.density << '\n';
            }
        if (mode == "acceleration"){particle.acceleration_updated = true;}
        if (particle.id == 204){
            cout << "Id: "<< particle.id << ", Density: " << particle.density << '\n';
            cout << "Pos: " << particle.pos[0] << ", " << particle.pos[1] << ", " << particle.pos[2] << '\n';
            cout << "ACC: " << particle.acceleration[0] << ", " << particle.acceleration[1] << ", " << particle.acceleration[2] << '\n';
        }
        //if (mode == "acceleration"){cout << "Id: "<< particle.id << ", Acceleration: " << particle.acceleration[0] << "," << particle.acceleration[1] << "," << particle.acceleration[2]  << '\n';}
    }
}*/
/*
// ACTUALIZADA
int modifyBlock(Block& current_block, std::map<std::vector<int>, vector<vector<int>>>& adjacent_blocks_map, SimulationData& data) {

    // En el mapa de bloques adyacentes tenemos pares clave-valor de la forma bloque {i,j,k} = vector bloques adyacentes{{i,j,k}, {i,j,k},...}
    // Tomamos un vector de bloques adyacentes al bloque en el que está la partícula que se recibe como parámetro
    // Así tomamos un submapa del mapa de partículas de manera que nos quedamos con los pares {i,j,k} = {P1,P2,..} para todos los bloques adyacentes
    // Ahora queremos meter el vector de bloques referenciados en el SubMapa
    SubGrid particleSubMap(current_block);
    createSubMap(data.grid, adjacent_blocks_map[current_block.block_index], particleSubMap);
    Teniendo el submapa con claves indices de bloques adyacentes al bloque de la partícula parámetro , y valores las partículas de dichos bloques adyacentes
    debemos iterar por los bloques del submapa, y para cada partícula de cada bloque adyacente calcular la densidad con respecto a la partícula parámetro
    if (!data.all_particles_density_updated){
        //updateBlock(current_block_particles, particleSubMap, data, "density");
        updateBlock2(particleSubMap, data, "density");
    }
    else{
        //updateBlock(current_block_particles, particleSubMap, data, "acceleration");
        updateBlock2(particleSubMap, data, "acceleration");
    }

    return 0;
}*/


double calcCord(Particle& particle, int index){
    double cord;
    cord = particle.pos[index] + particle.hv[index]*TIMESTAMP;
    return cord;
}
double calcVariation(int index_block, double cordParticle, const SimulationData data, int index){
    double var = 0.0;
    if(index_block == 0){
        var = PARTICLE_SIZE - (cordParticle - LOWER_LIMIT[index]);
    }
    else if(index_block == data.grid.grid_dimensions[index] - 1){
        var = PARTICLE_SIZE - (UPPER_LIMIT[index] - cordParticle);
    }
    return var;
}
void calcAcceleration(Particle& particle, double var, SimulationData data, int index){
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

void updatePosition(Particle& particle, int index){
    particle.pos[index] += particle.hv[index]*TIMESTAMP + particle.acceleration[index]*TIMESTAMP_2;
}
void updateVelocity(Particle& particle, int index){
    particle.vel[index] = particle.hv[index] + (particle.acceleration[index]*TIMESTAMP)/2;
}
void updateHv(Particle& particle, int index){
    particle.hv[index] += particle.acceleration[index]*TIMESTAMP;
}
void checkBorderLimits(Particle& particle, SimulationData& data, int index){
    double d = 0.0;
    int actual_block_index = calcParticleIndex(particle, data.grid)[index];
    if (actual_block_index==0){
        d = particle.pos[index] - LOWER_LIMIT[index];
    }
    else if (actual_block_index == data.grid.grid_dimensions[index] - 1){
        d = UPPER_LIMIT[index] - particle.pos[index];
    }
    /*if (particle.id == 204){
        cout << d << '\n';
        cout << "Pos xxx: "  << particle.pos[0] << ", " << particle.pos[1] << ", " << particle.pos[2] << '\n';}*/
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
            new_block_particle_index = calcParticleIndex(particle, data.grid);
            if (new_block_particle_index != old_block_index){
                //if (particle.id == 204){cout << "Id Block 204: " << new_block_particle_index[0] << ", " << new_block_particle_index[1] << ", " << new_block_particle_index[2] << '\n';}
                //cout << "here " << particle.id <<'\n';
                // añadir a new_block_particles
                for(auto& new_block:data.grid.grid_blocks){
                    if (new_block_particle_index == new_block.second.block_index){
                        new_block.second.block_particles.push_back(particle);
                    }
                }
                // eliminar de old_block_particles
                particles_to_remove.push_back(particle);
            }
        }
        if(!particles_to_remove.empty()){removeParticlesFromBlock(current_block.second, particles_to_remove);}
    }
}
void updateParticle(std::vector<int>block_index, Particle& particle, SimulationData& data){
    double cord, var;
    std::vector<double> ngrids{data.grid.grid_dimensions[0], data.grid.grid_dimensions[1], data.grid.grid_dimensions[2]};
    /*if (particle.id == 204) {
        cout << "---------------------------------" << '\n';
        cout << "Pos before update: "  << particle.pos[0] << ", " << particle.pos[1] << ", " << particle.pos[2] << '\n';
        cout << "Hv before update: "  << particle.hv[0] << ", " << particle.hv[1] << ", " << particle.hv[2] << '\n';
        cout << "V before update: "  << particle.vel[0] << ", " << particle.vel[1] << ", " << particle.vel[2] << '\n';
        cout << "Acc before update: "  << particle.acceleration[0] << ", " << particle.acceleration[1] << ", " << particle.acceleration[2] << '\n';
        cout << "---------------------------------" << '\n';
    }*/
    for(int i=0; i < 3 ; ++i){
        cord = calcCord(particle, i);
        var = calcVariation(block_index[i], cord, data, i);
        if (var>MINIMUM_VARIATION) {calcAcceleration(particle, var, data, i);}
        updatePosition(particle, i);
        updateVelocity(particle, i);
        updateHv(particle, i);

        checkBorderLimits(particle, data, i);
    }
    /*if (particle.id == 204) {
        cout << "Pos after update: "  << particle.pos[0] << ", " << particle.pos[1] << ", " << particle.pos[2] << '\n';
        cout << "Hv after update: "  << particle.hv[0] << ", " << particle.hv[1] << ", " << particle.hv[2] << '\n';
        cout << "V after update: "  << particle.vel[0] << ", " << particle.vel[1] << ", " << particle.vel[2] << '\n';
        cout << "Acc after update: "  << particle.acceleration[0] << ", " << particle.acceleration[1] << ", " << particle.acceleration[2] << '\n';
        cout << "---------------------------------" << '\n';
    }*/
}

void establishParticleFunctionality(SimulationData data){
    std::vector<int>block_index;
    std::shared_ptr<std::vector<Particle>> block_particles;
    for(auto& block:data.grid.grid_blocks){
        for (Particle& particle: block.second.block_particles){
            //if (particle.id == 204){cout << 204 << '\n';}
            updateParticle(block.second.block_index, particle, data);
        }
    }
    updateParticleBlockBelonging(data);
    //cout << "here-4" << '\n';
}


int processSimulation(SimulationData& data){
    //map<vector<int>,vector<vector<int>>> adjacent_blocks = createAdjacentBlocks(data.grid);

    //cout << "check" << '\n';
    updateBlocks(data, "density");
    //std::cout << "-----------END DENSITY---------------" << '\n';
    data.all_particles_density_updated = true;
    updateBlocks(data, "acceleration");

    //std::cout << "-----------END ACCELERATION---------------" << '\n';

    data.all_particles_density_updated = false;
    establishParticleFunctionality(data);
    initializeDensityAcceleration(data);
    //std::cout << "-----------END ITERATION---------------" << '\n';
    return 0;
}
