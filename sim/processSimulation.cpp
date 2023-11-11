//

#include "particle.hpp"
#include "grid.hpp"
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

bool blockExists(int i, int j, int k, gridSize& grid){
  if ((i < 0) or (i > grid.nx-1) or (j < 0) or (i > grid.ny-1) or (k < 0) or (k >grid.nz-1)){return false;}
  return true;
}

map<vector<int>,vector<vector<int>>> createAdjacentBlocks(gridSize& grid) {
  map<vector<int>, vector<vector<int>>> adjacentBlocks;
  for (int i =0; i < grid.nx; ++i) {
    for (int j=0; j < grid.ny; ++j) {
      for (int k=0; k < grid.nz; ++k){
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
// ACTUALIZADA
std::tuple<std::shared_ptr<std::vector<Particle>>, std::map<std::vector<int>, std::shared_ptr<std::vector<Particle>>> > createSubMap(std::vector<int> current_block_key, std::map<std::vector<int>, std::shared_ptr<std::vector<Particle>>>& particleMap, vector<vector<int>> adjacent_blocks){
    std::map<std::vector<int>, std::shared_ptr<std::vector<Particle>>> particleSubMap;
    vector<int> block_key;
    std::shared_ptr<std::vector<Particle>> block_adj_particles;
    std::shared_ptr<std::vector<Particle>> current_block_particles;
    //cout << "-------------------------" << '\n';
    //cout << current_block_key[0] << ", " << current_block_key[1] << ", " << current_block_key[2] << '\n';
    //cout << "-------------------------" << '\n';
    for (auto block: particleMap){
        block_key = block.first;
        block_adj_particles = block.second;
        if (block_key == current_block_key){
            current_block_particles = block.second;
        }
        for (auto adj_block: adjacent_blocks){
            if(block_key == adj_block){
                //cout << block_key[0] << ", " << block_key[1] << ", " << block_key[2] << '\n';
                particleSubMap[adj_block] = block_adj_particles;
            }
        }
    }
    return std::make_tuple(current_block_particles, particleSubMap);
}

double calculateNorm(const std::vector<double>& particlei, const std::vector<double>& particlej) {
    double sumOfSquares = 0.0;
    for (int i = 0; i < 3; i++) {
        double diff = particlei[i] - particlej[i];
        sumOfSquares += diff * diff;
    }
    double norm = std::sqrt(sumOfSquares);
    return norm;
}
double transformDensity(double density, const SimulationData& data){
  density = (density + pow(data.smoothing_length,6)) * (315/(64*PI*pow(data.smoothing_length,9))) * data.particle_mass;
  return density;
}

std::vector<double> transferAcceleration(Particle& particlei, Particle& particlej, double dist, const SimulationData& data){
  std::vector<double> variation_acc{};
  double escalar_pos, escalar_vel;
  vector<double> diff_pos = {particlei.px - particlej.px, particlei.py - particlej.py, particlei.pz - particlej.pz};
  vector<double> diff_vel = {particlej.vx - particlei.vx, particlej.vy - particlei.vy, particlej.vz - particlei.vz};
  escalar_pos = ((15/(PI*pow(data.smoothing_length,6)))*((3*data.particle_mass*STIFFNESS_PRESSURE)/2)*((pow(data.smoothing_length-dist,2))/dist)*(particlei.density+particlej.density-2*FLUID_DENSITY));
  escalar_vel = (45/(PI*pow(data.smoothing_length,6)))*VISCOSITY*data.particle_mass;
  variation_acc = {(diff_pos[0]*(escalar_pos)+diff_vel[0]*(escalar_vel))/(particlei.density*particlej.density),
                   (diff_pos[1]*(escalar_pos)+diff_vel[1]*(escalar_vel))/(particlei.density*particlej.density),
                   (diff_pos[2]*(escalar_pos)+diff_vel[2]*(escalar_vel))/(particlei.density*particlej.density)};
  return variation_acc;
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
void initializeDensityAcceleration(std::map<std::vector<int>, std::shared_ptr<std::vector<Particle>>>& particleMap){
    std::vector<int> block_index;
    std::shared_ptr<std::vector<Particle>> block_particles;
    for (auto block: particleMap){
        block_particles = block.second;
        for(Particle& particle: *block_particles){
            particle.density = 0.0;
            particle.density_updated = false;
            particle.acceleration[0] = 0.0; particle.acceleration[1] = -9.8; particle.acceleration[2] = 0.0;
            particle.acceleration_updated = false;
        }
    }
}
// ACTUALIZADA
void updateBlock2(std::shared_ptr<std::vector<Particle>> current_block_particles, std::map<std::vector<int>, std::shared_ptr<std::vector<Particle>>>& particleSubMap, SimulationData& data, std::string mode){
    // REFACTORIZAR CON VECTORES
    //cout << mode << '\n';
    double h = data.smoothing_length;
    double norm, dist;
    //check_iter=0-> USAMOS density_updated, check_iter=1-> USAMOS density_updated2
    //el using_update_adj, usará en la ADJUNTA density_updated si estamos en check_iter=0, y usará density_updated2 si estamos en check_iter=1
    for(Particle& particle : *current_block_particles){
        for(auto block: particleSubMap){
            std::shared_ptr<std::vector<Particle>> adj_particles = block.second;
            for(Particle& adj_particle: *adj_particles){
                norm = calculateNorm({particle.px,particle.py,particle.pz}, {adj_particle.px, adj_particle.py, adj_particle.pz});
                if (mode =="density"){
                    if ((adj_particle.density_updated == false) && particle != adj_particle){
                        if(pow(norm,2)<pow(h,2)){
                            //if (particle.id == 204 ){cout << "Particle: " << particle.density << '\n';}
                            //if (adj_particle.id == 204){cout << "Adj_Particle: " << adj_particle.density << '\n';}
                            particle.density += pow((pow(h,2)-pow(norm,2)),3), adj_particle.density += pow((pow(h,2)-pow(norm,2)),3);
                        }
                        else{particle.density += 0;}
                    }
                }
                else if (mode=="acceleration"){
                    if (particle.density_updated && adj_particle.density_updated && adj_particle.acceleration_updated == false && particle != adj_particle) {
                        if (pow(norm, 2) < pow(h, 2)) {
                            //if (particle.id ==204){cout << "ACC Particle before: " << particle.acceleration[0] << ", " << particle.acceleration[1] << ", " << particle.acceleration[2] << '\n';}
                            //if (adj_particle.id == 204){cout << "ACC Adj_Particle before: " << adj_particle.acceleration[0] << ", " << adj_particle.acceleration[1] << ", " << adj_particle.acceleration[2] << '\n';}
                            if (particle.id ==2002 or adj_particle.id == 2002){
                                //cout << "BEFORE Id particle: " << particle.id << " density: " << particle.density << " acc: " << particle.acceleration[0] << ", " << particle.acceleration[1] << ", " << particle.acceleration[2] << '\n';
                                //cout << "BEFORE Id adj_particle: " << adj_particle.id << " density: " << adj_particle.density << " acc: " << adj_particle.acceleration[0] << ", " << adj_particle.acceleration[1] << ", " << adj_particle.acceleration[2] << '\n';
                            }
                            dist = sqrt(max(pow(norm, 2), pow(10, -12)));
                            std::vector<double> transfered_acceleration = transferAcceleration(particle, adj_particle, dist, data);
                            particle.acceleration = {particle.acceleration[0] + transfered_acceleration[0],particle.acceleration[1] + transfered_acceleration[1],particle.acceleration[2] + transfered_acceleration[2]};
                            adj_particle.acceleration = {adj_particle.acceleration[0] - transfered_acceleration[0],adj_particle.acceleration[1] - transfered_acceleration[1],adj_particle.acceleration[2] - transfered_acceleration[2]};
                            if (particle.id ==2002 or adj_particle.id == 2002){
                                //cout << "AFTER Id particle: " << particle.id << " density: " << particle.density << " acc: " << particle.acceleration[0] << ", " << particle.acceleration[1] << ", " << particle.acceleration[2] << '\n';
                                //cout << "AFTER Id adj_particle: " << adj_particle.id << " density: " << adj_particle.density << " acc: " << adj_particle.acceleration[0] << ", " << adj_particle.acceleration[1] << ", " << adj_particle.acceleration[2] << '\n';
                            }
                            //if (particle.id ==204){cout << "ACC Particle after: " << particle.acceleration[0] << ", " << particle.acceleration[1] << ", " << particle.acceleration[2] << '\n';}
                            //if (adj_particle.id == 204){cout << "ACC Adj_Particle after: " << adj_particle.acceleration[0] << ", " << adj_particle.acceleration[1] << ", " << adj_particle.acceleration[2] << '\n';}
                        }
                    }
                }
            }
        }
        if (mode == "density"){
            particle.density = transformDensity(particle.density, data);
            particle.density_updated = true;
            //cout << "Id: "<< particle.id << ", Density: " << particle.density << '\n';
            }
        particle.acceleration_updated = (mode=="acceleration")?true:false;
        if (particle.id == 2002){
            cout << "Id: "<< particle.id << ", Density: " << particle.density << '\n';
            cout << "Pos: " << particle.px << ", " << particle.py << ", " << particle.pz << '\n';
            cout << "ACC: " << particle.acceleration[0] << ", " << particle.acceleration[1] << ", " << particle.acceleration[2] << '\n';
        }
        //if (mode == "acceleration"){cout << "Id: "<< particle.id << ", Acceleration: " << particle.acceleration[0] << "," << particle.acceleration[1] << "," << particle.acceleration[2]  << '\n';}
    }
}

// ACTUALIZADA
int modifyBlock(std::vector<int> current_block_key, std::map<std::vector<int>, std::shared_ptr<std::vector<Particle>>>& particleMap, std::map<std::vector<int>,std::vector<std::vector<int>>>& adjacent_blocks_map, SimulationData& data) {

  // En el mapa de bloques adyacentes tenemos pares clave-valor de la forma bloque {i,j,k} = vector bloques adyacentes{{i,j,k}, {i,j,k},...}
  // Tomamos un vector de bloques adyacentes al bloque en el que está la partícula que se recibe como parámetro
  std::vector<vector<int>> adjacent_blocks = adjacent_blocks_map[current_block_key];
  // Así tomamos un submapa del mapa de partículas de manera que nos quedamos con los pares {i,j,k} = {P1,P2,..} para todos los bloques adyacentes
  std::map<std::vector<int>, std::shared_ptr<std::vector<Particle>>> particleSubMap;
  std::shared_ptr<std::vector<Particle>> current_block_particles;
  std::tie(current_block_particles, particleSubMap) = createSubMap(current_block_key, particleMap, adjacent_blocks);
  /*Teniendo el submapa con claves indices de bloques adyacentes al bloque de la partícula parámetro , y valores las partículas de dichos bloques adyacentes
   debemos iterar por los bloques del submapa, y para cada partícula de cada bloque adyacente calcular la densidad con respecto a la partícula parámetro*/
  if (!data.all_particles_density_updated){
      //updateBlock(current_block_particles, particleSubMap, data, "density");
      updateBlock2(current_block_particles, particleSubMap, data, "density");
  }
  else{
      //updateBlock(current_block_particles, particleSubMap, data, "acceleration");
      updateBlock2(current_block_particles, particleSubMap, data, "acceleration");
  }
  return 0;
}


double calcCord(Particle& particle, int index){
    double cord = 0.0;
    if (index==0){cord = particle.px + particle.hvx*TIMESTAMP;}
    else if (index==1){cord = particle.py + particle.hvy*TIMESTAMP;}
    else if (index==2){cord = particle.pz + particle.hvz*TIMESTAMP;}
    return cord;
}
double calcVariation(int index_block, double cordParticle, double ngrid, int index){
    double var = 0.0;

    if(index_block == 0){
        var = PARTICLE_SIZE - (cordParticle - LOWER_LIMIT[index]);
    }
    else if(index_block == ngrid - 1){
        var = PARTICLE_SIZE - (UPPER_LIMIT[index] - cordParticle);
    }
    return var;
}
double calcAcceleration(Particle& particle, double var, double ngrid, int index){
    double v=0.0; double a=0.0;
    int cordBlock = 0;
    a = particle.acceleration[index];
    if (index==0){v = particle.vx;}
    else if (index==1){v = particle.vy;}
    else if (index==2){v = particle.vz;}

    // CAMBIAR ESTO
    if (index==0){cordBlock = particle.i;}
    else if (index==1){cordBlock = particle.j;}
    else if (index==2){cordBlock = particle.k;}
    if(cordBlock == 0){
        a += (STIFFNESS_COLLISIONS*var - DAMPING*v);
    }
    else if(cordBlock == ngrid - 1){
        a -= (STIFFNESS_COLLISIONS*var + DAMPING*v);
    }
    return a;
}

double updatePosition(Particle& particle, int index){
    double p=0.0; double hv=0.0; double a=0.0;
    if (index==0){p = particle.px; hv = particle.hvx;}
    else if (index==1){p = particle.py; hv = particle.hvy;}
    else if (index==2){p = particle.pz; hv = particle.hvz;}
    a = particle.acceleration[index];
    p += hv*TIMESTAMP + a*pow(TIMESTAMP, 2);
    return p;
}
double updateVelocity(Particle& particle, int index){
    double v=0.0; double hv=0.0; double a=0.0;
    if (index==0){hv = particle.hvx;}
    else if (index==1){hv = particle.hvy;}
    else if (index==2){hv = particle.hvz;}
    a = particle.acceleration[index];
    v = hv + (a*TIMESTAMP)/2;
    return v;
}
double updateHv(Particle& particle, int index){
    double hv=0.0; double a=0.0;
    if (index==0){hv = particle.hvx;}
    else if (index==1){hv = particle.hvy;}
    else if (index==2){hv = particle.hvz;}
    a = particle.acceleration[index];
    hv += a*TIMESTAMP;
    return hv;
}
double checkBorderLimits(Particle& particle, double ngrid, int index, int actual_block_index){
    double d = 0.0;
    double p = 0.0;
    if (index==0){p = particle.px;}
    else if (index==1){p = particle.py;}
    else if (index==2){p = particle.pz;}
    if (actual_block_index==0){
        d = p - LOWER_LIMIT[index];
    }
    else if (actual_block_index == ngrid - 1){
        d = UPPER_LIMIT[index] - p;
    }
    if (d<0){
        double particleAxisPos = (actual_block_index==0)?LOWER_LIMIT[index] - d: (actual_block_index == ngrid - 1)?UPPER_LIMIT[index] + d : p;
        if (index==0){particle.px = particleAxisPos, particle.vx = -particle.vx, particle.hvx=-particle.hvx;}
        else if (index==1){particle.py = particleAxisPos, particle.vy = -particle.vy, particle.hvy=-particle.hvy;}
        else if (index==2){particle.pz = particleAxisPos, particle.vz = -particle.vz, particle.hvz=-particle.hvz;}
    }
    return d;
}
void removeParticlesFromBlock(std::shared_ptr<std::vector<Particle>>& old_block_particles, std::vector<Particle>& particles_to_remove){
    old_block_particles->erase(
            std::remove_if(
                    old_block_particles->begin(),
                    old_block_particles->end(),
                    [&particles_to_remove](const Particle& particle) {
                        // Verificar si la partícula está en la lista de partículas a ser eliminadas
                        return std::find(particles_to_remove.begin(), particles_to_remove.end(), particle) != particles_to_remove.end();
                    }
            ),
            old_block_particles->end());
}
void updateParticleBlockBelonging(std::map<std::vector<int>, std::shared_ptr<std::vector<Particle>>>& particleMap, SimulationData& data){
    vector<int> new_block_particle_index;
    std::vector<int>old_block_index;
    std::vector<int>check_block_index{1,0,7};
    std::shared_ptr<std::vector<Particle>> old_block_particles;
    for(auto& current_block:particleMap){
        old_block_index = current_block.first;
        old_block_particles = current_block.second;
        std::vector<Particle> particles_to_remove;
        for (Particle& particle: *old_block_particles){
            new_block_particle_index = calcParticleIndex(particle, data.grid.block_dimensions);
            if (new_block_particle_index != old_block_index){
                cout << "here " << particle.id <<'\n';
                std::shared_ptr<std::vector<Particle>> new_block_particles;
                new_block_particles = particleMap[new_block_particle_index];
                // añadir a new_block_particles
                new_block_particles->emplace_back(particle);
                // eliminar de old_block_particles
                particles_to_remove.push_back(particle);
            }
        }
        removeParticlesFromBlock(old_block_particles, particles_to_remove);
    }
}
void updateParticle(std::vector<int>block_index, Particle& particle, SimulationData data){
    double cord, var;
    std::vector<double> pos{0.0, 0.0, 0.0}, v{0.0, 0.0, 0.0}, hv{0.0, 0.0, 0.0};
    double minVar = pow(10,-10);
    std::vector<double> ngrids{data.grid.nx, data.grid.ny, data.grid.nz};
    if (particle.id == 2002) {
        cout << "---------------------------------" << '\n';
        cout << "Pos before update: "  << particle.px << ", " << particle.py << ", " << particle.pz << '\n';
        cout << "Hv before update: "  << particle.hvx << ", " << particle.hvy << ", " << particle.hvz << '\n';
        cout << "V before update: "  << particle.vx << ", " << particle.vy << ", " << particle.vz << '\n';
        cout << "Acc before update: "  << particle.acceleration[0] << ", " << particle.acceleration[1] << ", " << particle.acceleration[2] << '\n';
        cout << "---------------------------------" << '\n';
    }
    for(int i=0; i < 3 ; ++i){
        cord = calcCord(particle, i);
        var = calcVariation(block_index[i], cord, ngrids[i], i);
        if (var>minVar) {particle.acceleration[i] = calcAcceleration(particle, var, ngrids[i], i);}
        pos[i] = updatePosition(particle, i);
        v[i] = updateVelocity(particle, i);
        hv[i] = updateHv(particle, i);
        if (i==0){particle.px = pos[i],particle.vx = v[i],particle.hvx = hv[i];}
        else if (i==1){particle.py = pos[i],particle.vy = v[i],particle.hvy = hv[i];}
        else if (i==2){particle.pz = pos[i],particle.vz = v[i],particle.hvz = hv[i];}
        checkBorderLimits(particle, ngrids[i], i, block_index[i]);
    }
    if (particle.id == 2002) {
        cout << "Pos after update: "  << particle.px << ", " << particle.py << ", " << particle.pz << '\n';
        cout << "Hv after update: "  << particle.hvx << ", " << particle.hvy << ", " << particle.hvz << '\n';
        cout << "V after update: "  << particle.vx << ", " << particle.vy << ", " << particle.vz << '\n';
        cout << "Acc after update: "  << particle.acceleration[0] << ", " << particle.acceleration[1] << ", " << particle.acceleration[2] << '\n';
        cout << "---------------------------------" << '\n';
    }
}

void establishParticleFunctionality(std::map<std::vector<int>, std::shared_ptr<std::vector<Particle>>>& particleMap, SimulationData data){
    std::vector<int>block_index;
    std::shared_ptr<std::vector<Particle>> block_particles;
    for(auto block:particleMap){
        block_index = block.first;
        block_particles = block.second;
        for (Particle& particle: *block_particles){
            updateParticle(block_index, particle, data);
        }
    }
    updateParticleBlockBelonging(particleMap, data);
    cout << "here-4" << '\n';
}


int processSimulation(std::map<std::vector<int>, std::shared_ptr<std::vector<Particle>>>& particleMap, SimulationData& data){
    map<vector<int>,vector<vector<int>>> adjacent_blocks;
    adjacent_blocks = createAdjacentBlocks(data.grid);

    vector<int> current_block_key;
    //int c=0;
    //el updateDensityFalse se usará si decidimos dar prioridad a la memoria
    //updateDensityFalse(particleMap)
    std::cout << "before modify density" << '\n';
    for(auto current_block : particleMap){
        current_block_key = current_block.first;
        modifyBlock(current_block_key, particleMap, adjacent_blocks, data);
        //std::cout << c << '\n';
        //c+=1;
    }
    std::cout << "after modify density" << '\n';
    data.all_particles_density_updated = true;
    for(auto current_block: particleMap){
        current_block_key = current_block.first;
        modifyBlock(current_block_key, particleMap, adjacent_blocks, data);
    }
    data.all_particles_density_updated = false;
    establishParticleFunctionality(particleMap, data);
    initializeDensityAcceleration(particleMap);

    return 0;
}
