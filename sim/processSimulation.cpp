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

bool blockExists(int i, int j, int k, gridSize grid){
  if ((i < 0 and i > grid.nx-1) or (j < 0 and  grid.ny-1) or (k < 0 and grid.nz-1)){return false;}
  return true;
};

map<vector<int>,vector<vector<int>>> createAdjacentBlocks(gridSize grid) {
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

              if (blockExists(i, j, k, grid)) {
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
std::tuple<vector<Particle>, std::map<std::vector<int>, std::vector<Particle>>> createSubMap(std::vector<int> current_block_key, std::map<std::vector<int>, std::vector<Particle>> particleMap, vector<vector<int>> adjacent_blocks){
    std::map<std::vector<int>, std::vector<Particle>> particleSubMap;
    vector<int> block_key;
    vector<Particle> block_adj_particles;
    vector<Particle> current_block_particles;
    for (auto block: particleMap){
        block_key = block.first;
        block_adj_particles = block.second;
        if (block_key == current_block_key){
            current_block_particles = block_adj_particles;
        }
        for (auto adj_block: adjacent_blocks){
            if(block_key == adj_block){
                particleSubMap[adj_block] = block_adj_particles;
            }
        }
    }
    return std::make_tuple(current_block_particles, particleSubMap);
}

double calculateNorm(vector<double>particlei,vector<double>particlej){
  std::vector<double> dist{0.0,0.0,0.0};
  for (int i = 0; i < 2; i++) {
    dist[i] = particlei[i] - particlej[i];
  }
  double norm = std::sqrt(std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0));
  return norm;

}
double transformDensity(double density, SimulationData data){
  density = (density + pow(data.smoothing_length,6)) * (315/(64*PI*pow(data.smoothing_length,9))) * data.particle_mass;
  return density;
}

std::vector<double> transferAcceleration(Particle particlei, Particle particlej, double dist, SimulationData data){
  std::vector<double> variation_acc{};
  double escalar_pos, escalar_vel;
  vector<double> diff_pos = {particlei.px - particlej.px, particlei.py - particlej.py, particlei.pz - particlej.pz};
  vector<double> diff_vel = {particlej.vx - particlei.vx, particlej.vy - particlei.vy, particlej.vz - particlei.vz};
  escalar_pos = ((15/PI*pow(data.smoothing_length,6))*((3*data.particle_mass*STIFFNESS_PRESSURE)/2)*((pow(data.smoothing_length-dist,2))/dist)*(particlei.density+particlej.density-2*FLUID_DENSITY));
  escalar_vel = (45/PI*pow(data.smoothing_length,6))*VISCOSITY*data.particle_mass;
  variation_acc = {(diff_pos[0]*(escalar_pos)+diff_vel[0]*(escalar_vel))/(particlei.density*particlej.density),
                   (diff_pos[1]*(escalar_pos)+diff_vel[1]*(escalar_vel))/(particlei.density*particlej.density),
                   (diff_pos[2]*(escalar_pos)+diff_vel[2]*(escalar_vel))/(particlei.density*particlej.density)};
  return variation_acc;
}
int updateDensityFalse(std::map<std::vector<int>, std::vector<Particle>> particleMap){
    std::vector<int> block_key;
    std::vector<Particle> particles_block;
    for(auto block: particleMap){
        block_key = block.first;
        particles_block = block.second;
        for(Particle particle: particles_block){
            particle.density_updated = false;
        }
    }
    return 0;
}
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
void updateBlock2(std::vector<Particle> current_block_particles, std::map<std::vector<int>, std::vector<Particle>> particleSubMap, SimulationData data, string mode){
    // REFACTORIZAR CON VECTORES
    double h = data.smoothing_length;
    double norm, dist;
    //check_iter=0-> USAMOS density_updated, check_iter=1-> USAMOS density_updated2
    //el using_update_adj, usará en la ADJUNTA density_updated si estamos en check_iter=0, y usará density_updated2 si estamos en check_iter=1
    bool check_iter;bool using_update_adj;
    for(Particle particle: current_block_particles){
        for(auto block: particleSubMap){
            vector<Particle> adj_particles = block.second;
            check_iter = (mode=="density" && particle.density_updated==true)?1:0;
            for(auto adj_particle: adj_particles){
                using_update_adj = (check_iter==0)?adj_particle.density_updated:adj_particle.density_updated2;
                if (mode =="density"){
                    if ((!using_update_adj) && particle != adj_particle){
                        norm = calculateNorm({particle.px,particle.py,particle.pz}, {adj_particle.px, adj_particle.py, adj_particle.pz});
                        if(pow(norm,2)<pow(h,2)){
                            particle.density += (pow(h,2)-pow(norm,2)), adj_particle.density += (pow(h,2)-pow(norm,2));
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
            particle.density_updated = (mode=="density" && check_iter==0)?true:false;
            particle.density_updated2 = (mode=="density" && check_iter==1)?true:false;
            particle.acceleration_updated = (mode=="acceleration")?true:false;
        }
    }
}
int modifyBlock(std::vector<int> current_block_key, std::map<std::vector<int>, std::vector<Particle>> particleMap, map<std::vector<int>,std::vector<std::vector<int>>> adjacent_blocks_map, SimulationData data) {

  // En el mapa de bloques adyacentes tenemos pares clave-valor de la forma bloque {i,j,k} = vector bloques adyacentes{{i,j,k}, {i,j,k},...}
  // Tomamos un vector de bloques adyacentes al bloque en el que está la partícula que se recibe como parámetro
  std::vector<vector<int>> adjacent_blocks = adjacent_blocks_map[current_block_key];
  // Así tomamos un submapa del mapa de partículas de manera que nos quedamos con los pares {i,j,k} = {P1,P2,..} para todos los bloques adyacentes
  std::map<std::vector<int>, std::vector<Particle>> particleSubMap;
  std::vector<Particle> current_block_particles;
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


double calcCord(Particle particle, int index){
    double cord;
    cord = (index == 0)?particle.px + particle.hvx*TIMESTAMP :(index == 1)?particle.py + particle.hvy*TIMESTAMP :particle.pz + particle.hvz*TIMESTAMP;
    return cord;
}
double calcVariation(Particle particle, double cordParticle, double ngrid, int index){
    double var;
    int cordBlock;
    cordBlock = (index == 0)? particle.i : (index == 1)? particle.j : particle.k;
    if(cordBlock == 0){
        var = PARTICLE_SIZE - (cordParticle - LOWER_LIMIT[index]);
    }
    else if(cordBlock == ngrid - 1){
        var = PARTICLE_SIZE - (UPPER_LIMIT[index] - cordParticle);
    }
    return var;
}
double calcAcceleration(Particle particle, double var, double ngrid, int index){
    double v, a;
    int cordBlock;
    a = particle.acceleration[index];
    v = (index == 0)? particle.vx : (index == 1)? particle.vy : particle.vz;
    cordBlock = (index == 0)? particle.i : (index == 1)? particle.j : particle.k;
    if(cordBlock == 0){
        a += (STIFFNESS_COLLISIONS*var - DAMPING*v);
    }
    else if(cordBlock == ngrid - 1){
        a -= (STIFFNESS_COLLISIONS*var + DAMPING*v);
    }
    return a;
}
double updatePosition(Particle particle, int index, vector<int> block){
    double p, hv, a;
    p = (index == 0)? particle.px : (index == 1)? particle.py : particle.pz;
    hv = (index == 0)? particle.hvx : (index == 1)? particle.hvy : particle.hvz;
    a = particle.acceleration[index];
    p += hv*TIMESTAMP + a*pow(TIMESTAMP, 2);
    particle.i = (index==0)?std::floor((particle.px - LOWER_LIMIT[0])/block[0]):particle.i;
    particle.j = (index==1)?std::floor((particle.py - LOWER_LIMIT[1])/block[1]):particle.j;
    particle.k = (index==2)?std::floor((particle.pz - LOWER_LIMIT[2])/block[2]):particle.k;
    return p;
}
double updateVelocity(Particle particle, int index){
    double v, hv, a;
    hv = (index == 0)? particle.hvx : (index == 1)? particle.hvy : particle.hvz;
    a = particle.acceleration[index];
    v = hv + (a*TIMESTAMP)/2;
    return v;
}
double updateHv(Particle particle, int index){
    double hv, a;
    hv = (index == 0)? particle.hvx : (index == 1)? particle.hvy : particle.hvz;
    a = particle.acceleration[index];

    hv += a*TIMESTAMP;

    return hv;
}
double checkBorderLimits(Particle particle, double ngrid, int index){
    double d, p;
    int cordBlock;
    p = (index == 0)? particle.px : (index == 1)? particle.py : particle.pz;
    cordBlock = (index == 0)? particle.i : (index == 1)? particle.j : particle.k;
    if (cordBlock==0){
        d = p - LOWER_LIMIT[index];
    }
    else if (cordBlock == ngrid - 1){
        d = UPPER_LIMIT[index] - p;
    }
    if (d<0){
        double particleAxisPos = (cordBlock==0)?LOWER_LIMIT[index] - d: (cordBlock == ngrid - 1)?UPPER_LIMIT[index] + d : p;
        if (index==0){particle.px = particleAxisPos, particle.vx = -particle.vx, particle.hvx=-particle.hvx;}
        else if (index==1){particle.py = particleAxisPos, particle.vy = -particle.vy, particle.hvy=-particle.hvy;}
        else if (index==2){particle.pz = particleAxisPos, particle.vz = -particle.vz, particle.hvz=-particle.hvz;}
    }
    return d;
}
void establishParticleFunctionality(std::map<std::vector<int>, std::vector<Particle>> particleMap, gridSize grid){
    std::vector<int>block_index;
    std::vector<Particle> block_particles;
    std::vector<double> cords{}, vars{};
    std::vector<double> pos{}, v{}, hv{};
    std::vector<double> ngrids{grid.nx, grid.ny, grid.nz};
    double minVar = pow(10,-10);
    for(auto block:particleMap){
        block_index = block.first;
        block_particles = block.second;
        for (Particle particle: block_particles){
            for(int i=0; i < 3 ; ++i){
                cords[i] = calcCord(particle, i);
                vars[i] = calcVariation(particle, cords[i], ngrids[i], i);
                particle.acceleration[i] = (vars[i]>minVar) ? calcAcceleration(particle, vars[i], ngrids[i], i) : particle.acceleration[i];
                pos[i] = updatePosition(particle, i, block_index);
                v[i] = updateVelocity(particle, i);
                hv[i] = updateHv(particle, i);
                if (i==0){particle.px = pos[i],particle.vx = v[i],particle.hvx = hv[i];}
                else if (i==1){particle.py = pos[i],particle.vy = v[i],particle.hvy = hv[i];}
                else if (i==2){particle.pz = pos[i],particle.vz = v[i],particle.hvz = hv[i];}
                checkBorderLimits(particle, ngrids[i], i);
            }
        }
    }
}

int processSimulation(std::map<std::vector<int>, std::vector<Particle>> particleMap, SimulationData data){
    map<vector<int>,vector<vector<int>>> adjacent_blocks;
    std::cout << "before adjacents" << '\n';
    adjacent_blocks = createAdjacentBlocks(data.grid);
    std::cout << "after adjacents" << '\n';
    vector<int> current_block_key;
    int c=0;
    //el updateDensityFalse se usará si decidimos dar prioridad a la memoria
    //updateDensityFalse(particleMap)
    for(auto current_block: particleMap){
        current_block_key = current_block.first;
        modifyBlock(current_block_key, particleMap, adjacent_blocks, data);
        std::cout << c << '\n';
        c+=1;
    }
    std::cout << "after modify density" << '\n';
    data.all_particles_density_updated = true;
    for(auto current_block: particleMap){
        current_block_key = current_block.first;
        modifyBlock(current_block_key, particleMap, adjacent_blocks, data);
    }
    data.all_particles_density_updated = false;
    establishParticleFunctionality(particleMap, data.grid);
    return 0;
}
