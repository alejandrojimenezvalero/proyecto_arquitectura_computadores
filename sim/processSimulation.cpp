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
std::map<std::vector<int>, std::vector<Particle>> createSubMap(std::map<std::vector<int>, std::vector<Particle>> particleMap, vector<vector<int>>adjacent_blocks){
  std::map<std::vector<int>, std::vector<Particle>> particleSubMap;
  vector<int> block_key;
  vector<Particle> block_adj_particles;
  for (auto block: particleMap){
    block_key = block.first;
    block_adj_particles = block.second;
    for (auto adj_block: adjacent_blocks){
      if(block_key == adj_block){
        particleSubMap[adj_block] = block_adj_particles;
      }
    }
  }
  return particleSubMap;
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

vector<double> transferAcceleration(Particle particlei, Particle particlej, double dist, SimulationData data){
  vector<double> variation_acc{0.0,0.0,0.0};
  double escalar_pos;
  escalar_pos = ((15/PI*pow(data.smoothing_length,6))*((3*data.particle_mass*STIFFNESS_PRESSURE)/2)*((pow(data.smoothing_length-dist,2))/dist)*(particlei.density+particlej.density-2*FLUID_DENSITY));

  return variation_acc;


}
int modifyDensity(Particle particle, std::map<std::vector<int>, std::vector<Particle>> particleMap, map<std::vector<int>,std::vector<std::vector<int>>> adjacent_blocks_map, SimulationData data) {

  // En el mapa de bloques adyacentes tenemos pares clave-valor de la forma bloque {i,j,k} = vector bloques adyacentes{{i,j,k}, {i,j,k},...}
  // Tomamos un vector de bloques adyacentes al bloque en el que está la partícula que se recibe como parámetro
  std::vector<vector<int>> adjacent_blocks = adjacent_blocks_map[{particle.i, particle.j, particle.k}];
  // Así tomamos un submapa del mapa de partículas de manera que nos quedamos con los pares {i,j,k} = {P1,P2,..} para todos los bloques adyacentes
  std::map<std::vector<int>, std::vector<Particle>> particleSubMap;
  particleSubMap = createSubMap(particleMap, adjacent_blocks);
  /*Teniendo el submapa con claves indices de bloques adyacentes al bloque de la partícula parámetro , y valores las partículas de dichos bloques adyacentes
   debemos iterar por los bloques del submapa, y para cada partícula de cada bloque adyacente calcular la densidad con respecto a la partícula parámetro*/
  double h = data.smoothing_length;
  double norm, dist;
  for(auto block: particleSubMap){
    vector<Particle> adj_particles = block.second;
    particle.density_updated = false, particle.acceleration_updated = false;
    for(auto adj_particle: adj_particles){
      if (!adj_particle.density_updated){
        norm = calculateNorm({particle.px,particle.py,particle.pz}, {adj_particle.px, adj_particle.py, adj_particle.pz});
        if(pow(norm,2)<pow(h,2)){
          particle.density += (pow(h,2)-pow(norm,2)), adj_particle.density += (pow(h,2)-pow(norm,2));
          particle.density = transformDensity(particle.density, data), adj_particle.density = transformDensity(adj_particle.density, data);
        }
      }
      if (!adj_particle.acceleration_updated){
        if(pow(norm,2)<pow(h,2)){
          dist = sqrt(max(pow(norm,2),pow(10,-12)));
          //particle.acceleration = transferAcceleration(particle,adj_particle, dist, data);
        }
      }
    }
    particle.density_updated = true, particle.acceleration_updated = true;




  }
  return 0;
}



/* (px-1, py-1, pz-1),
(px-1, py-1, pz),
(px-1, py-1, pz+1),
(px-1, py, pz-1),
(px-1, py, pz),
(px-1, py, pz+1),
(px-1, py+1, pz-1),
(px-1, py+1, pz),
(px-1, py+1, pz+1),
(px, py-1, pz-1),
(px, py-1, pz),
(px, py-1, pz+1),
(px, py, pz-1),
(px, py, pz+1),
(px, py+1, pz-1),
(px, py+1, pz),
(px, py+1, pz+1),
(px+1, py-1, pz-1),
(px+1, py-1, pz),
(px+1, py-1, pz+1),
(px+1, py, pz-1),
(px+1, py, pz),
(px+1, py, pz+1),
(px+1, py+1, pz-1),
(px+1, py+1, pz),
(px+1, py+1, pz+1). */