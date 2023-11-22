//

#include "particle.hpp"
#include "grid.hpp"
#include "block.hpp"
#include "initSimulation.hpp"
#include "constants.hpp"


#include <cmath>
#include <algorithm>
#include <tuple>

using namespace std;
using namespace simulationConstants;


double calculateNorm(const std::vector<double> &particlei, const std::vector<double> &particlej) {
    double sumOfSquares = 0.0;
    for (int i = 0; i < 3; i++) {
        sumOfSquares += (particlei[i] - particlej[i]) * (particlei[i] - particlej[i]);
    }
    return sumOfSquares;
}

void transformDensity(Particle &particle, const double &smoothing_length_6, const double &escalar_density) {
    particle.density = (particle.density + smoothing_length_6) * escalar_density;
}

void transferAcceleration(Particle &particlei, Particle &particlej, const double &dist, const SimulationData &data) {
    const double escalar_pos = data.escalar_pos *
                               ((std::pow(data.smoothing_length - dist, 2)) / dist) *
                               (particlei.density + particlej.density - 2 * simulationConstants::FLUID_DENSITY);
    const double escalar_vel = data.escalar_vel * simulationConstants::VISCOSITY * data.particle_mass;
    const double density_product = particlei.density * particlej.density;

    std::vector<double> diff_pos = {particlei.pos[0] - particlej.pos[0], particlei.pos[1] - particlej.pos[1],
                                    particlei.pos[2] - particlej.pos[2]};
    std::vector<double> diff_vel = {particlej.vel[0] - particlei.vel[0], particlej.vel[1] - particlei.vel[1],
                                    particlej.vel[2] - particlei.vel[2]};

    std::vector<double> variation_acc = {
            (diff_pos[0] * escalar_pos + diff_vel[0] * escalar_vel) / (density_product),
            (diff_pos[1] * escalar_pos + diff_vel[1] * escalar_vel) / (density_product),
            (diff_pos[2] * escalar_pos + diff_vel[2] * escalar_vel) / (density_product)
    };

    particlei.acceleration = {particlei.acceleration[0] + variation_acc[0],
                              particlei.acceleration[1] + variation_acc[1],
                              particlei.acceleration[2] + variation_acc[2]};
    particlej.acceleration = {particlej.acceleration[0] - variation_acc[0],
                              particlej.acceleration[1] - variation_acc[1],
                              particlej.acceleration[2] - variation_acc[2]};
}

void updateBlocksDensity(SimulationData &data) {
    for (Block &block: data.grid.grid_blocks) {
        for (Particle &particle: block.block_particles_v2) {
            for (const int index_adj: block.adj_index_vector_blocks) {
                Block &adj_block = data.grid.grid_blocks[index_adj];
                if (!adj_block.updated_density) {
                    for (Particle &adj_particle: adj_block.block_particles_v2) {
                        if (!adj_particle.density_updated && particle.id != adj_particle.id) {
                            const double norm_2 = calculateNorm(particle.pos, adj_particle.pos);
                            if (norm_2 < data.smoothing_length_2) {
                                const double added_norm = std::pow((data.smoothing_length_2 - norm_2), 3);
                                particle.density += added_norm, adj_particle.density += added_norm;
                            }
                        }
                    }
                }

            }
            transformDensity(particle, data.smoothing_length_6, data.escalar_density);
            particle.density_updated = true;
        }
        block.updated_density = true;
    }
}

void updateBlocksAcceleration(SimulationData &data) {
    for (Block &block: data.grid.grid_blocks) {
        for (Particle &particle: block.block_particles_v2) {
            for (const int index_adj: block.adj_index_vector_blocks) {
                Block &adj_block = data.grid.grid_blocks[index_adj];
                if (!adj_block.updated_acceleration) {
                    for (Particle &adj_particle: adj_block.block_particles_v2) {
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
        }
        block.updated_acceleration = true;
    }
}

double calcCord(Particle &particle, int index) {
    const double cord = particle.pos[index] + particle.hv[index] * TIMESTAMP;
    return cord;
}

double calcVariation(int &index_block, const double &cordParticle, const double &grid_dimensions, int &index) {
    double var = 0.0;
    if (index_block == 0) {
        var = PARTICLE_SIZE - (cordParticle - LOWER_LIMIT[index]);
    } else if (index_block == grid_dimensions - 1) {
        var = PARTICLE_SIZE - (UPPER_LIMIT[index] - cordParticle);
    }
    return var;
}

void calcAcceleration(Particle &particle, const double &var, SimulationData &data, int &index) {
    const double axis_particle_velocity = particle.vel[index];
    const int cordBlock = calcParticleIndex(particle, data.grid)[index];
    if (cordBlock == 0) {
        particle.acceleration[index] += (STIFFNESS_COLLISIONS * var - DAMPING * axis_particle_velocity);
    } else if (cordBlock == data.grid.grid_dimensions[index] - 1) {
        particle.acceleration[index] -= (STIFFNESS_COLLISIONS * var + DAMPING * axis_particle_velocity);
    }
}

void updatePosition(Particle &particle, int &index) {
    particle.pos[index] += particle.hv[index] * TIMESTAMP + particle.acceleration[index] * TIMESTAMP_2;
}

void updateVelocity(Particle &particle, int &index) {
    particle.vel[index] = particle.hv[index] + (particle.acceleration[index] * TIMESTAMP) / 2;
}

void updateHv(Particle &particle, int &index) {
    particle.hv[index] += particle.acceleration[index] * TIMESTAMP;
}

void checkBorderLimits(Particle &particle, const double &grid_dimensions, int index, int actual_block_index) {
    double d_variable = 0.0;
    const double lower_limit = LOWER_LIMIT[index];
    const double upper_limit = UPPER_LIMIT[index];
    if (actual_block_index == 0) {
        d_variable = particle.pos[index] - lower_limit;
    } else if (actual_block_index == grid_dimensions - 1) {
        d_variable = upper_limit - particle.pos[index];
    }
    if (d_variable < 0) {
        double particleAxisPos = particle.pos[index];
        if (actual_block_index == 0) {
            particleAxisPos = lower_limit - d_variable;
        } else if (actual_block_index == grid_dimensions - 1) {
            particleAxisPos = upper_limit + d_variable;
        }
        particle.pos[index] = particleAxisPos;
        particle.vel[index] = -particle.vel[index];
        particle.hv[index] = -particle.hv[index];
    }
}

void updateParticleBlockBelonging(SimulationData &data) {
    for (std::size_t i = 0; i < data.grid.grid_blocks.size(); i += 1) {
        Block &current_block = data.grid.grid_blocks[i];
        const vector<int> &current_block_index = current_block.block_index;
        for (Particle &particle: current_block.block_particles_v1) {
            const vector<int> &new_block_particle_index = calcParticleIndex(particle, data.grid);
            particle.density = 0.0;
            particle.density_updated = false;
            particle.acceleration = GRAVITY;
            particle.acceleration_updated = false;
            if (new_block_particle_index != current_block_index) {
                const int index_in_vector = calcParticleIndexVector(data.grid, new_block_particle_index);
                data.grid.grid_blocks[index_in_vector].block_particles_v2.push_back(particle);
            } else {
                current_block.block_particles_v2.push_back(particle);
            }
        }
        current_block.updated_density = false;
        current_block.updated_acceleration = false;
        current_block.block_particles_v1.clear();
    }
}

void updateParticle(std::vector<int> &block_index, Particle &particle, SimulationData &data) {
    for (int i = 0; i < 3; ++i) {
        const double cord = calcCord(particle, i);
        const double var = calcVariation(block_index[i], cord, data.grid.grid_dimensions[i], i);
        if (var > MINIMUM_VARIATION) { calcAcceleration(particle, var, data, i); }
        updatePosition(particle, i);
        updateVelocity(particle, i);
        updateHv(particle, i);
        checkBorderLimits(particle, data.grid.grid_dimensions[i], i, block_index[i]);
    }
    // Las propiedades asociadas a una particula al final de una iteración
    /*if (particle.id == 2229) {
        cout << "Density after limit process: " << particle.density << '\n';
        cout << "Pos after limit process: "  << particle.pos[0] << ", " << particle.pos[1] << ", " << particle.pos[2] << '\n';
        cout << "Hv after limit process: "  << particle.hv[0] << ", " << particle.hv[1] << ", " << particle.hv[2] << '\n';
        cout << "V after limit process: "  << particle.vel[0] << ", " << particle.vel[1] << ", " << particle.vel[2] << '\n';
        cout << "Acc after limit process: "  << particle.acceleration[0] << ", " << particle.acceleration[1] << ", " << particle.acceleration[2] << '\n';
        cout << "---------------------------------" << '\n';
    }*/
    const int index_in_vector = calcParticleIndexVector(data.grid, block_index);
    data.grid.grid_blocks[index_in_vector].block_particles_v1.push_back(particle);
}

void establishParticleFunctionality(SimulationData &data) {
    for (Block &block: data.grid.grid_blocks) {
        for (Particle &particle: block.block_particles_v2) {
            updateParticle(block.block_index, particle, data);
        }
        block.block_particles_v2.clear();
    }
}

int processSimulation(SimulationData &data) {
    updateParticleBlockBelonging(data);
    updateBlocksDensity(data);
    updateBlocksAcceleration(data);
    establishParticleFunctionality(data);
    return 0;
}
