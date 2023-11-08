#ifndef FLUID_CONSTANTS_HPP
#define FLUID_CONSTANTS_HPP

#include <vector>
#include <cmath>
namespace simulationConstants {
  const float RADIO_MULTIPLICATOR = 1.695;
  const int FLUID_DENSITY         = 1000;
  const float STIFFNESS_PRESSURE = 3.0;
  const double STIFFNESS_COLLISIONS = 30000;
  const double DAMPING = 128.0;
  const float VISCOSITY = 0.4;
  const float PARTICLE_SIZE = 0.0002;
  const float TIMESTAMP = 0.001;
  const double PI = M_PI;
  const std::vector<double> GRAVITY = {0.0, -9.8, 0.0};
  const std::vector<float> UPPER_LIMIT{0.065, 0.1, 0.065};
  const std::vector<float> LOWER_LIMIT{-0.065,-0.08, -0.065};
}

#endif  // FLUID_CONSTANTS_HPP
