#ifndef FLUID_CONSTANTS_HPP
#define FLUID_CONSTANTS_HPP

#include <vector>
#include <cmath>
namespace simulationConstants {
  const double RADIO_MULTIPLICATOR = 1.695;
  const double FLUID_DENSITY         = 1000;
  const double STIFFNESS_PRESSURE = 3.0;
  const double STIFFNESS_COLLISIONS = 30000;
  const double DAMPING = 128.0;
  const double VISCOSITY = 0.4;
  const double PARTICLE_SIZE = 0.0002;
  const double TIMESTAMP = 0.001;
  const double TIMESTAMP_2 = pow(TIMESTAMP, 2);
  const double PI = M_PI;
  const double MINIMUM_VARIATION = pow(10,-10);
  const double MINIMUM_DISTANCE = pow(10, -12);
  const std::vector<double> GRAVITY = {0.0, -9.8, 0.0};
  const std::vector<double> UPPER_LIMIT{0.065, 0.1, 0.065};
  const std::vector<double> LOWER_LIMIT{-0.065,-0.08, -0.065};
}

#endif  // FLUID_CONSTANTS_HPP
