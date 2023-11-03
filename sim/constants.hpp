#ifndef FLUID_CONSTANTS_HPP
#define FLUID_CONSTANTS_HPP

#include <vector>

namespace simulationConstants {
  const float RADIO_MULTIPLICATOR = 1.695;
  const int FLUID_DENSITY         = 1000;
  const float STIFFNESS_PRESSURE = 3.0;
  const int STIFFNESS_COLLISIONS = 30000;
  const float DAMPING = 128.0;
  const float VISCOSITY = 0.4;
  const float PARTICLE_SIZE = 0.0002;
  const float TIMESTAMP = 0.001;
  const std::vector<float> UPPER_LIMIT{0.065, 0.1, 0.065};
  const std::vector<float> LOWER_LIMIT{-0.065,-0.08, -0.065};
}

#endif  // FLUID_CONSTANTS_HPP
