//

#include "gtest/gtest.h"
#include "../sim/constants.hpp"
#include "../sim/grid.hpp"
#include "../sim/block.hpp"


using namespace simulationConstants;

// Test: Cálculo de tamaño de bloque con una cuadrícula no trivial
TEST(CalculateBlockSizeTest, NonTrivialGrid) {
  Grid grid {10, 20, 30};
  std::vector<double> result = calculateBlockSize(grid);

  double expected_sx = (UPPER_LIMIT[0] - LOWER_LIMIT[0]) / grid.nx;
  double expected_sy = (UPPER_LIMIT[1] - LOWER_LIMIT[1]) / grid.ny;
  double expected_sz = (UPPER_LIMIT[2] - LOWER_LIMIT[2]) / grid.nz;

  EXPECT_DOUBLE_EQ(result[0], expected_sx);
  EXPECT_DOUBLE_EQ(result[1], expected_sy);
  EXPECT_DOUBLE_EQ(result[2], expected_sz);
}

// Test: Cálculo de tamaño de bloque con una cuadrícula de tamaño cero
TEST(CalculateBlockSizeTest, ZeroSizeGrid) {
  Grid grid {0, 0, 0};
  std::vector<double> result = calculateBlockSize(grid);

  // En este caso, los resultados deberían ser cero
  EXPECT_DOUBLE_EQ(result[0], 0.0);
  EXPECT_DOUBLE_EQ(result[1], 0.0);
  EXPECT_DOUBLE_EQ(result[2], 0.0);
}

// Test: Cálculo de tamaño de bloque con una cuadrícula de tamaño negativo
TEST(CalculateBlockSizeTest, NegativeSizeGrid) {
  Grid grid {-5, -10, -15};
  std::vector<double> result = calculateBlockSize(grid);

  // En este caso, los resultados deberían ser cero
  EXPECT_DOUBLE_EQ(result[0], 0.0);
  EXPECT_DOUBLE_EQ(result[1], 0.0);
  EXPECT_DOUBLE_EQ(result[2], 0.0);
}

// Test: Cálculo de tamaño de bloque con una cuadrícula de tamaño positivo
TEST(CalculateBlockSizeTest, PositiveSizeGrid) {
  Grid grid {5, 10, 15};
  std::vector<double> result = calculateBlockSize(grid);

  // Calcula los valores esperados
  double expected_sx = (UPPER_LIMIT[0] - LOWER_LIMIT[0]) / grid.nx;
  double expected_sy = (UPPER_LIMIT[1] - LOWER_LIMIT[1]) / grid.ny;
  double expected_sz = (UPPER_LIMIT[2] - LOWER_LIMIT[2]) / grid.nz;

  EXPECT_DOUBLE_EQ(result[0], expected_sx);
  EXPECT_DOUBLE_EQ(result[1], expected_sy);
  EXPECT_DOUBLE_EQ(result[2], expected_sz);
}
