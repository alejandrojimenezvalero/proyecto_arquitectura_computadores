//
#include "../sim/constants.hpp"
#include "../sim/grid.hpp"

#include "cmath"
#include "gtest/gtest.h"
using namespace simulationConstants;

//Valid test with positive smoothing_length
TEST(CalculateGridSizeTest, PositiveSmoothingLength) {
  double smoothing_length = 1.0;
  gridSize result = calculateGridSize(smoothing_length);
  EXPECT_EQ(result.nx, static_cast<int>((UPPER_LIMIT[0] - LOWER_LIMIT[0]) / smoothing_length));
  EXPECT_EQ(result.ny, static_cast<int>((UPPER_LIMIT[1] - LOWER_LIMIT[1]) / smoothing_length));
  EXPECT_EQ(result.nz, static_cast<int>((UPPER_LIMIT[2] - LOWER_LIMIT[2]) / smoothing_length));
}
//Valid test with zero smoothing_length
TEST(CalculateGridSizeTest, ZeroSmoothingLength) {
  double smoothing_length = 0.0;
  gridSize result = calculateGridSize(smoothing_length);
  EXPECT_EQ(result.nx, 0);
  EXPECT_EQ(result.ny, 0);
  EXPECT_EQ(result.nz, 0);
}
//Valid test with negative smoothing_length
TEST(CalculateGridSizeTest, NegativeSmoothingLength) {
  double smoothing_length = -1.0;
  gridSize result = calculateGridSize(smoothing_length);
  EXPECT_EQ(result.nx, 0);
  EXPECT_EQ(result.ny, 0);
  EXPECT_EQ(result.nz, 0);
}
