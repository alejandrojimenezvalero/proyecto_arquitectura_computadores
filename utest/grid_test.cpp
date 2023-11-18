//
#include "../sim/constants.hpp"
#include "../sim/grid.hpp"

#include "cmath"
#include "gtest/gtest.h"
using namespace simulationConstants;

TEST(CalculateBlockSizeTest, PositiveSizeGrid) {
  Grid grid {{10, 20, 30},{},{}, {}};
  calculateBlockSize(grid);
  double expected_sx = (UPPER_LIMIT[0] - LOWER_LIMIT[0]) / grid.grid_dimensions[0];
  double expected_sy = (UPPER_LIMIT[1] - LOWER_LIMIT[1]) / grid.grid_dimensions[1];
  double expected_sz = (UPPER_LIMIT[2] - LOWER_LIMIT[2]) / grid.grid_dimensions[2];

  EXPECT_DOUBLE_EQ(grid.block_dimensions[0], expected_sx);
  EXPECT_DOUBLE_EQ(grid.block_dimensions[1], expected_sy);
  EXPECT_DOUBLE_EQ(grid.block_dimensions[2], expected_sz);
}


TEST(CalculateBlockSizeTest, ZeroSizeGrid) {
  Grid grid {{0, 0, 0},{},{}, {}};
  calculateBlockSize(grid);

  EXPECT_DOUBLE_EQ(grid.block_dimensions[0], 0);
  EXPECT_DOUBLE_EQ(grid.block_dimensions[1], 0);
  EXPECT_DOUBLE_EQ(grid.block_dimensions[2], 0);
}


TEST(CalculateBlockSizeTest, NegativeSizeGrid) {
  Grid grid {{-5, -10, -15},{},{}, {}};
  calculateBlockSize(grid);

  EXPECT_DOUBLE_EQ(grid.block_dimensions[0], 0);
  EXPECT_DOUBLE_EQ(grid.block_dimensions[1], 0);
  EXPECT_DOUBLE_EQ(grid.block_dimensions[2], 0);
}


//Valid test with positive smoothing_length
TEST(InitGridSizeTest, PositiveSmoothingLength) {
  Grid grid {{},{},{}, {}};
  double smoothing_length = 1.0;
  initGrid(grid, smoothing_length);

  EXPECT_EQ(grid.grid_dimensions[0], static_cast<int>((UPPER_LIMIT[0] - LOWER_LIMIT[0]) / smoothing_length));
  EXPECT_EQ(grid.grid_dimensions[1], static_cast<int>((UPPER_LIMIT[1] - LOWER_LIMIT[1]) / smoothing_length));
  EXPECT_EQ(grid.grid_dimensions[2], static_cast<int>((UPPER_LIMIT[2] - LOWER_LIMIT[2]) / smoothing_length));
}
//Valid test with zero smoothing_length
TEST(InitGridSizeTest, ZeroSmoothingLength) {
  Grid grid {{},{},{}, {}};
  double smoothing_length = 0.0;
  initGrid(grid, smoothing_length);

  EXPECT_EQ(grid.grid_dimensions[0], 0);
  EXPECT_EQ(grid.grid_dimensions[1], 0);
  EXPECT_EQ(grid.grid_dimensions[2], 0);
}
//Valid test with negative smoothing_length
TEST(InitGridSizeTest, NegativeSmoothingLength) {
  Grid grid {{},{},{}, {}};
  double smoothing_length = -1.0;
  initGrid(grid, smoothing_length);

  EXPECT_EQ(grid.grid_dimensions[0], 0);
  EXPECT_EQ(grid.grid_dimensions[1], 0);
  EXPECT_EQ(grid.grid_dimensions[2], 0);
}

TEST(CalcParticleIndexTest, ParticleInsideGrid) {
  Particle particle;
  Grid grid {{10.0, 10.0, 10.0},{},{},{}};
  particle.pos = {20.0, 40.0, 50.0};

  std::vector<int> result = calcParticleIndex(particle, grid);

  EXPECT_EQ(result[0],2);
  EXPECT_EQ(result[1],4);
  EXPECT_EQ(result[2],5);
}

TEST(CalcParticleIndexTest, ParticleIndexLowerThanLimit) {
  Particle particle;
  Grid grid {{10.0, 10.0, 10.0},{},{},{}};
  particle.pos = {-0.075, -0.09, -0.075};

  std::vector<int> result = calcParticleIndex(particle, grid);

  EXPECT_EQ(result[0], 0);
  EXPECT_EQ(result[1], 0);
  EXPECT_EQ(result[2], 0);
}

TEST(CalcParticleIndexTest, ParticleIndexGreaterThanLimit) {
  Particle particle;
  Grid grid {{10.0, 10.0, 10.0},{},{},{}};
  particle.pos = {300.0, 300.0, 300.0};

  std::vector<int> result = calcParticleIndex(particle, grid);

  EXPECT_EQ(result[0], 9);
  EXPECT_EQ(result[1], 9);
  EXPECT_EQ(result[2], 9);
}