//
#include "../sim/constants.hpp"
#include "../sim/grid.hpp"

#include "cmath"
#include "gtest/gtest.h"
using namespace simulationConstants;

TEST(initGridTest, PositiveSmoothinglength) {
  Grid grid {{},{},{}, {},};
  double smoothing_length =  0.00830882;
  initGrid(grid, smoothing_length);

  EXPECT_DOUBLE_EQ(grid.grid_dimensions[0],15.0);
  EXPECT_DOUBLE_EQ(grid.grid_dimensions[1],21.0);
  EXPECT_DOUBLE_EQ(grid.grid_dimensions[2],15.0);

  double s_x = (UPPER_LIMIT[0] - LOWER_LIMIT[0])/grid.grid_dimensions[0];
  double s_y = (UPPER_LIMIT[1] - LOWER_LIMIT[1])/grid.grid_dimensions[1];
  double s_z = (UPPER_LIMIT[2] - LOWER_LIMIT[2])/grid.grid_dimensions[2];

  EXPECT_DOUBLE_EQ(grid.block_dimensions[0], s_x);
  EXPECT_DOUBLE_EQ(grid.block_dimensions[1], s_y);
  EXPECT_DOUBLE_EQ(grid.block_dimensions[2], s_z);
}


TEST(initGridTest, NegativeSmoothinglength) {
  Grid grid {{},{},{}, {},};
  double smoothing_length =  -0.00830882;
  initGrid(grid, smoothing_length);

  EXPECT_DOUBLE_EQ(grid.grid_dimensions[0],0.0);
  EXPECT_DOUBLE_EQ(grid.grid_dimensions[1],0.0);
  EXPECT_DOUBLE_EQ(grid.grid_dimensions[2],0.0);

  EXPECT_DOUBLE_EQ(grid.block_dimensions[0], 0.0);
  EXPECT_DOUBLE_EQ(grid.block_dimensions[1], 0.0);
  EXPECT_DOUBLE_EQ(grid.block_dimensions[2], 0.0);
}

TEST(initGridTest, ZeroSmoothinglength) {
  Grid grid {{},{},{}, {},};
  double smoothing_length =  0.0;
  initGrid(grid, smoothing_length);

  EXPECT_DOUBLE_EQ(grid.grid_dimensions[0],0.0);
  EXPECT_DOUBLE_EQ(grid.grid_dimensions[1],0.0);
  EXPECT_DOUBLE_EQ(grid.grid_dimensions[2],0.0);

  EXPECT_DOUBLE_EQ(grid.block_dimensions[0], 0.0);
  EXPECT_DOUBLE_EQ(grid.block_dimensions[1], 0.0);
  EXPECT_DOUBLE_EQ(grid.block_dimensions[2], 0.0);
}


TEST(CalcParticleIndexTest, ParticleInsideGrid) {
  Particle particle;
  Grid grid {{15.0, 21.0, 15.0},{0.00866667,0.00857143,0.00866667},{},{}};
  particle.pos = {20.0, 40.0, 50.0};

  std::vector<int> result = calcParticleIndex(particle, grid);



  // Ahora puedes realizar las comparaciones
  //EXPECT_EQ(result[0], 2);
  //EXPECT_EQ(result[1], 4);
  //EXPECT_EQ(result[2], 5);
}

TEST(ParticleIndexCalculationTest, TestParticleInsideLowerLimit) {
  Particle particle;
  Grid grid {{15.0, 21.0, 15.0},{0.00866667,0.00857143,0.00866667},{},{}};

  // Configurar valores específicos para particle.pos, grid.block_dimensions, etc.
  particle.pos = {-0.04, -0.05, -0.04};

  std::vector<int> result = calcParticleIndex(particle, grid);


  EXPECT_EQ(result[0], 2);
  EXPECT_EQ(result[1], 3);
  EXPECT_EQ(result[2], 2);
}

TEST(ParticleIndexCalculationTest, TestParticleOutsideXLowerLimit) {
  Particle particle;
  Grid grid{{15.0, 21.0, 15.0},{0.00866667,0.00857143,0.00866667},{},{}};

  particle.pos = {-0.07, -0.07, -0.05};
  std::vector<int> result = calcParticleIndex(particle, grid);
  EXPECT_EQ(result[0], 0);
  EXPECT_EQ(result[1], 1);
  EXPECT_EQ(result[2], 1);

}

TEST(ParticleIndexCalculationTest, TestParticleOutsideYLowerLimit) {
  Particle particle;
  Grid grid{{15.0, 21.0, 15.0},{0.00866667,0.00857143,0.00866667},{},{}};

  particle.pos = {-0.05, -0.09, -0.05};
  std::vector<int> result = calcParticleIndex(particle, grid);
  EXPECT_EQ(result[0], 1);
  EXPECT_EQ(result[1], 0);
  EXPECT_EQ(result[2], 1);

}
TEST(ParticleIndexCalculationTest, TestParticleOutsideZLowerLimit) {
  Particle particle;
  Grid grid{{15.0, 21.0, 15.0},{0.00866667,0.00857143,0.00866667},{},{}};

  particle.pos = {-0.05, -0.07, -0.09};
  std::vector<int> result = calcParticleIndex(particle, grid);
  EXPECT_EQ(result[0], 1);
  EXPECT_EQ(result[1], 1);
  EXPECT_EQ(result[2], 0);

}


TEST(blockExistsTest, IndexesInRange) {
  Grid grid {{10.0, 10.0, 10.0},{},{},{}};
  int i = 1; int j = 1; int k = 1;
  bool result = blockExists(i, j, k, grid);

  ASSERT_TRUE(result);
}

TEST(blockExistsTest, iNotInRange) {
  Grid grid {{10.0, 10.0, 10.0},{},{},{}};
  int i = -2; int j = 1; int k = 1;
  bool result = blockExists(i, j, k, grid);

  ASSERT_FALSE(result);
}

TEST(blockExistsTest, jNotInRange) {
  Grid grid {{10.0, 10.0, 10.0},{},{},{}};
  int i = 1; int j = -2; int k = 1;
  bool result = blockExists(i, j, k, grid);

  ASSERT_FALSE(result);
}

TEST(blockExistsTest, kNotInRange) {
  Grid grid {{10.0, 10.0, 10.0},{},{},{}};
  int i = 1; int j = 1; int k = -2;
  bool result = blockExists(i, j, k, grid);

  ASSERT_FALSE(result);
}