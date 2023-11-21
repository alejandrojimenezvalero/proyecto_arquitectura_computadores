#include <gtest/gtest.h>
#include "../sim/block.hpp"

TEST(CreateBlockTest, InitializesBlockCorrectly) {
    const int x_axis = 1;
    const int y_axis = 2;
    const int z_axis = 3;
    Block result = createBlock(x_axis, y_axis, z_axis);

    EXPECT_EQ(result.block_index[0], x_axis);
    EXPECT_EQ(result.block_index[1], y_axis);
    EXPECT_EQ(result.block_index[2], z_axis);
}