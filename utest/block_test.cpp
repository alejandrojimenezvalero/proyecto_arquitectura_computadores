#include <gtest/gtest.h>
#include "../sim/block.hpp"

TEST(CreateBlockTest, InitializesBlockCorrectly) {
    int i = 1, j = 2, k = 3;
    Block result = createBlock(i, j, k);

    EXPECT_EQ(result.block_index[0], i);
    EXPECT_EQ(result.block_index[1], j);
    EXPECT_EQ(result.block_index[2], k);
}