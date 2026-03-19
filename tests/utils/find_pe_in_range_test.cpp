#include "kagen/tools/utils.h"

#include <gtest/gtest.h>
#include <mpi.h>

using namespace kagen;

static std::vector<VertexRange> MakeContiguousRanges(std::initializer_list<SInt> bounds) {
    std::vector<SInt>        b(bounds);
    std::vector<VertexRange> ranges;
    for (std::size_t i = 0; i + 1 < b.size(); ++i) {
        ranges.emplace_back(b[i], b[i + 1]);
    }
    return ranges;
}

class FindPEInRangeTest : public ::testing::Test {
protected:
    // Ranges: [0,5), [5,10), [10,15)
    std::vector<VertexRange> ranges = MakeContiguousRanges({0, 5, 10, 15});
};

TEST_F(FindPEInRangeTest, FindsCorrectPE) {
    EXPECT_EQ(FindPEInRange(0, ranges), 0);
    EXPECT_EQ(FindPEInRange(4, ranges), 0);
    EXPECT_EQ(FindPEInRange(5, ranges), 1);
    EXPECT_EQ(FindPEInRange(9, ranges), 1);
    EXPECT_EQ(FindPEInRange(10, ranges), 2);
    EXPECT_EQ(FindPEInRange(14, ranges), 2);
}

TEST_F(FindPEInRangeTest, ReturnsMinusOneForOutOfRange) {
    EXPECT_EQ(FindPEInRange(15, ranges), -1);
    EXPECT_EQ(FindPEInRange(100, ranges), -1);
}

TEST_F(FindPEInRangeTest, BinarySearchMatchesLinear) {
    for (SInt node = 0; node <= 20; ++node) {
        EXPECT_EQ(FindPEInRangeWithBinarySearch(node, ranges), FindPEInRange(node, ranges))
            << "Mismatch for node " << node;
    }
}

TEST(FindPEInRangeBinarySearchTest, SingleRange) {
    std::vector<VertexRange> ranges = {{0, 100}};
    EXPECT_EQ(FindPEInRangeWithBinarySearch(0, ranges), 0);
    EXPECT_EQ(FindPEInRangeWithBinarySearch(50, ranges), 0);
    EXPECT_EQ(FindPEInRangeWithBinarySearch(99, ranges), 0);
    EXPECT_EQ(FindPEInRangeWithBinarySearch(100, ranges), -1);
}

TEST(FindPEInRangeBinarySearchTest, EmptyAndNonEmptyRanges) {
    std::vector<VertexRange> ranges = {{0, 1}, {1, 1}, {1, 1}, {1, 100}, {100, 101}};
    EXPECT_EQ(FindPEInRangeWithBinarySearch(0, ranges), 0);
    EXPECT_EQ(FindPEInRangeWithBinarySearch(1, ranges), 3);
    EXPECT_EQ(FindPEInRangeWithBinarySearch(99, ranges), 3);
    EXPECT_EQ(FindPEInRangeWithBinarySearch(100, ranges), 4);
}

TEST(FindPEInRangeLinearTest, EmptyRanges) {
    std::vector<VertexRange> ranges = {{0, 1}, {1, 1}, {1, 1}, {1, 100}, {100, 101}};
    EXPECT_EQ(FindPEInRange(0, ranges), 0);
    EXPECT_EQ(FindPEInRange(1, ranges), 3);
    EXPECT_EQ(FindPEInRange(99, ranges), 3);
    EXPECT_EQ(FindPEInRange(100, ranges), 4);
}

TEST(FindPEInRangeBinarySearchTest, ManyRanges) {
    // 100 ranges: [0,10), [10,20), ..., [990,1000)
    std::vector<VertexRange> ranges;
    for (SInt i = 0; i < 100; ++i) {
        ranges.emplace_back(i * 10, (i + 1) * 10);
    }
    for (SInt node = 0; node <= 1005; ++node) {
        EXPECT_EQ(FindPEInRangeWithBinarySearch(node, ranges), FindPEInRange(node, ranges))
            << "Mismatch for node " << node;
    }
}
