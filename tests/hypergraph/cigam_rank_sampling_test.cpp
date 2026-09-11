#include "kagen/generators/hyper/h_erdos/hyper_er_common.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

namespace kagen {
namespace {

TEST(CIGAMRankSamplingTest, RejectsInvalidSizesAndPositions) {
    EXPECT_THROW(SortedUniformOrderStatistics(0, 1), ConfigurationError);

    const SortedUniformOrderStatistics sample(8, 1);
    EXPECT_THROW(sample.Ascending(-1), ConfigurationError);
    EXPECT_THROW(sample.Ascending(8), ConfigurationError);
    EXPECT_THROW(sample.Descending(-1), ConfigurationError);
    EXPECT_THROW(sample.Descending(8), ConfigurationError);
}

TEST(CIGAMRankSamplingTest, ProducesStrictlyOrderedBoundedUniforms) {
    constexpr SInt                     n = 257;
    const SortedUniformOrderStatistics sample(n, 17);

    long double previous_ascending  = 0.0L;
    long double previous_descending = 1.0L;

    for (SInt position = 0; position < n; ++position) {
        const long double ascending  = sample.Ascending(position);
        const long double descending = sample.Descending(position);

        EXPECT_GT(ascending, 0.0L);
        EXPECT_LT(ascending, 1.0L);
        EXPECT_GT(descending, 0.0L);
        EXPECT_LT(descending, 1.0L);

        EXPECT_GT(ascending, previous_ascending) << "at ascending position " << position;
        EXPECT_LT(descending, previous_descending) << "at descending position " << position;

        previous_ascending  = ascending;
        previous_descending = descending;
    }
}

TEST(CIGAMRankSamplingTest, AscendingAndDescendingViewsAgree) {
    constexpr SInt                     n = 128;
    const SortedUniformOrderStatistics sample(n, 1234567);

    for (SInt position = 0; position < n; ++position) {
        EXPECT_EQ(sample.Descending(position), sample.Ascending(n - 1 - position));
    }
}

TEST(CIGAMRankSamplingTest, ReconstructionIsDeterministicAndStateless) {
    constexpr SInt                     n = 513;
    const SortedUniformOrderStatistics first(n, 42);
    const SortedUniformOrderStatistics second(n, 42);

    // Deliberately query in different orders. Reconstruction must be stateless.
    std::vector<SInt> positions{512, 0, 127, 400, 1, 256, 23};

    for (const SInt position: positions) {
        EXPECT_EQ(first.Ascending(position), second.Ascending(position));
    }

    std::reverse(positions.begin(), positions.end());

    for (const SInt position: positions) {
        EXPECT_EQ(first.Ascending(position), second.Ascending(position));
    }
}

TEST(CIGAMRankSamplingTest, ReconstructionDependsOnSeed) {
    constexpr SInt                     n = 64;
    const SortedUniformOrderStatistics first(n, 1);
    const SortedUniformOrderStatistics second(n, 2);

    bool found_difference = false;
    for (SInt position = 0; position < n; ++position) {
        found_difference |= first.Ascending(position) != second.Ascending(position);
    }

    EXPECT_TRUE(found_difference);
}

TEST(CIGAMRankSamplingTest, SampledOrderStatisticsHaveCorrectMeans) {
    constexpr SInt          n             = 7;
    constexpr std::uint64_t samples       = 4096;
    constexpr long double   absolute_slop = 0.015L;

    std::vector<long double> sum(static_cast<std::size_t>(n), 0.0L);

    for (std::uint64_t seed = 0; seed < samples; ++seed) {
        const SortedUniformOrderStatistics sample(n, seed);

        for (SInt position = 0; position < n; ++position) {
            sum[static_cast<std::size_t>(position)] += sample.Ascending(position);
        }
    }

    for (SInt position = 0; position < n; ++position) {
        const long double observed = sum[static_cast<std::size_t>(position)] / samples;
        const long double expected = static_cast<long double>(position + 1) / static_cast<long double>(n + 1);

        EXPECT_NEAR(static_cast<double>(observed), static_cast<double>(expected), static_cast<double>(absolute_slop))
            << "at ascending position " << position;
    }
}

TEST(CIGAMRankSamplingTest, BulkAscendingMatchesScalarReconstruction) {
    constexpr SInt                     n = 257;
    const SortedUniformOrderStatistics sample(n, 42);

    const std::vector<std::pair<SInt, SInt>> intervals{
        {0, 1}, {0, 64}, {17, 93}, {128, 257}, {256, 257},
    };

    for (const auto [begin, end]: intervals) {
        std::vector<long double> values;
        sample.FillAscending(begin, end, values);

        ASSERT_EQ(values.size(), static_cast<std::size_t>(end - begin));

        for (SInt position = begin; position < end; ++position) {
            EXPECT_EQ(values[static_cast<std::size_t>(position - begin)], sample.Ascending(position));
        }
    }
}

TEST(CIGAMRankSamplingTest, BulkDescendingMatchesScalarReconstruction) {
    constexpr SInt                     n = 257;
    const SortedUniformOrderStatistics sample(n, 42);

    std::vector<long double> values;
    sample.FillDescending(31, 201, values);

    ASSERT_EQ(values.size(), std::size_t{170});

    for (SInt position = 31; position < 201; ++position) {
        EXPECT_EQ(values[static_cast<std::size_t>(position - 31)], sample.Descending(position));
    }
}

TEST(CIGAMRankSamplingTest, BulkLookupAcceptsAnEmptyInterval) {
    const SortedUniformOrderStatistics sample(32, 7);
    std::vector<long double>           values{1.0L};

    sample.FillAscending(12, 12, values);

    EXPECT_TRUE(values.empty());
}

TEST(CIGAMRankSamplingTest, BulkLookupRejectsInvalidIntervals) {
    const SortedUniformOrderStatistics sample(32, 7);
    std::vector<long double>           values;

    EXPECT_THROW(sample.FillAscending(-1, 2, values), ConfigurationError);
    EXPECT_THROW(sample.FillAscending(3, 2, values), ConfigurationError);
    EXPECT_THROW(sample.FillAscending(0, 33, values), ConfigurationError);
    EXPECT_THROW(sample.FillDescending(-1, 2, values), ConfigurationError);
    EXPECT_THROW(sample.FillDescending(3, 2, values), ConfigurationError);
    EXPECT_THROW(sample.FillDescending(0, 33, values), ConfigurationError);
}

} // namespace
} // namespace kagen
