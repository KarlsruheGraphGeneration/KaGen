#include <gtest/gtest.h>

#include <cmath>

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/generators/geometric/geometric_2d.h"
#include "kagen/generators/geometric/rgg.h"

using namespace kagen;

bool sortByFirstElement(const std::pair<SInt, SInt>& lhs, const std::pair<SInt, SInt>& rhs) {
    if (lhs.first == rhs.first) {
        return lhs.second < rhs.second;
    }
    return lhs.first < rhs.first;
}

TEST(KARGB, generates_graph_on_one_PE) {
    PGeneratorConfig config;
    config.n = 16;
    config.r = 0.5;
    config.coordinates = true;

    RGG2DFactory factory;
    config = factory.NormalizeParameters(config, 0, 1, false);
    auto generator = factory.Create(config, 0, 1);
    generator->Generate(GraphRepresentation::EDGE_LIST);
    const auto result = generator->Take();

    // Creating the correct edge list as a test instance
    std::vector<std::pair<SInt, SInt>> edgeList;
    for (SInt i = 0; i < config.n; i++) {
        for (SInt j = 0; j < config.n; j++) {
            auto [x1, y1] = result.coordinates.first[i];
            auto [x2, y2] = result.coordinates.first[j];
            // Comparing all coordinates
            if(i != j && std::hypot( x1- x2, y1 - y2) <= config.r) {
                edgeList.push_back(std::pair(i, j));
            }
        }
    }
    // Sorting both lists before comparing them
    std::vector<std::pair<SInt, SInt>> resultList = result.edges;
    std::sort(resultList.begin(), resultList.end(), sortByFirstElement);
    std::sort(edgeList.begin(), edgeList.end(), sortByFirstElement);
    ASSERT_EQ(resultList, edgeList);
}
