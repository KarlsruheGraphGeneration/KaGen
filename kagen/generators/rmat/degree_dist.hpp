/*******************************************************************************
 * rmat/degree_dist.hpp
 *
 * Graph utils for R-MAT graph generator
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef RMAT_DEGREE_DIST_HEADER
#define RMAT_DEGREE_DIST_HEADER

#include <tlx/logger.hpp>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <unordered_map>
#include <vector>

// Degree distribution of a (generated) graph
// out_deg: whether to use out-degree or in-degree of vertices
template <typename node_t, bool out_deg = false>
class degree_dist {
    static constexpr bool debug = false;
public:
    degree_dist() = default;

    void add_edge(const node_t &src, const node_t &dst) {
        if constexpr (out_deg) {
            degrees[src]++;
        } else {
            degrees[dst]++;
        }
    }

    degree_dist operator + (const degree_dist &other) const {
        degree_dist dist;
        dist.degrees = degrees;
        dist += other;
    }

    degree_dist& operator += (const degree_dist &other) {
        for (const auto &pair : other.degrees) {
            degrees[pair.first] += pair.second;
        }
        return *this;
    }

    uint32_t score(const node_t &node) const {
        return degrees[node];
    }

    // return mapping degree -> frequency
    std::unordered_map<uint32_t, uint32_t> histogram() const {
        std::unordered_map<uint32_t, uint32_t> degdist;
        // for each node, increase its degrees counter
        for (const auto& pair : degrees) {
            degdist[pair.second]++;
        }
        return degdist;
    }

    // return vector of frequencies
    std::vector<uint32_t> histogram_vec() const {
        auto deg_dist = histogram(); // degree -> #nodes
        uint32_t max_deg = std::max_element(
            deg_dist.begin(), deg_dist.end(),
            [](const auto &p1, const auto &p2) {
                return p1.first < p2.first;
            })->first;
        LOG << "max degree: " << max_deg;
        std::vector<uint32_t> hist_vec(max_deg + 1, 0);
        for (const auto &pair : deg_dist) {
            hist_vec[pair.first] = pair.second;
        }
        return hist_vec;
    }

    void write_histogram(const std::string &outfilename) const {
        std::ofstream out(outfilename);
        assert(out.is_open());
        auto hist_vec = histogram_vec();
        for (auto it = hist_vec.begin(); it != hist_vec.end(); ++it) {
            if (*it == 0) continue;
            out << it - hist_vec.begin() << " " << *it << std::endl;
        }
    }

protected:
    std::unordered_map<node_t, uint32_t> degrees;
};

#endif // RMAT_DEGREE_DIST_HEADER
