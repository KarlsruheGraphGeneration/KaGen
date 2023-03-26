/*******************************************************************************
 * rmat/rmat.hpp
 *
 * R-MAT graph generator using weighted random sampling
 *
 * The method implemented in this file is described in the following paper:
 * @article{HubSan2019RMAT,
 *   title={Linear Work Generation of {R-MAT} Graphs},
 *   author={H{\"u}bschle-Schneider, Lorenz and Sanders, Peter},
 *   journal={arXiv preprint arXiv:1905.03525},
 *   year={2019}
 * }
 *
 * Copyright (C) 2018-2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/
#pragma once

#ifndef RMAT_RMAT_HEADER
    #define RMAT_RMAT_HEADER

    #include "kagen/generators/rmat/alias_key.hpp"
    #include "kagen/generators/rmat/graph500.hpp"
    #include "kagen/tlx/aggregate.hpp"
    #include "kagen/tlx/clz.hpp"

    #include <array>
    #include <ios>
    #include <queue>
    #include <utility>
    #include <vector>

namespace rmat {
template <bool enabled>
struct rmat_stats_t;

template <>
struct rmat_stats_t<false> {
    void record_sample_drawn() {}
    void record_edge_done() {}
    void record_sample_bits(unsigned) {}
    // TODO get rid of these?
    tlx::Aggregate<double> samples_per_edge;
    tlx::Aggregate<double> bits_per_sample;
};

template <>
struct rmat_stats_t<true> {
    void record_sample_drawn() {
        lvl++;
    }
    void record_edge_done() {
        samples_per_edge.add(lvl);
        lvl = 0;
    }
    void record_sample_bits(unsigned bits) {
        bits_per_sample.add(bits);
    }

    unsigned               lvl = 0;
    tlx::Aggregate<double> samples_per_edge;
    tlx::Aggregate<double> bits_per_sample;
};

template <bool Scramble_IDs = false>
class rmat {
public:
    using node   = int64_t;
    using prefix = uint32_t;
    using path   = uint64_t;
    using edge   = std::pair<node, node>;
    using entry  = std::pair<prefix, double>;

    static constexpr unsigned prefix_half_bits = 4 * sizeof(prefix);
    static constexpr unsigned path_half_bits   = 4 * sizeof(path);
    static constexpr prefix   dst_mask         = (prefix{1} << prefix_half_bits) - 1;

    static constexpr bool        debug         = false;
    static constexpr bool        collect_stats = false;
    static constexpr bool        scramble_ids  = Scramble_IDs;
    static constexpr const char* name          = "fixdepth";

    template <typename RNG>
    rmat(RNG& rng, int log_n_, double a_, double b_, double c_)
        : quad({a_, b_, c_, 1.0 - (a_ + b_ + c_)}),
          log_n(log_n_),
          node_mask((path{1} << log_n_) - 1),
          scramble_state(graph500::init_scramble_state(rng, log_n)) {
        // sLOG << "RMAT: need" << 2 * log_n << "path bits";
        // sLOG << "RMAT: node extraction mask" << std::hex << node_mask;
    }

    void init(int max_depth) {
        timer t;

        table_depth         = max_depth;
        marker_pos          = prefix_half_bits + table_depth;
        marker_removal_mask = ~(prefix{1} << marker_pos);

        entries.clear();

        enum_items(0, 1.0, 0);

        // sLOG1 << "generated" << entries.size() << "path prefixes in" << t.get_and_reset() << "ms";

        table.init(entries.size());
        // sLOG1 << "init table in" << t.get_and_reset() << "ms";

        table.construct(entries.begin(), entries.end(), /* is_dist */ true);
        // sLOG1 << "construct table in" << t.get_and_reset() << "ms";
    }

    template <typename RNG>
    std::pair<node, node> get_edge(size_t /* ei = ignored */, RNG& rng) const {
        return get_edge(rng);
    }

    template <typename RNG>
    std::pair<node, node> get_edge(RNG& rng) const {
        path result        = 0;
        int  bits_per_half = 0;

        do {
            prefix next = table.sample(rng.next());
            stats.record_sample_drawn();

            path     src_part, dst_part;
            unsigned bits_in_prefix = split_prefix(next, src_part, dst_part);

            /*
            sLOG << "Prefix" << std::hex << next << "with" << std::dec << bits_in_prefix
                 << "bits per half; src half:" << std::hex << src_part << "dst:" << dst_part
                 << "old tentative result:" << result;
            */

            bits_per_half += bits_in_prefix;
            result <<= bits_in_prefix;
            result |= (src_part << path_half_bits);
            result |= dst_part;
        } while (bits_per_half < log_n);
        // remove unneeded bits
        unsigned shift = bits_per_half - log_n;
        /*
        sLOG << "got enough bits in" << std::hex << result << "-- removing" << std::dec << shift << "of"
             << bits_per_half;
        */
        result >>= shift;
        node src = static_cast<node>(result >> path_half_bits);
        node dst = static_cast<node>(result & node_mask);
        // sLOG << "Extracted nodes" << std::hex << src << dst << "from" << result;
        stats.record_edge_done();

        if constexpr (scramble_ids) {
            // scrambles in-place
            graph500::scramble_two(src, dst, scramble_state);
        }

        return std::make_pair(src, dst);
    }

    template <typename RNG, typename Callback>
    void get_edges(Callback&& callback, size_t min, size_t max, RNG& rng) const {
        get_edges(callback, max - min, rng);
    }

    template <typename RNG, typename Callback>
    void get_edges(Callback&& callback, size_t num_edges, RNG& rng) const {
        path result        = 0;
        int  bits_per_half = 0;

        for (size_t i = 0; i < num_edges; /* nothing! */) {
            prefix next = table.sample(rng.next());
            stats.record_sample_drawn();

            path     src_part, dst_part;
            unsigned bits_in_prefix = split_prefix(next, src_part, dst_part);

            /*
            sLOG << "Prefix" << std::hex << next << "with" << std::dec << bits_in_prefix
                 << "bits per half; src half:" << std::hex << src_part << "dst:" << dst_part
                 << "old tentative result:" << result;
            */

            bits_per_half += bits_in_prefix;
            result <<= bits_in_prefix;
            result |= (src_part << path_half_bits);
            result |= dst_part;

            if (bits_per_half >= log_n) {
                /*
                sLOG << "got enough bits in" << std::hex << result << "-- have" << std::dec << bits_per_half << "need"
                     << log_n;
                */
                // we have enough bits for an edge, extract it
                int  excess = bits_per_half - log_n;
                path tmp = (result >> excess), removal_mask = (path{1} << excess) - 1;
                node src = (tmp >> path_half_bits), dst = tmp & node_mask;

                // sLOG << "Extracted nodes" << std::hex << src << dst << "from" << result;
                //  implement clip-and-flip
    #ifdef RMAT_CLIPFLIP
                if (src > dst)
                    std::swap(src, dst);
    #endif // RMAT_CLIPFLIP
           // emit the edge
                if constexpr (scramble_ids) {
                    // scrambles in-place
                    graph500::scramble_two(src, dst, scramble_state);
                }
                stats.record_edge_done();
                callback(src, dst);
                // this is placed here because compilers' optimisers are daft
                removal_mask |= (removal_mask << path_half_bits);
                i++;
                // Reuse the remaining `leftover` bits
                bits_per_half -= log_n;
                result &= removal_mask;
                // sLOG << bits_per_half << "bits remaining, tentative result:" << std::hex << result;
            }
        }
    }

    size_t table_size() const {
        return table.size();
    }

    tlx::Aggregate<double> get_depth_stats() const {
        return stats.samples_per_edge;
    }

    tlx::Aggregate<double> get_sample_stats() const {
        return stats.bits_per_sample;
    }

protected:
    template <typename integral>
    constexpr unsigned msb(integral i) const {
        return 8 * sizeof(integral) - tlx::clz(i) - 1;
    }

    TLX_ATTRIBUTE_ALWAYS_INLINE
    constexpr unsigned split_prefix(const prefix& in, path& out_src, path& out_dst) const {
        assert(in < (path{1} << (table_depth + prefix_half_bits)));
        // We know how many bits are in each prefix, no need to calculate it
        out_dst = (in & dst_mask);
        assert(out_dst < (path{1} << table_depth));
        out_src = (in >> prefix_half_bits);
        stats.record_sample_bits(table_depth);
        return table_depth;
    }

    void enum_items(prefix pref, double prob, int depth) {
        for (size_t d = 0; d < quad.size(); ++d) {
            prefix   new_prefix = pref;
            unsigned dst_pos    = depth;
            unsigned src_pos    = depth + prefix_half_bits;

            /*
            sLOG << "Considering bits" << src_pos << "and" << dst_pos << "of" << pref
                 << "src half:" << (new_prefix >> prefix_half_bits) << "dst half:" << (new_prefix & dst_mask)
                 << "direction" << d;
            */

            // Then add this position and the new marker bit
            new_prefix |= ((d >> 1) << src_pos);
            new_prefix |= ((d & 1) << dst_pos);

            /*
            sLOG << "\tresult:" << new_prefix << "src half:" << (new_prefix >> prefix_half_bits)
                 << "dst half:" << (new_prefix & dst_mask);
            */

            const double new_prob = prob * quad[d];
            if (depth + 1 < table_depth) {
                enum_items(new_prefix, new_prob, depth + 1);
            } else {
                entries.emplace_back(new_prefix, new_prob);
            }
        }
    }

    alias_key<prefix>     table;
    std::array<double, 4> quad;
    std::vector<entry>    entries;

    const int  log_n;
    const path node_mask;
    unsigned   marker_pos;
    prefix     marker_removal_mask;
    int        table_depth;

    graph500::scramble_state scramble_state;

    mutable rmat_stats_t<collect_stats> stats;
};
} // namespace rmat
#endif // RMAT_RMAT_HEADER
