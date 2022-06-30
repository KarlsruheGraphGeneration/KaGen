/*******************************************************************************
 * rmat/graph_generator.hpp
 *
 * R-MAT graph generator using weighted random sampling
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef RMAT_GRAPH_GENERATOR_HEADER
#define RMAT_GRAPH_GENERATOR_HEADER

#include <rmat/degree_dist.hpp>
#include <rmat/parallel_do.hpp>
#include <rmat/timer.hpp>

#include <tlx/logger.hpp>

#include <atomic>
#include <mutex>

// namespace rmat {

template <typename rmat_t, typename RNG, bool stats = false>
class graph_generator {
    static constexpr bool debug = true;
    static constexpr bool verbose = false;
public:
    using node_t = typename rmat_t::node;
    graph_generator(const rmat_t &r, size_t block_size)
        : r_(r)
        , block_size_(block_size) {}

    void get_edges_static(size_t num_edges, size_t seed, const bool debug = true) const {
        tlx::Aggregate<double> thread_stats;
        std::mutex stats_lock;
        auto worker = [&](size_t min, size_t max, int thread) {
            double duration = process_range(seed, min, max, thread);
            LOG << "RESULT type=staticthread"
                << " total=" << duration
                << " thread=" << thread
                << " edges=" << num_edges
                << " scramble=" << rmat_t::scramble_ids
                << " method=" << rmat_t::name
                << " paths=" << r_.table_size()
                << " threads=" << rmat::get_num_threads();
            stats_lock.lock();
            thread_stats.add(duration);
            stats_lock.unlock();
        };
        rmat::parallel_do_range(worker, num_edges);
        sLOG << "Edge generation with static load balancing stats:"
             << thread_stats;
    }

    void get_edges(size_t num_edges, size_t seed, const bool debug = true) const {
        const size_t num_blocks = (num_edges + block_size_ - 1) / block_size_;
        std::atomic<size_t> next_block(0);
        tlx::Aggregate<double> thread_stats;
        std::mutex stats_lock;
        auto worker = [&](size_t, size_t, int thread)
        {
            rmat::timer t;
            tlx::Aggregate<double> block_stats;
            degree_dist<node_t> thread_deg_stats;
            size_t my_next_block = 0, my_blocks = 0;

            while ((my_next_block = next_block.fetch_add(1)) < num_blocks) {
                sLOGC(debug && verbose)
                    << "Thread" << thread << "handling block"
                    << my_next_block << "of" << num_blocks;
                ++my_blocks;
                size_t min = my_next_block * block_size_,
                    max = std::min(num_edges, min + block_size_);
                double duration = process_range(seed, min, max, my_next_block,
                                                thread_deg_stats);
                block_stats.add(duration);
            }
            double total = t.get();
            LOG << "RESULT type=dynthread"
                << " total=" << total
                << " blocks=" << my_blocks
                << " blockavg=" << block_stats.avg()
                << " blockdev=" << block_stats.stdev()
                << " thread=" << thread
                << " edges=" << num_edges
                << " scramble=" << rmat_t::scramble_ids
                << " method=" << rmat_t::name
                << " paths=" << r_.table_size()
                << " threads=" << rmat::get_num_threads();
            stats_lock.lock();
            thread_stats.add(total);
            deg_stats += thread_deg_stats;
            stats_lock.unlock();
        };
        rmat::parallel_do_range(worker, num_edges);
        sLOG << "Edge generation with dynamic load balancing stats:"
             << thread_stats;
    }

    degree_dist<node_t>& get_stats() {
        return deg_stats;
    }

protected:
    double process_range(size_t seed, size_t min, size_t max, int block,
                         degree_dist<node_t> &thread_deg_stats) const {
        rmat::timer t;
        RNG gen(seed + block);
        std::vector<std::pair<node_t, node_t>> dummy(1);
        auto cb = [&](const node_t &src, const node_t &dst) {
            LOG0 << "got edge (" << std::hex << src << ", " << dst << ")";
            if constexpr (stats) {
                thread_deg_stats.add_edge(src, dst);
            } else {
                dummy[0] = std::make_pair(src, dst);
            }
        };
        r_.get_edges(cb, min, max, gen);
        double duration = t.get();
        sLOGC(verbose) << "Block" << block << "time" << duration;
        return duration;
    }

    const rmat_t &r_;
    mutable degree_dist<node_t> deg_stats;
    size_t block_size_;
};

#endif // RMAT_GRAPH_GENERATOR_HEADER
