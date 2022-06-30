/*******************************************************************************
 * rmat/parallel_do.cpp
 *
 * Pthreads parallelization helper
 *
 * Copyright (C) 2018-2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#include "kagen/generators/rmat/parallel_do.hpp"

namespace rmat {

int              g_num_numa_nodes;
int              g_total_threads;
tlx::ThreadPool* global_pool;

void init_threads(int threads) {
    g_total_threads  = threads;
    g_num_numa_nodes = 1;
    global_pool      = new tlx::ThreadPool(threads);
}

void release_threads() {
    global_pool->loop_until_empty();
    global_pool->terminate();
    delete global_pool;
    global_pool      = nullptr;
    g_num_numa_nodes = 0;
    g_total_threads  = 0;
}

int get_num_threads() {
    return g_total_threads;
}

tlx::ThreadPool* get_pool() {
    return global_pool;
}

} // namespace rmat
