/*******************************************************************************
 * rmat/parallel_do.hpp
 *
 * Pthreads parallelization helper
 *
 * Copyright (C) 2018-2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef RMAT_PARALLEL_DO_HEADER
#define RMAT_PARALLEL_DO_HEADER

#include <tlx/thread_pool.hpp>

#include <numeric>
#include <thread>
#include <vector>

namespace rmat {

void init_threads(int threads);

void release_threads();

int get_num_threads();

tlx::ThreadPool* get_pool();

template <typename F, typename size_type>
void parallel_do_range(F&& callback, size_type max_, size_type min_ = 0) {
    static_assert(std::is_integral_v<size_type>, "size_type must be integral");
    int num_threads = get_num_threads();
    size_type elems_per_thread = (max_ - min_ + num_threads - 1) / num_threads;
    size_type min = min_, max = min_ + elems_per_thread;

    for (int i = 0; i < num_threads; i++) {
        get_pool()->enqueue(
            [callback, min, max, i]() { callback(min, max, i); });
        min = max;
        max = std::min(max_, max + elems_per_thread);
        // don't spawn empty threads
        if (min == max_) break;
    }
    get_pool()->loop_until_empty();
}


} // namespace rmat

#endif // RMAT_PARALLEL_DO_HEADER
