/*******************************************************************************
 * rmat/util.hpp
 *
 * Utilities
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef RMAT_UTIL_HEADER
#define RMAT_UTIL_HEADER

#include <tlx/logger.hpp>

#include <string>
#include <vector>

namespace rmat {

template <typename Array>
void log_arr(const bool debug, const Array &arr, size_t size,
             const std::string &desc, size_t max_log_size = 100) {
    using value_type = std::remove_reference_t<decltype(arr[0])>;
    if (size > max_log_size) return;
    LOG << desc << ": "
        << std::vector<value_type>(arr, arr + std::min(max_log_size, size))
        << (max_log_size < size ? " (truncated)" : "");
}

#define LOG_ARR(arr, desc) log_arr(debug, (arr), size_, (desc))
#define LOG_ARRs(arr, desc, size) log_arr(debug, (arr), (size), (desc))

}  // namespace rmat

#endif // RMAT_UTILI_HEADER
