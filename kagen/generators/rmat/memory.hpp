/*******************************************************************************
 * rmat/memory.hpp
 *
 * Memory helpers
 *
 * Copyright (C) 2018-2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef RMAT_MEMORY_HEADER
#define RMAT_MEMORY_HEADER

#include <tlx/logger.hpp>

#include <sys/mman.h> // madvise
#include <cstdlib>
#include <cstring> // memcpy
#include <memory>

namespace rmat {

constexpr size_t align_size(size_t size, size_t alignment) {
    return ((size + alignment - 1) / alignment) * alignment;
}

// Allocate memory using huge pages
void* alloc_hugepage(size_t size);

// Allocate memory, using huge pages for allocations larger than 1MB
void* allocate(size_t size);

struct deallocator {
    template <typename T>
    void operator()(T *ptr) {
        free((void*)ptr);
    }
};

// to avoid alloc-dealloc-mismatches
template <typename T>
using alloc_arr_ptr = std::unique_ptr<T[], deallocator>;

template <typename T>
alloc_arr_ptr<T> make_alloc_arr(size_t num_elems) {
    T* ptr = static_cast<T*>(allocate(num_elems * sizeof(T)));
    return alloc_arr_ptr<T>(ptr);
}

template <typename T>
alloc_arr_ptr<T> copy_alloc_arr(T* other, size_t num_elems) {
    T* ptr = static_cast<T*>(allocate(num_elems * sizeof(T)));
    memcpy(ptr, other, num_elems * sizeof(T));
    return alloc_arr_ptr<T>(ptr);
}

} // namespace rmat

#endif // RMAT_MEMORY_HEADER
