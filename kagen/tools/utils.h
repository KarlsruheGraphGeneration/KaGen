#pragma once

#include "kagen/kagen.h"

#include <mpi.h>

namespace kagen {
inline std::pair<SInt, SInt> ComputeRange(const SInt n, const PEID size, const PEID rank) {
    const SInt chunk = n / size;
    const SInt rem   = n % size;
    const SInt from  = rank * chunk + std::min<SInt>(rank, rem);
    const SInt to    = std::min<SInt>(from + ((static_cast<SInt>(rank) < rem) ? chunk + 1 : chunk), n);
    return {from, to};
}

inline SInt FindNumberOfVerticesInEdgelist(const Edgelist& edges, MPI_Comm comm) {
    SInt n = 0;
    for (const auto& [u, v]: edges) {
        n = std::max(n, std::max(u, v));
    }
    MPI_Allreduce(MPI_IN_PLACE, &n, 1, KAGEN_MPI_SINT, MPI_MAX, comm);
    return n + 1;
}

inline PEID GetCommRank(MPI_Comm comm) {
    PEID rank;
    MPI_Comm_rank(comm, &rank);
    return rank;
}

inline PEID GetCommSize(MPI_Comm comm) {
    PEID size;
    MPI_Comm_size(comm, &size);
    return size;
}

template <typename T>
T FloorLog2(const T arg) {
    constexpr std::size_t arg_width = std::numeric_limits<T>::digits;
    static_assert(
        arg_width == std::numeric_limits<unsigned long>::digits
            || arg_width == std::numeric_limits<unsigned int>::digits,
        "unsupported data type width");

    T log2 = static_cast<T>(arg_width);
    if constexpr (arg_width == std::numeric_limits<unsigned int>::digits) { // 32 bit
        log2 -= __builtin_clz(arg);
    } else { // 64 bit
        log2 -= __builtin_clzl(arg);
    }

    return log2 - 1;
}
} // namespace kagen
