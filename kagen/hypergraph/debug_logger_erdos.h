#pragma once

#include "kagen/hypergraph/debug_logger.h"
#include "kagen/kagen.h"

#include <string>

namespace kagen {

class ErdosHypergraphDebugLogger final : public CsvDebugLogger {
public:
    explicit ErdosHypergraphDebugLogger(const std::string& filename, const bool append = false)
        : CsvDebugLogger(
              filename, append,
              "rank,k,stage,strategy,block_id,lo,hi,universe,assigned_m,"
              "remaining_m,remaining_universe,seed,duration_ns,extra") {}

    void
    LogSize(const SInt rank, const SInt k, const SInt assigned_m, const SInt remaining_m, const std::string& strategy) {
        out() << rank << ',' << k << ',' << "size," << strategy << ',' << 0 << ',' << 0 << ',' << 0 << ',' << 0 << ','
              << assigned_m << ',' << remaining_m << ',' << 0 << ',' << 0 << ',' << 0 << ',' << '\n';
    }

    void LogBlock(
        const SInt rank, const SInt k, const std::string& stage, const std::string& strategy, const SInt block_id,
        const SInt lo, const SInt hi, const std::string& universe, const SInt assigned_m, const SInt remaining_m,
        const std::string& remaining_universe, const SInt seed, const long long duration_ns,
        const std::string& extra = "") {
        out() << rank << ',' << k << ',' << stage << ',' << strategy << ',' << block_id << ',' << lo << ',' << hi << ','
              << universe << ',' << assigned_m << ',' << remaining_m << ',' << remaining_universe << ',' << seed << ','
              << duration_ns << ',' << extra << '\n';
    }

    void LogFallback(
        const SInt rank, const SInt k, const std::string& strategy, const SInt lo, const SInt hi,
        const SInt requested_m, const SInt emitted, const SInt attempts) {
        out() << rank << ',' << k << ',' << "fallback," << strategy << ',' << 0 << ',' << lo << ',' << hi << ',' << 0
              << ',' << requested_m << ',' << 0 << ',' << 0 << ',' << 0 << ',' << 0 << ',' << "emitted=" << emitted
              << ";attempts=" << attempts << '\n';
    }
};

} // namespace kagen