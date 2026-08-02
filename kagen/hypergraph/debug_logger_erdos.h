#pragma once

#include "kagen/hypergraph/debug_logger.h"
#include "kagen/kagen.h"

#include <cstdint>
#include <string>

namespace kagen {

struct ErdosHyperedgeDebugEvent {
    SInt hyperedge_id         = 0;
    SInt hyperedge_size       = 0;
    SInt minimum_vertex       = 0;
    SInt sampling_attempts    = 0;
    SInt duplicate_rejections = 0;
    SInt minimum_search_steps = 0;
    SInt minimum_cache_gets   = 0;

    std::int64_t duration_ns = 0;
};

class ErdosHypergraphDebugLogger final : public CsvDebugLogger {
public:
    explicit ErdosHypergraphDebugLogger(const std::string& filename, const bool append = false)
        : CsvDebugLogger(
              filename, append,
              "hyperedge_id,hyperedge_size,minimum_vertex,sampling_attempts,"
              "duplicate_rejections,minimum_search_steps,minimum_cache_gets,duration_ns") {}

    void LogHyperedge(const ErdosHyperedgeDebugEvent& event) {
        out() << event.hyperedge_id << ',' << event.hyperedge_size << ',' << event.minimum_vertex << ','
              << event.sampling_attempts << ',' << event.duplicate_rejections << ',' << event.minimum_search_steps
              << ',' << event.minimum_cache_gets << ',' << event.duration_ns << '\n';
    }
};

} // namespace kagen