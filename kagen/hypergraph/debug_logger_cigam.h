#pragma once

#include "kagen/hypergraph/debug_logger.h"
#include "kagen/kagen.h"

#include <cstdint>
#include <string>

namespace kagen {

struct CIGAMHyperedgeDebugEvent {
    SInt hyperedge_id    = 0;
    SInt mode            = 0;
    SInt hyperedge_size  = 0;
    SInt dominant_vertex = 0;
    SInt layer           = 0;
    SInt endpoint_vertex = 0;
    SInt block_j_min     = 0;
    SInt block_j_max     = 0;

    long double block_log_population = 0.0L;
    long double log_probability      = 0.0L;

    SInt sampling_attempts    = 0;
    SInt duplicate_rejections = 0;

    std::int64_t duration_ns = 0;
};

class CIGAMHypergraphDebugLogger final : public CsvDebugLogger {
public:
    explicit CIGAMHypergraphDebugLogger(const std::string& filename, const bool append = false)
        : CsvDebugLogger(
              filename, append,
              "hyperedge_id,mode,hyperedge_size,dominant_vertex,layer,"
              "endpoint_vertex,block_j_min,block_j_max,block_log_population,"
              "log_probability,sampling_attempts,duplicate_rejections,duration_ns") {}

    void LogHyperedge(const CIGAMHyperedgeDebugEvent& event) {
        out() << event.hyperedge_id << ',' << event.mode << ',' << event.hyperedge_size << ',' << event.dominant_vertex
              << ',' << event.layer << ',' << event.endpoint_vertex << ',' << event.block_j_min << ','
              << event.block_j_max << ',' << event.block_log_population << ',' << event.log_probability << ','
              << event.sampling_attempts << ',' << event.duplicate_rejections << ',' << event.duration_ns << '\n';
    }
};

} // namespace kagen