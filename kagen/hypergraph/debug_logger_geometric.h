#pragma once

#include "kagen/hypergraph/debug_logger.h"
#include "kagen/kagen.h"

namespace kagen {

class GeometricHypergraphDebugLogger final : public CsvDebugLogger {
public:
    explicit GeometricHypergraphDebugLogger(const std::string& filename, const bool append = false)
        : CsvDebugLogger(
              filename, append,
              "hyperedge_id,hyperedge_center,radius,candidate_cells,inside_cells,partial_cells,outside_cells,"
              "emitted_pins,emitted_ranges,estimated_size,duration_ns,inside_estimated_size,"
              "partial_estimated_size") {}

    void LogHyperedge(
        const SInt hyperedge_id, const std::string& hyperedge_center, const double radius, const SInt candidate_cells,
        const SInt inside_cells, const SInt partial_cells, const SInt outside_cells, const SInt emitted_pins,
        const SInt emitted_ranges, const SInt estimated_size, const long long duration_ns,
        const SInt inside_estimated_size, const SInt partial_estimated_size) {
        out() << hyperedge_id << ',' << hyperedge_center << ',' << radius << ',' << candidate_cells << ','
              << inside_cells << ',' << partial_cells << ',' << outside_cells << ',' << emitted_pins << ','
              << emitted_ranges << ',' << estimated_size << ',' << duration_ns << ',' << inside_estimated_size << ','
              << partial_estimated_size << '\n';
    }
};

} // namespace kagen