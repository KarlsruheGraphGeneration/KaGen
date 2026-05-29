#pragma once

#include "kagen/kagen.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>

namespace kagen {

class HypergraphDebugLogger {
public:
    explicit HypergraphDebugLogger(const std::string& filename, const bool append = false)
        : out_(filename, append ? (std::ios::out | std::ios::app) : (std::ios::out | std::ios::trunc)) {
        const bool should_write_header =
            !append || !std::filesystem::exists(filename) || std::filesystem::file_size(filename) == 0;
        if (should_write_header) {
            out_ << "hyperedge_id,hyperedge_center,radius,candidate_cells,inside_cells,partial_cells,outside_cells,"
                    "emitted_pins,emitted_ranges,estimated_size,duration_ns\n";
        }
    }

    void LogHyperedge(
        const SInt hyperedge_id, const std::string hyperedge_center, const double radius, const SInt candidate_cells,
        const SInt inside_cells, const SInt partial_cells, const SInt outside_cells, const SInt emitted_pins,
        const SInt emitted_ranges, const SInt estimated_size, const long long duration_ns) {
        out_ << hyperedge_id << /*',' << hyperedge_center <<*/ ',' << radius << ',' << candidate_cells << ','
             << inside_cells << ',' << partial_cells << ',' << outside_cells << ',' << emitted_pins << ','
             << emitted_ranges << ',' << estimated_size << ',' << duration_ns << '\n';
    }

private:
    std::ofstream out_;
};

} // namespace kagen