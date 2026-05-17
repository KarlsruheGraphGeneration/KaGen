#include "kagen/io/hmetis_hypergraph.h"

#include "kagen/io/buffered_writer.h"

#include <mpi.h>

#include <algorithm>
#include <stdexcept>

namespace kagen {
namespace {
SInt NumHyperedgesFromOffsets(const std::vector<SInt>& offsets) {
    return offsets.empty() ? 0 : static_cast<SInt>(offsets.size() - 1);
}
} // namespace

void HmetisHypergraphWriter::WriteHeader(const std::string& filename, const SInt n, const SInt /*m*/) {
    IgnoresEdgeWeights();

    const SInt num_pin_hyperedges   = NumHyperedgesFromOffsets(graph_.hyperedge_offsets);
    const SInt num_range_hyperedges = NumHyperedgesFromOffsets(graph_.hyperedge_range_offsets);

    if (num_pin_hyperedges != 0 && num_range_hyperedges != 0 && num_pin_hyperedges != num_range_hyperedges) {
        throw IOError("Hyperedge pin offsets and range offsets describe different numbers of hyperedges");
    }

    const SInt num_hyperedges = std::max(num_pin_hyperedges, num_range_hyperedges);

    BufferedTextOutput<> out(tag::append, filename);
    out.WriteInt(num_hyperedges).WriteChar(' ').WriteInt(n).WriteChar('\n').Flush();
}

bool HmetisHypergraphWriter::WriteBody(const std::string& filename) {
    IgnoresEdgeWeights();

    BufferedTextOutput<> out(tag::append, filename);

    const auto& offsets       = graph_.hyperedge_offsets;
    const auto& pins          = graph_.hyperedge_pins;
    const auto& range_offsets = graph_.hyperedge_range_offsets;
    const auto& ranges        = graph_.hyperedge_ranges;

    const SInt num_pin_hyperedges   = NumHyperedgesFromOffsets(offsets);
    const SInt num_range_hyperedges = NumHyperedgesFromOffsets(range_offsets);

    if (num_pin_hyperedges == 0 && num_range_hyperedges == 0) {
        return info_.has_vertex_weights;
    }

    if (num_pin_hyperedges != 0 && num_range_hyperedges != 0 && num_pin_hyperedges != num_range_hyperedges) {
        throw IOError("Hyperedge pin offsets and range offsets describe different numbers of hyperedges");
    }

    const SInt num_hyperedges = std::max(num_pin_hyperedges, num_range_hyperedges);

    for (SInt e = 0; e < num_hyperedges; ++e) {
        bool first = true;

        if (num_pin_hyperedges != 0) {
            const SInt begin = offsets[e];
            const SInt end   = offsets[e + 1];

            for (SInt i = begin; i < end; ++i) {
                if (!first) {
                    out.WriteChar(' ');
                }

                out.WriteInt(pins[i] + 1); // hMetis is 1-based
                first = false;
            }
        }

        if (num_range_hyperedges != 0) {
            const SInt begin = range_offsets[e];
            const SInt end   = range_offsets[e + 1];

            for (SInt i = begin; i < end; ++i) {
                const PinRange& range = ranges[i];

                for (SInt v = range.begin; v < range.end; ++v) {
                    if (!first) {
                        out.WriteChar(' ');
                    }

                    out.WriteInt(v + 1); // hMetis is 1-based
                    first = false;
                }
            }
        }

        out.WriteChar('\n').Flush();
    }

    return info_.has_vertex_weights;
}

void HmetisHypergraphWriter::WriteFooter(const std::string& filename) {
    BufferedTextOutput<> out(tag::append, filename);
    for (const SSInt weight: graph_.vertex_weights) {
        out.WriteInt(weight).WriteChar('\n').Flush();
    }
}

} // namespace kagen