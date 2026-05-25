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

void HmetisHypergraphWriter::WriteHeader(const std::string& filename, const SInt n, const SInt m) {
    IgnoresEdgeWeights();

    const SInt num_pin_hyperedges   = NumHyperedgesFromOffsets(graph_.hyperedge_offsets);
    const SInt num_range_hyperedges = NumHyperedgesFromOffsets(graph_.hyperedge_range_offsets);

    if (num_pin_hyperedges != 0 && num_range_hyperedges != 0 && num_pin_hyperedges != num_range_hyperedges) {
        throw IOError("Hyperedge pin offsets and range offsets describe different numbers of hyperedges");
    }

    BufferedTextOutput<> out(tag::append, filename);
    out.WriteInt(m).WriteChar(' ').WriteInt(n).WriteChar('\n').Flush();
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

    auto WritePin = [&](const SInt v, bool& first) {
        if (!first) {
            out.WriteChar(' ');
        }

        out.WriteInt(v + 1); // hMetis is 1-based
        first = false;
    };

    for (SInt e = 0; e < num_hyperedges; ++e) {
        bool first = true;

        SInt       p     = num_pin_hyperedges != 0 ? offsets[e] : 0;
        const SInt p_end = num_pin_hyperedges != 0 ? offsets[e + 1] : 0;

        SInt       r     = num_range_hyperedges != 0 ? range_offsets[e] : 0;
        const SInt r_end = num_range_hyperedges != 0 ? range_offsets[e + 1] : 0;

        SInt last_written = -1;

        auto WriteUniquePin = [&](const SInt v) {
            if (v != last_written) {
                WritePin(v, first);
                last_written = v;
            }
        };

        while (p < p_end || r < r_end) {
            if (r >= r_end) {
                WriteUniquePin(pins[p++]);
            } else if (p >= p_end) {
                const PinRange& range = ranges[r++];

                for (SInt v = range.begin; v < range.end; ++v) {
                    WriteUniquePin(v);
                }
            } else {
                const SInt      pin   = pins[p];
                const PinRange& range = ranges[r];

                if (pin < range.begin) {
                    WriteUniquePin(pin);
                    ++p;
                } else if (pin < range.end) {
                    // Pin is already covered by the range.
                    ++p;
                } else {
                    for (SInt v = range.begin; v < range.end; ++v) {
                        WriteUniquePin(v);
                    }
                    ++r;
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