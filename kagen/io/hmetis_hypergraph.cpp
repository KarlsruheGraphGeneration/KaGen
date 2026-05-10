#include "kagen/io/hmetis_hypergraph.h"

#include "kagen/io/buffered_writer.h"

#include <mpi.h>

namespace kagen {
void HmetisHypergraphWriter::WriteHeader(const std::string& filename, const SInt n, const SInt m) {
    IgnoresEdgeWeights();

    BufferedTextOutput<> out(tag::append, filename);
    out.WriteInt(m).WriteChar(' ').WriteInt(n).WriteChar('\n').Flush();
};

bool HmetisHypergraphWriter::WriteBody(const std::string& filename) {
    IgnoresEdgeWeights();

    BufferedTextOutput<> out(tag::append, filename);

    const auto& offsets = graph_.hyperedge_offsets;
    const auto& pins    = graph_.hyperedge_pins;

    if (offsets.empty()) {
        return info_.has_vertex_weights;
    }

    for (SInt e = 0; e + 1 < static_cast<SInt>(offsets.size()); ++e) {
        const SInt begin = offsets[e];
        const SInt end   = offsets[e + 1];

        for (SInt i = begin; i < end; ++i) {
            if (i != begin) {
                out.WriteChar(' ');
            }
            out.WriteInt(pins[i] + 1); // hMetis is 1-based
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
