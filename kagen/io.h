#pragma once

#include "kagen/context.h"
#include "kagen/io/graph_format.h"
#include "kagen/kagen.h"

#include <mpi.h>

#include <string>
#include <unordered_map>

namespace kagen {
class IOError : public std::exception {
public:
    IOError(std::string what) : _what(std::move(what)) {}

    const char* what() const noexcept override {
        return _what.c_str();
    }

private:
    std::string _what;
};

const std::unordered_map<FileFormat, std::unique_ptr<FileFormatFactory>>& GetGraphFormatFactories();

const std::unique_ptr<FileFormatFactory>& GetGraphFormatFactory(FileFormat format);

std::unique_ptr<GraphReader>
CreateGraphReader(const std::string& filename, const InputGraphConfig& config, PEID rank, PEID size);

std::unique_ptr<GraphReader>
CreateGraphReader(const FileFormat format, const InputGraphConfig& config, PEID rank, PEID size);

struct GraphFragment {
    Graph graph;
    int   deficits;

    //! True if this fragment was produced by a direct strict edge-balanced read (GraphReader::ReadStrictEdgeRange)
    //! rather than the vertex-balanced read that feeds the redistribution postprocessing pass.
    bool strict_edge_direct = false;
    //! True if the direct read must be validated against global tail-sortedness before it can be trusted
    //! (edge-list readers); if the check fails, FinalizeGraphFragment falls back to the redistribution path.
    bool strict_needs_sort_check = false;
    //! True for --distribution=balance-edges (as opposed to balance-edges-strict): after the direct slice read,
    //! heal the partial boundary vertices to whole-vertex (vertex-atomic) ownership via a neighbor exchange, so
    //! no vertex ends up split. If healing is not possible (empty PEs or a vertex spanning >2 PEs), fall back.
    bool heal_to_vertex_atomic = false;
    //! Total number of vertices in the graph when known exactly (0 => derive from the edge list). Only used by
    //! the direct strict-read path.
    SInt n = 0;
};

GraphFragment ReadGraphFragment(
    GraphReader& reader, GraphRepresentation representation, const InputGraphConfig& config, PEID rank, PEID size);

Graph FinalizeGraphFragment(GraphFragment fragment, const InputGraphConfig& config, bool output, MPI_Comm comm);

void WriteGraph(GraphWriter& writer, const OutputGraphConfig& config, bool output, MPI_Comm comm);
} // namespace kagen
