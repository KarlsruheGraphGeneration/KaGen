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
};

GraphFragment ReadGraphFragment(
    GraphReader& reader, GraphRepresentation representation, const InputGraphConfig& config, PEID rank, PEID size);

Graph FinalizeGraphFragment(GraphFragment fragment, bool output, MPI_Comm comm);

void WriteGraph(GraphWriter& writer, const OutputGraphConfig& config, bool output, MPI_Comm comm);
} // namespace kagen
