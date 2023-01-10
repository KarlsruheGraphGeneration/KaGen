#pragma once

#include "kagen/generators/static/graph_reader.h"
#include "kagen/generators/static/mmap_toker.h"

namespace kagen::staticgraph {
class MetisReader : public GraphReader {
public:
    MetisReader(const std::string& filename);

    GraphSize ReadSize() final;

    Graph Read(SInt from, SInt to, SInt num_edges) final;

    SInt FindNodeByEdge(SInt edge) final;

private:
    MappedFileToker toker_;

    SInt cached_first_node_     = 0;
    SInt cached_first_edge_     = 0;
    SInt cached_first_node_pos_ = 0;
};
} // namespace kagen::staticgraph
