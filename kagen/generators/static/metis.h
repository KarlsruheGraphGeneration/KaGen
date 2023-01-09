#pragma once

#include <iostream>

#include "kagen/generators/static/graph_reader.h"
#include "kagen/generators/static/mmap_toker.h"

namespace kagen::staticgraph {
class MetisReader : public GraphReader {
public:
    MetisReader(const std::string& filename);

    GraphSize ReadSize() final;

    Graph Read(SInt from, SInt to) final;

    SInt FindNodeByEdge(SInt edge) final;

private:
    MappedFileToker toker_;

    SInt cached_first_node_     = std::numeric_limits<SInt>::max();
    SInt cached_first_node_pos_ = std::numeric_limits<SInt>::max();
};
} // namespace kagen::staticgraph
