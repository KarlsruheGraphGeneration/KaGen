#pragma once

#include <fstream>

#include "kagen/generators/static/graph_reader.h"

namespace kagen::staticgraph {
class BinaryParhipReader : public GraphReader {
public:
    BinaryParhipReader(const std::string& filename);

    GraphSize ReadSize() final;

    Graph Read(SInt from, SInt to_node, SInt to_edge) final;

    SInt FindNodeByEdge(SInt edge) final;

private:
    std::ifstream in_;

    SInt n_       = 0;
    SInt m_       = 0;
    SInt version_ = 0;
};
} // namespace kagen::staticgraph
