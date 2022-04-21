#pragma once

#include "kagen/definitions.h"

namespace kagen {
bool ValidateVertexRanges(const EdgeList& edge_list, VertexRange vertex_range, bool expect_consecutive = true);

bool ValidateSimpleGraph(EdgeList& edge_list, VertexRange vertex_range);
} // namespace kagen
