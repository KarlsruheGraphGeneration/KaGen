#pragma once

#include "kagen/definitions.h"

namespace kagen {
void SortEdges(EdgeList& edge_list);

void AddReverseEdges(EdgeList& edge_list, VertexRange vertex_range);

void AddReverseEdgesAndRedistribute(EdgeList& edge_list, VertexRange vertex_range);
} // namespace kagen
