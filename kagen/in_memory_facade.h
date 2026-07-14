#pragma once

#include "kagen/context.h"
#include "kagen/kagen.h"

#include <mpi.h>

namespace kagen {
void  GenerateInMemoryToDisk(PGeneratorConfig config, MPI_Comm comm);
Graph GenerateInMemory(const PGeneratorConfig& config, GraphRepresentation representation, MPI_Comm comm);

/*!
 * Throws ConfigurationError if a graph with split vertices (any_split, already collectively agreed on across
 * all PEs by the caller) is combined with something that requires single-PE vertex ownership: an
 * adjacency-grouped output format, --validate-simple-graph, a cross-PE-lookup edge-weight generator, or
 * statistics. Split into its own function (rather than inlined in GenerateInMemory) so the decision logic can be
 * unit-tested directly without driving the full MPI generation pipeline.
 *
 * @param any_split Whether any vertex anywhere in the graph was actually split, already reduced across all PEs.
 * @param config The (already-normalized) generator configuration.
 */
void CheckSplitVertexCompatibility(bool any_split, const PGeneratorConfig& config);
} // namespace kagen
