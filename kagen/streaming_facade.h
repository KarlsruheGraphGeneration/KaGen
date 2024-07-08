#pragma once

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/kagen.h"

#include <mpi.h>

#include <string>

namespace kagen {
class StreamingGenerator {
public:
    StreamingGenerator(const std::string& options, PEID chunks_per_pe, MPI_Comm comm);

    void Initialize();

    [[nodiscard]] StreamedGraph Next();
    [[nodiscard]] bool          Continue() const;

private:
    std::unique_ptr<Generator> CreateGenerator(PEID chunk);

    void     ExchangeNonlocalEdges(const std::vector<SInt>& vertex_distribution);
    Edgelist GenerateBadHyperbolicEdges(PEID chunk);

    PGeneratorConfig config_;

    PEID next_streaming_chunk_;
    PEID streaming_chunks_per_pe_;

    PEID     rank_;
    PEID     size_;
    MPI_Comm comm_;

    std::unique_ptr<GeneratorFactory> factory_;

    std::vector<VertexRange> my_vertex_ranges_;
    std::vector<Edgelist>    nonlocal_edges_;
};
} // namespace kagen
