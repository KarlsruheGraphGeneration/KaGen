#pragma once

#include "kagen/comm/comm.h"
#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/kagen.h"

#include <string>

namespace kagen {
class StreamingGenerator {
public:
    StreamingGenerator(const std::string& options, PEID chunks_per_pe, Comm& comm);

    [[nodiscard]] VertexRange EstimateVertexRange(PEID pe = -1) const;

    void Initialize();

    [[nodiscard]] StreamedGraph Next();
    [[nodiscard]] bool          Continue() const;

private:
    std::unique_ptr<Generator> CreateGenerator(PEID chunk);

    void     ExchangeNonlocalEdges();
    Edgelist GenerateBadHyperbolicEdges(PEID chunk);

    PGeneratorConfig config_;

    PEID next_streaming_chunk_;
    PEID streaming_chunks_per_pe_;

    PEID  rank_;
    PEID  size_;
    Comm& comm_;

    bool initialized_ = false;

    std::unique_ptr<GeneratorFactory> factory_;

    std::vector<SInt>        vertex_distribution_;
    std::vector<VertexRange> my_vertex_ranges_;
    std::vector<Edgelist>    nonlocal_edges_;
};

} // namespace kagen
