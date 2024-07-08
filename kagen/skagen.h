#pragma once

#include "kagen/context.h"
#include "kagen/generators/generator.h"
#include "kagen/kagen.h"

#include <mpi.h>

#include <string>

namespace kagen {
struct EdgelistChunk {
    VertexRange vertex_range;
    Edgelist    primary_edges;
    Edgelist    secondary_edges;

    template <typename EdgeConsumer>
    void ForEachEdge(EdgeConsumer&& consumer) {
        std::size_t primary_idx   = 0;
        std::size_t secondary_idx = 0;

        for (SInt u = vertex_range.first; u < vertex_range.second; ++u) {
            while (primary_idx < primary_edges.size() && primary_edges[primary_idx].first == u) {
                consumer(primary_edges[primary_idx].first, primary_edges[primary_idx].second);
                ++primary_idx;
            }

            while (secondary_idx < secondary_edges.size() && secondary_edges[secondary_idx].first == u) {
                consumer(secondary_edges[secondary_idx].first, secondary_edges[secondary_idx].second);
                ++secondary_idx;
            }
        }
    }
};

class sKaGen {
public:
    /**
     * @param options The options string to be passed to KaGen, e.g, `rhg;N=10;M=12`.
     * @param chunks Lower bound on the number of chunks *per PE* that generation will be split into. Larger chunk
     * counts are possible at the discretion of the generator.
     * @param comm The MPI communicator to be used.
     */
    sKaGen(const std::string& options, PEID chunks_per_pe, MPI_Comm comm);

    /**
     * This function must be called before the first call to Continue().
     * Depending on the generator, this function may run for a long time.
     */
    void Initialize();

    /**
     * Continue generating the next chunk of the graph.
     *
     * @param chunk The chunk to be filled with the next set of edges.
     * @return True if the generation is not finished, false otherwise.
     */
    [[nodiscard]] bool Continue(EdgelistChunk& chunk);

private:
    std::unique_ptr<Generator> CreateGenerator(PEID chunk);

    void     ExchangeNonlocalEdges(const std::vector<SInt> &vertex_distribution);
    Edgelist GenerateBadHyperbolicEdges(PEID chunk);

    PGeneratorConfig config_;

    PEID next_streaming_chunk_;
    PEID streaming_chunks_per_pe_;

    PEID     rank_;
    PEID     size_;
    MPI_Comm comm_;

    std::unique_ptr<GeneratorFactory> factory_;

    std::vector<VertexRange> my_vertex_ranges_;
    std::vector<Edgelist> nonlocal_edges_;
};
} // namespace kagen
