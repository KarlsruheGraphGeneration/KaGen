#include "kagen/skagen.h"

#include "kagen/factories.h"

#include <mpi.h>

namespace kagen {
sKaGen::sKaGen(const std::string& options, int chunks, MPI_Comm comm)
    : config_(CreateConfigFromString(options)),
      comm_(comm),
      factory_(CreateGeneratorFactory(config_.generator)) {
    MPI_Comm_size(comm_, &size_);
    MPI_Comm_rank(comm_, &rank_);

    // Let the generator decide the actual number of chunks:
    config_.k      = chunks * size_;
    config_        = factory_->NormalizeParameters(config_, rank_, size_, rank_ == 0);
    next_chunk_    = 0;
    chunks_per_pe_ = config_.k / size_;
    // This value should be >= chunks
}

void sKaGen::Initialize() {
    nonlocal_edges_.clear();
    nonlocal_edges_.resize(chunks_per_pe_);

    // Special treatment for RHG: we lack inward-facing inter-chunk edges between annulis
    // To handle this in the streaming setting, we therefore generate the graph twice, only tracking these difficult
    // edges during the first iteration, and then use communication to fix them up
    if (config_.generator == GeneratorType::RHG) {
        for (PEID chunk = 0; chunk < chunks_per_pe_; ++chunk) {
            auto generator = CreateGenerator(chunk);
            generator->Generate(GraphRepresentation::EDGE_LIST);
            nonlocal_edges_[chunk] = generator->TakeNonlocalEdges();
        }
    }

    ExchangeNonlocalEdges();

    next_chunk_ = 0;
}

void sKaGen::ExchangeNonlocalEdges() {}

std::unique_ptr<Generator> sKaGen::CreateGenerator(const PEID chunk) {
    return factory_->Create(config_, chunks_per_pe_ * rank_ + chunk, chunks_per_pe_ * size_);
}

bool sKaGen::Continue(EdgelistChunk& chunk) {
    if (next_chunk_ >= chunks_per_pe_) {
        return false;
    }

    Graph graph_chunk     = CreateGenerator(next_chunk_)->Generate(GraphRepresentation::EDGE_LIST)->Take();
    chunk.vertex_range    = graph_chunk.vertex_range;
    chunk.primary_edges   = std::move(graph_chunk.edges);
    chunk.secondary_edges = std::move(nonlocal_edges_[next_chunk_]);

    return (++next_chunk_) < chunks_per_pe_;
}
} // namespace kagen
