#include "kagen/streaming_facade.h"

#include "kagen/factories.h"

#include <mpi.h>

#include <numeric>

namespace kagen {
StreamingGenerator::StreamingGenerator(const std::string& options, const PEID chunks_per_pe, MPI_Comm comm)
    : config_(CreateConfigFromString(options)),
      comm_(comm),
      factory_(CreateGeneratorFactory(config_.generator)) {
    MPI_Comm_size(comm_, &size_);
    MPI_Comm_rank(comm_, &rank_);

    const PEID streaming_size = chunks_per_pe * size_;
    const PEID streaming_rank = chunks_per_pe * rank_;

    config_                  = factory_->NormalizeParameters(config_, streaming_rank, streaming_size, rank_ == 0);
    next_streaming_chunk_    = 0;
    streaming_chunks_per_pe_ = config_.k / size_;
}

void StreamingGenerator::Initialize() {
    nonlocal_edges_.clear();
    nonlocal_edges_.resize(streaming_chunks_per_pe_);

    my_vertex_ranges_.clear();
    my_vertex_ranges_.reserve(streaming_chunks_per_pe_);

    // Special treatment for RHG: we lack inward-facing inter-chunk edges between annulis
    // To handle this in the streaming setting, we therefore generate the graph twice, only tracking these difficult
    // edges during the first iteration, and then use communication to fix them up
    if (config_.generator == GeneratorType::RHG) {
        std::vector<SInt> vertex_distribution(size_ + 1);

        for (PEID chunk = 0; chunk < streaming_chunks_per_pe_; ++chunk) {
            auto generator = CreateGenerator(chunk);
            generator->Generate(GraphRepresentation::EDGE_LIST);

            auto  nonlocal_edges = generator->TakeNonlocalEdges();
            Graph graph          = generator->Take();

            std::sort(nonlocal_edges.begin(), nonlocal_edges.end(), [](const auto& lhs, const auto& rhs) {
                return lhs.first < rhs.first;
            });

            my_vertex_ranges_[chunk] = graph.vertex_range;
            nonlocal_edges_[chunk]   = std::move(nonlocal_edges);
        }

        vertex_distribution[rank_]     = my_vertex_ranges_.front().first;
        vertex_distribution[rank_ + 1] = my_vertex_ranges_.back().second;
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, vertex_distribution.data() + 1, 1, KAGEN_MPI_SINT, comm_);

        ExchangeNonlocalEdges(vertex_distribution);
    }

    next_streaming_chunk_ = 0;
}

void StreamingGenerator::ExchangeNonlocalEdges(const std::vector<SInt>& vertex_distribution) {
    std::vector<int> send_counts(size_);
    std::vector<int> send_displs(size_);

    for (const auto& edges: nonlocal_edges_) {
        std::size_t idx = 0;

        for (PEID pe = 0; pe < size_; ++pe) {
            while (idx < edges.size() && edges[idx].first < vertex_distribution[pe + 1]) {
                send_counts[pe] += 2;
                idx += 1;
            }
        }
    }

    std::exclusive_scan(send_counts.begin(), send_counts.end(), send_displs.begin(), 0);
    std::vector<SInt> send_bufs(send_displs.back() + send_counts.back());
    std::fill(send_counts.begin(), send_counts.end(), 0);

    for (auto& edges: nonlocal_edges_) {
        std::size_t idx = 0;

        for (PEID pe = 0; pe < size_; ++pe) {
            while (idx < edges.size() && edges[idx].first < vertex_distribution[pe + 1]) {
                send_bufs[send_displs[pe] + send_counts[pe]++] = edges[idx].first;
                send_bufs[send_displs[pe] + send_counts[pe]++] = edges[idx].second;
                idx += 1;
            }
        }

        edges.clear();
    }

    std::vector<int> recv_counts(size_);
    std::vector<int> recv_displs(size_);

    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, comm_);
    std::exclusive_scan(recv_counts.begin(), recv_counts.end(), recv_displs.begin(), 0);

    std::vector<SInt> recv_bufs(recv_displs.back() + recv_counts.back());

    MPI_Alltoallv(
        send_bufs.data(), send_counts.data(), send_displs.data(), KAGEN_MPI_SINT, recv_bufs.data(), recv_counts.data(),
        recv_displs.data(), KAGEN_MPI_SINT, comm_);

    for (PEID pe = 0; pe < size_; ++pe) {
        std::sort(recv_bufs.begin() + recv_displs[pe], recv_bufs.begin() + recv_displs[pe] + recv_counts[pe]);
        PEID chunk = 0;

        for (PEID idx = recv_displs[pe]; idx < recv_displs[pe] + recv_counts[pe]; idx += 2) {
            const SInt from = recv_bufs[idx];
            const SInt to   = recv_bufs[idx + 1];

            while (from >= my_vertex_ranges_[chunk].second) {
                chunk += 1;
            }

            nonlocal_edges_[chunk].emplace_back(from, to);
        }
    }
}

std::unique_ptr<Generator> StreamingGenerator::CreateGenerator(const PEID chunk) {
    return factory_->Create(config_, streaming_chunks_per_pe_ * rank_ + chunk, streaming_chunks_per_pe_ * size_);
}

StreamedGraph StreamingGenerator::Next() {
    if (next_streaming_chunk_ >= streaming_chunks_per_pe_) {
        return {
            .vertex_range    = {-1, -1},
            .primary_edges   = {},
            .secondary_edges = {},
        };
    }

    Graph         graph  = CreateGenerator(next_streaming_chunk_)->Generate(GraphRepresentation::EDGE_LIST)->Take();
    StreamedGraph sgraph = {
        .vertex_range    = graph.vertex_range,
        .primary_edges   = std::move(graph.edges),
        .secondary_edges = std::move(nonlocal_edges_[next_streaming_chunk_]),
    };
    sgraph.SortEdgelist();

    ++next_streaming_chunk_;

    return sgraph;
}

bool StreamingGenerator::Continue() const {
    return next_streaming_chunk_ < streaming_chunks_per_pe_;
}
} // namespace kagen
