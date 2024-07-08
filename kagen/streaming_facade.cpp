#include "kagen/streaming_facade.h"

#include "kagen/definitions.h"
#include "kagen/factories.h"

#include <mpi.h>

#include <cassert>
#include <iomanip>
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
    my_vertex_ranges_.resize(streaming_chunks_per_pe_);

    // Special treatment for RHG: we lack inward-facing inter-chunk edges between annulis
    // To handle this in the streaming setting, we therefore generate the graph twice, only tracking these difficult
    // edges during the first iteration, and then use communication to fix them up
    if (config_.generator == GeneratorType::RHG) {
        if (rank_ == ROOT && !config_.quiet) {
            std::cout << "Hyperbolic generator requires two passes to generate the graph" << std::endl;
            std::cout << "Initializaing " << std::flush;
        }

        std::vector<SInt> vertex_distribution(size_ + 1);

        SInt max_nonlocal_edges = 0; // PE-level max. nonlocal edges
        SInt num_local_edges    = 0; // Total number of local edges
        SInt max_local_edges    = 0; // Chunk-level max local edges

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
            max_nonlocal_edges += nonlocal_edges_[chunk].size();
            num_local_edges += graph.edges.size();
            max_local_edges = std::max(max_local_edges, static_cast<SInt>(graph.edges.size()));

            if (rank_ == ROOT && !config_.quiet) {
                std::cout << "." << std::flush;
            }
        }

        MPI_Allreduce(MPI_IN_PLACE, &max_nonlocal_edges, 1, KAGEN_MPI_SINT, MPI_MAX, comm_);
        MPI_Allreduce(MPI_IN_PLACE, &num_local_edges, 1, KAGEN_MPI_SINT, MPI_SUM, comm_);
        MPI_Allreduce(MPI_IN_PLACE, &max_local_edges, 1, KAGEN_MPI_SINT, MPI_MAX, comm_);

        vertex_distribution[rank_]     = my_vertex_ranges_.front().first;
        vertex_distribution[rank_ + 1] = my_vertex_ranges_.back().second;
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, vertex_distribution.data() + 1, 1, KAGEN_MPI_SINT, comm_);

        if (rank_ == ROOT && !config_.quiet) {
            std::cout << std::endl;
            std::cout << "Total number of local edges:      " << num_local_edges << std::endl;
            std::cout << "Maximum number of local edges:    " << max_local_edges << " = " << std::fixed
                      << std::setprecision(3) << 16 * max_local_edges / 1024.0 / 1024.0 << " MB" << std::endl;
            std::cout << "Maximum number of nonlocal edges: " << max_nonlocal_edges << " = " << std::fixed
                      << std::setprecision(3) << 16 * max_nonlocal_edges / 1024.0 / 1024.0 << " MB" << std::endl;
            std::cout << "Memory peak:                      roughly "
                      << 16 * (max_local_edges + max_nonlocal_edges) / 1024.0 / 1024.0 << " MB" << std::endl;
        }

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
                ++send_counts[pe];
                ++idx;
            }
        }
    }

    std::exclusive_scan(send_counts.begin(), send_counts.end(), send_displs.begin(), 0);
    Edgelist send_bufs(send_displs.back() + send_counts.back());
    std::fill(send_counts.begin(), send_counts.end(), 0);

    for (auto& edges: nonlocal_edges_) {
        std::size_t idx = 0;

        for (PEID pe = 0; pe < size_; ++pe) {
            while (idx < edges.size() && edges[idx].first < vertex_distribution[pe + 1]) {
                send_bufs[send_displs[pe] + send_counts[pe]] = edges[idx];
                ++send_counts[pe];
                ++idx;
            }
        }

        edges.clear();
    }

    std::vector<int> recv_counts(size_);
    std::vector<int> recv_displs(size_);

    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, comm_);
    std::exclusive_scan(recv_counts.begin(), recv_counts.end(), recv_displs.begin(), 0);

    MPI_Datatype sint_pair = MPI_DATATYPE_NULL;
    MPI_Type_contiguous(2, KAGEN_MPI_SINT, &sint_pair);
    MPI_Type_commit(&sint_pair);

    Edgelist recv_bufs(recv_displs.back() + recv_counts.back());
    MPI_Alltoallv(
        send_bufs.data(), send_counts.data(), send_displs.data(), sint_pair, recv_bufs.data(), recv_counts.data(),
        recv_displs.data(), sint_pair, comm_);

    MPI_Type_free(&sint_pair);

    for (PEID pe = 0; pe < size_; ++pe) {
        std::sort(recv_bufs.begin() + recv_displs[pe], recv_bufs.begin() + recv_displs[pe] + recv_counts[pe]);
        PEID chunk = 0;

        for (PEID idx = recv_displs[pe]; idx < recv_displs[pe] + recv_counts[pe]; ++idx) {
            const auto [from, to] = recv_bufs[idx];

            assert(from >= vertex_distribution[rank_] && from < vertex_distribution[rank_ + 1]);

            while (chunk < static_cast<int>(my_vertex_ranges_.size()) && from >= my_vertex_ranges_[chunk].second) {
                chunk += 1;
            }

            assert(chunk < static_cast<int>(my_vertex_ranges_.size()));
            assert(from >= my_vertex_ranges_[chunk].first);
            assert(from < my_vertex_ranges_[chunk].second);

            nonlocal_edges_.at(chunk).emplace_back(from, to);
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
