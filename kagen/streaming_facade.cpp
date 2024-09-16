#include "kagen/streaming_facade.h"

#include "kagen/definitions.h"
#include "kagen/factories.h"
#include "kagen/tools/utils.h"

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

VertexRange StreamingGenerator::EstimateVertexRange(PEID pe) const {
    if (pe < 0) {
        pe = rank_;
    }

    if (!vertex_distribution_.empty()) {
        return {vertex_distribution_[pe], vertex_distribution_[pe + 1]};
    }

    return ComputeRange(config_.n, size_, pe);
}

void StreamingGenerator::Initialize() {
    nonlocal_edges_.clear();
    nonlocal_edges_.resize(streaming_chunks_per_pe_);

    my_vertex_ranges_.clear();
    my_vertex_ranges_.resize(streaming_chunks_per_pe_);

    if (!initialized_) {
        if (rank_ == ROOT && !config_.quiet) {
            std::cout << "Initializaing " << std::flush;
        }

        initialized_ = true;
        vertex_distribution_.resize(size_ + 1);

        SInt max_nonlocal_edges = 0; // PE-level max. nonlocal edges
        SInt num_nonlocal_edges = 0; // Total number of nonlocal edges
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

            // Some generators only report a meaningful vertex range if there is at least one vertex in the chunk.
            // Otherwise, it reports SInt max() for both first and second. For a nicer interface, fix the range.
            // @todo: this assumes that there is at least one none-empty chunk on each PE ...
            if (chunk > 0 && my_vertex_ranges_[chunk].first == std::numeric_limits<SInt>::max()) {
                my_vertex_ranges_[chunk].first = my_vertex_ranges_[chunk].second = my_vertex_ranges_[chunk - 1].second;
            } else if (chunk > 0 && my_vertex_ranges_[chunk].first != std::numeric_limits<SInt>::max()) {
                for (PEID prev_chunk = chunk - 1;
                     prev_chunk >= 0 && my_vertex_ranges_[prev_chunk].first == std::numeric_limits<SInt>::max();
                     --prev_chunk) {
                    my_vertex_ranges_[prev_chunk].first = my_vertex_ranges_[prev_chunk].second =
                        my_vertex_ranges_[chunk].first;
                }
            }

            nonlocal_edges_[chunk] = std::move(nonlocal_edges);
            max_nonlocal_edges += nonlocal_edges_[chunk].size();
            num_nonlocal_edges += nonlocal_edges_[chunk].size();
            num_local_edges += graph.edges.size();
            max_local_edges = std::max(max_local_edges, static_cast<SInt>(graph.edges.size()));

            if (rank_ == ROOT && !config_.quiet) {
                std::cout << "." << std::flush;
            }
        }

        MPI_Allreduce(MPI_IN_PLACE, &max_nonlocal_edges, 1, KAGEN_MPI_SINT, MPI_MAX, comm_);
        MPI_Allreduce(MPI_IN_PLACE, &num_nonlocal_edges, 1, KAGEN_MPI_SINT, MPI_SUM, comm_);
        MPI_Allreduce(MPI_IN_PLACE, &num_local_edges, 1, KAGEN_MPI_SINT, MPI_SUM, comm_);
        MPI_Allreduce(MPI_IN_PLACE, &max_local_edges, 1, KAGEN_MPI_SINT, MPI_MAX, comm_);

        vertex_distribution_[rank_]     = my_vertex_ranges_.front().first;
        vertex_distribution_[rank_ + 1] = my_vertex_ranges_.back().second;
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, vertex_distribution_.data() + 1, 1, KAGEN_MPI_SINT, comm_);

        if (rank_ == ROOT && !config_.quiet) {
            std::cout << std::endl;
            std::cout << "Total number of local edges:      " << num_local_edges << std::endl;
            std::cout << "Maximum number of local edges:    " << max_local_edges << " = " << std::fixed
                      << std::setprecision(3) << 16 * max_local_edges / 1024.0 / 1024.0 << " MB" << std::endl;
            std::cout << "Total number of nonlocal edges:   " << num_nonlocal_edges << " = " << std::fixed
                      << std::setprecision(3) << 16 * num_nonlocal_edges / 1024.0 / 1024.0 << " MB" << std::endl;
            std::cout << "Maximum number of nonlocal edges: " << max_nonlocal_edges << " = " << std::fixed
                      << std::setprecision(3) << 16 * max_nonlocal_edges / 1024.0 / 1024.0 << " MB" << std::endl;
            std::cout << "Memory peak:                      roughly "
                      << 16 * (max_local_edges + max_nonlocal_edges) / 1024.0 / 1024.0 << " MB" << std::endl;
        }

        ExchangeNonlocalEdges();
    }

    next_streaming_chunk_ = 0;
}

void StreamingGenerator::ExchangeNonlocalEdges() {
    std::vector<int> send_counts(size_);
    std::vector<int> send_displs(size_);

    for (const auto& edges: nonlocal_edges_) {
        std::size_t idx = 0;

        for (PEID pe = 0; pe < size_; ++pe) {
            while (idx < edges.size() && edges[idx].first < vertex_distribution_[pe + 1]) {
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
            while (idx < edges.size() && edges[idx].first < vertex_distribution_[pe + 1]) {
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

            assert(from >= vertex_distribution_[rank_] && from < vertex_distribution_[rank_ + 1]);

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

    Graph graph = CreateGenerator(next_streaming_chunk_)->Generate(GraphRepresentation::EDGE_LIST)->Take();
    if (graph.edges.empty()) {
        ++next_streaming_chunk_;
        return Next();
    }

    StreamedGraph sgraph = {
        .vertex_range    = graph.vertex_range,
        .primary_edges   = std::move(graph.edges),
        .secondary_edges = std::move(nonlocal_edges_[next_streaming_chunk_]),
    };

    if (!std::is_sorted(sgraph.primary_edges.begin(), sgraph.primary_edges.end())) {
        std::sort(sgraph.primary_edges.begin(), sgraph.primary_edges.end());
    }
    if (!std::is_sorted(sgraph.secondary_edges.begin(), sgraph.secondary_edges.end())) {
        std::sort(sgraph.secondary_edges.begin(), sgraph.secondary_edges.end());
    }

    ++next_streaming_chunk_;

    return sgraph;
}

bool StreamingGenerator::Continue() const {
    return next_streaming_chunk_ < streaming_chunks_per_pe_;
}
} // namespace kagen
