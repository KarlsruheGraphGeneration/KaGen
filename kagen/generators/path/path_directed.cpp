#include "kagen/generators/path/path_directed.h"

#include "kagen/sampling/hash.hpp"
#include "kagen/tools/utils.h"

#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <exception>
#include <numeric>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef KAGEN_XXHASH_FOUND
    #include "kagen/tools/random_permutation.h"
#endif // KAGEN_XXHASH_FOUND

namespace kagen {
namespace {

class Distribution {
public:
    Distribution() = default;
    // Construct from local size; will perform an MPI_Allgather to build offsets
    explicit Distribution(SInt local_size, MPI_Comm comm = MPI_COMM_WORLD) {
        int rank = 0, comm_size = 0;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &comm_size);

        offset_.assign(static_cast<std::size_t>(comm_size) + 1u, SInt{0});
        offset_[static_cast<std::size_t>(rank)] = local_size;

        MPI_Allgather(
            &offset_[static_cast<std::size_t>(rank)], 1, MPI_UNSIGNED_LONG_LONG, offset_.data(), 1,
            MPI_UNSIGNED_LONG_LONG, comm);

        std::exclusive_scan(offset_.begin(), offset_.end(), offset_.begin(), SInt{0});
    }

    [[nodiscard]] auto get_owner(SInt idx) const -> std::size_t {
        SInt v = static_cast<SInt>(idx);
        // upper_bound returns iterator to first element > v
        auto it = std::upper_bound(offset_.begin(), offset_.end(), v);
        if (it == offset_.begin() || it == offset_.end()) {
            throw std::out_of_range("Index out of bounds in Distribution::get_owner");
        }
        return static_cast<std::size_t>(std::distance(offset_.begin(), it)) - 1u;
    }

    [[nodiscard]] auto get_owner_signed(SInt idx) const -> int {
        return static_cast<int>(get_owner(idx));
    }

    [[nodiscard]] auto counts() const -> std::vector<SInt> {
        std::vector<SInt> result;
        if (offset_.size() < 2u)
            return result;
        result.reserve(offset_.size() - 1u);
        for (std::size_t i = 0u; i + 1u < offset_.size(); ++i) {
            result.push_back(offset_[i + 1u] - offset_[i]);
        }
        return result;
    }

    [[nodiscard]] auto index_range_begin(std::size_t rank) const -> SInt {
        if (rank + 1u >= offset_.size()) {
            throw std::out_of_range("rank out of range in index_range_begin");
        }
        return offset_[rank];
    }

    [[nodiscard]] auto index_range_end(std::size_t rank) const -> SInt {
        if (rank + 1u >= offset_.size()) {
            throw std::out_of_range("rank out of range in index_range_end");
        }
        return offset_[rank + 1u];
    }

    [[nodiscard]] auto is_local(SInt idx, std::size_t rank) const -> bool {
        SInt v = static_cast<SInt>(idx);
        return index_range_begin(rank) <= v && v < index_range_end(rank);
    }

    [[nodiscard]] auto get_local_idx(SInt idx, std::size_t rank) const -> SInt {
        if (!is_local(idx, rank)) {
            throw std::out_of_range("get_local_idx: idx is not local to rank");
        }
        return static_cast<SInt>(static_cast<SInt>(idx) - index_range_begin(rank));
    }

    [[nodiscard]] auto get_global_idx(SInt local_idx, std::size_t rank) const -> SInt {
        return static_cast<SInt>(static_cast<SInt>(local_idx) + index_range_begin(rank));
    }

    [[nodiscard]] auto get_global_size() const -> SInt {
        return offset_.back();
    }

    [[nodiscard]] auto get_local_size(std::size_t rank) const -> SInt {
        return index_range_end(rank) - index_range_begin(rank);
    }

    [[nodiscard]] auto local_index_range(std::size_t rank) const -> std::pair<SInt, SInt> {
        SInt b = SInt{0};
        SInt e = static_cast<SInt>(get_local_size(rank));
        return {b, e};
    }

    [[nodiscard]] auto global_index_range(std::size_t rank) const -> std::pair<SInt, SInt> {
        SInt b = static_cast<SInt>(index_range_begin(rank));
        SInt e = static_cast<SInt>(index_range_end(rank));
        return {b, e};
    }

private:
    std::vector<SInt> offset_;
};

} // namespace

std::unique_ptr<Generator>
PathDirectedFactory::Create(const PGeneratorConfig& config, const PEID rank, const PEID size) const {
    return std::make_unique<PathDirected>(config, rank, size);
}

PGeneratorConfig PathDirectedFactory::NormalizeParameters(PGeneratorConfig config, PEID, PEID, bool) const {
#ifndef KAGEN_XXHASH_FOUND
    if (config.permute) {
        throw ConfigurationError("path permutation requires xxHash, but build was configured without xxHash");
    }
#endif // KAGEN_XXHASH_FOUND
    if (config.p == 0.0) {
        config.p = 1.0;
    }
    if (config.p > 1.0) {
        throw ConfigurationError("edge probability p must be in [0, 1]");
    }
    return config;
}

template <typename Request, typename Reply, typename SourceRankFn, typename MakeReplyFn>
auto request_reply(
    std::unordered_map<PEID, std::vector<Request>>&& send_bufs_requests, SourceRankFn const& get_source_rank,
    MakeReplyFn const& make_reply, MPI_Comm comm) {
    MPI_Datatype request_mpi_type;
    MPI_Datatype reply_mpi_type;
    MPI_Type_contiguous(sizeof(Request), MPI_BYTE, &request_mpi_type);
    MPI_Type_contiguous(sizeof(Reply), MPI_BYTE, &reply_mpi_type);
    MPI_Type_commit(&request_mpi_type);
    MPI_Type_commit(&reply_mpi_type);
    auto recv_requests = ExchangeMessageBuffers(std::move(send_bufs_requests), request_mpi_type, comm);
    std::unordered_map<PEID, std::vector<Reply>> send_bufs_replies;
    for (const auto& req: recv_requests) {
        PEID source_rank = get_source_rank(req);
        send_bufs_replies[source_rank].push_back(make_reply(req));
    }
    auto recv_replies = ExchangeMessageBuffers(std::move(send_bufs_replies), reply_mpi_type, comm);
    MPI_Type_free(&request_mpi_type);
    MPI_Type_free(&reply_mpi_type);
    return recv_replies;
}

struct PartialPermutator {
    SInt to_permutation_global(SInt global_index) {
        assert(path_distribution_.is_local(global_index, rank_) && is_permuted(global_index));
        auto it = to_local_permutation_index_map_.find(global_index);
        assert(it != to_local_permutation_index_map_.end());
        return permutation_distribution_.get_global_idx(it->second, rank_);
    }
    SInt from_permutation_local_to_global(SInt permutation_local_index) {
        assert(permutation_local_index < permutation_vertices_.size());
        return permutation_vertices_[permutation_local_index];
    }
    SInt from_permutation_to_global(SInt permutation_index) {
        assert(permutation_distribution_.is_local(permutation_index, rank_));
        auto permutation_local_index = permutation_distribution_.get_local_idx(permutation_index, rank_);
        return from_permutation_local_to_global(permutation_local_index);
    }
    bool is_permuted(SInt global_index) {
        assert(path_distribution_.is_local(global_index, rank_));
        return is_in_permutation.find(global_index) != is_in_permutation.end();
    }
    template <typename Fn>
    void fill_map(std::unordered_map<SInt, SInt>& map, Fn const& fn) {
        struct Request {
            SInt global_index_argument;
            SInt permutation_index_value;
            PEID origin;
        };
        struct Reply {
            SInt global_index_argument;
            SInt permutation_index_value;
            SInt global_index_value;
        };
        std::unordered_map<PEID, std::vector<Request>> send_bufs_requests;
        for (SInt i = range_.first; i < range_.second; ++i) {
            if (is_permuted(i)) {
                Request request;
                request.global_index_argument   = i;
                request.origin                  = rank_;
                request.permutation_index_value = fn(to_permutation_global(i));
                auto target                     = permutation_distribution_.get_owner(request.permutation_index_value);
                send_bufs_requests[target].push_back(request);
            }
        }
        auto get_source_rank = [](const Request& request) {
            return request.origin;
        };
        auto make_reply = [&](const Request& request) {
            Reply reply;
            reply.global_index_argument   = request.global_index_argument;
            reply.permutation_index_value = request.permutation_index_value;
            reply.global_index_value      = from_permutation_to_global(request.permutation_index_value);
            return reply;
        };
        auto recv_replies =
            request_reply<Request, Reply>(std::move(send_bufs_requests), get_source_rank, make_reply, comm_);
        for (const auto [global_index_argument, _, global_index_value]: recv_replies) {
            map[global_index_argument] = global_index_value;
        }
    }
    void fill_f_permuted(const auto& permutator) {
        fill_map(f_map, [&](SInt x) { return permutator.f(x); });
    }
    void fill_finv_permuted(const auto& permutator) {
        fill_map(finv_map, [&](SInt x) { return permutator.finv(x); });
    }
    PartialPermutator(SInt n, double permutation_ratio, int seed, VertexRange range, MPI_Comm comm)
        : comm_{comm},
          range_{range} {
        MPI_Comm_rank(comm, &rank_);
        auto              ranges     = AllgatherVertexRange(range, comm);
        SInt              local_size = range.second - range.first;
        std::vector<SInt> local_chunk(local_size);
        std::iota(local_chunk.begin(), local_chunk.end(), range.first);
        SInt            local_permutation_size = static_cast<SInt>(permutation_ratio * static_cast<double>(local_size));
        std::mt19937_64 gen;
        std::sample(
            local_chunk.begin(), local_chunk.end(), std::back_inserter(permutation_vertices_), local_permutation_size,
            gen);
        for (const auto& v: permutation_vertices_) {
            is_in_permutation.insert(v);
        }
        path_distribution_        = Distribution{local_size, comm};
        permutation_distribution_ = Distribution{local_permutation_size, comm};

        for (SInt i = 0; i < permutation_vertices_.size(); ++i) {
            to_local_permutation_index_map_.emplace(permutation_vertices_[i], i);
        }

        auto permutator = random_permutation::FeistelPseudoRandomPermutation::buildPermutation(
            permutation_distribution_.get_global_size() - 1, static_cast<std::uint64_t>(seed));
        fill_f_permuted(permutator);
        fill_finv_permuted(permutator);
        for (SInt i = range_.first; i < range_.second; ++i) {
            if (!is_permuted(i)) {
                assert(f_map.find(i) == f_map.end());
                assert(finv_map.find(i) == finv_map.end());
                f_map[i]    = i;
                finv_map[i] = i;
            }
        }

        struct Request {
            SInt request;
            PEID origin;
        };
        struct Reply {
            SInt request;
            SInt f_request;
        };
        std::unordered_map<PEID, std::vector<Request>> send_bufs_requests;
        for (const auto& [x, finv]: finv_map) {
            auto request = (finv + 1) % n;
            if (f_map.find(request) == f_map.end()) {
                auto target = FindPEInRange(request, ranges);
                send_bufs_requests[target].emplace_back(request, rank_);
            }
        }
        auto recv_replies = request_reply<Request, Reply>(
            std::move(send_bufs_requests), [](const Request& req) { return req.origin; },
            [&](const Request& req) {
                auto it = f_map.find(req.request);
                assert(path_distribution_.is_local(req.request, rank_));
                assert(it != f_map.end());
                return Reply{req.request, it->second};
            },
            comm_);

        for (const auto& [request, f_request]: recv_replies) {
            f_map.emplace(request, f_request);
        }
    }
    SInt f(SInt x) {
        auto it = f_map.find(x);
        if (it == f_map.end()) {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            throw std::runtime_error("wrong lookup in f " + std::to_string(x) + " " + std::to_string(rank));
        }
        return it->second;
    }
    SInt finv(SInt x) {
        auto it = finv_map.find(x);
        if (it == finv_map.end()) {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            throw std::runtime_error("wrong lookup in finv " + std::to_string(x) + " " + std::to_string(rank));
        }
        return it->second;
    }
    MPI_Comm                       comm_;
    Distribution                   path_distribution_;
    Distribution                   permutation_distribution_;
    VertexRange                    range_;
    VertexRange                    permutation_range_;
    std::vector<VertexRange>       permutation_ranges_;
    PEID                           rank_;
    std::unordered_set<SInt>       is_in_permutation;
    std::unordered_map<SInt, SInt> to_local_permutation_index_map_;
    std::vector<SInt>              permutation_vertices_;
    std::unordered_map<SInt, SInt> f_map;
    std::unordered_map<SInt, SInt> finv_map;
};

class Permutator {
public:
    Permutator(PGeneratorConfig const& config, VertexRange range, MPI_Comm) {
        if (config.perm_p < 1.0) {
            partial_permutator  = PartialPermutator(config.n, config.perm_p, config.seed, range, MPI_COMM_WORLD);
            use_partial_permute = true;
        } else {
            complete_permutator = random_permutation::FeistelPseudoRandomPermutation::buildPermutation(
                config.n - 1, static_cast<std::uint64_t>(config.seed));
        }
    }
    SInt f(SInt x) {
        if (use_partial_permute) {
            return partial_permutator->f(x);
        } else {
            return complete_permutator->f(x);
        }
    }
    SInt finv(SInt x) {
        if (use_partial_permute) {
            return partial_permutator->finv(x);
        } else {
            return complete_permutator->finv(x);
        }
    }

private:
    bool                                                              use_partial_permute{false};
    std::optional<PartialPermutator>                                  partial_permutator;
    std::optional<random_permutation::FeistelPseudoRandomPermutation> complete_permutator;
};

PathDirected::PathDirected(const PGeneratorConfig& config, const PEID rank, const PEID size)
    : config_(config),
      rank_(rank),
      size_(size),
      rng_(config) {}

void PathDirected::GenerateEdgeList() {
    if (config_.n <= 1) {
        return;
    }

    // all PE get n / size nodes
    // the first n modulo size PEs obtain one additional node.
    const SInt nodes_per_pe                = config_.n / size_;
    const PEID num_pe_with_additional_node = config_.n % size_;
    const bool has_pe_additional_node      = rank_ < num_pe_with_additional_node;
    const SInt begin_nodes                 = std::min(num_pe_with_additional_node, rank_) + rank_ * nodes_per_pe;
    const SInt end_nodes                   = begin_nodes + nodes_per_pe + has_pe_additional_node;
    // this is an exclusive range
    SetVertexRange(begin_nodes, end_nodes);

#ifdef KAGEN_XXHASH_FOUND
    // auto permutator = random_permutation::FeistelPseudoRandomPermutation::buildPermutation(config_.n - 1, 0);
    Permutator permutator(config_, graph_.vertex_range, MPI_COMM_WORLD);
#endif // KAGEN_XXHASH_FOUND

    for (SInt i = begin_nodes; i < end_nodes; ++i) {
        const auto [j, is_valid] = [&]() -> std::pair<SInt, bool> {
#ifdef KAGEN_XXHASH_FOUND
            if (config_.permute) {
                SInt j = (permutator.finv(i) + 1) % config_.n;
                return std::make_pair(permutator.f(j), j != 0 || config_.periodic);
            }
#endif // KAGEN_XXHASH_FOUND

            const SInt j = (i + 1) % config_.n;
            return std::make_pair(j, j != 0 || config_.periodic);
        }();

        if (is_valid) {
            SInt edge_seed = std::min(i, j) * config_.n + std::max(i, j);
            SInt h         = sampling::Spooky::hash(config_.seed + edge_seed);
            if (rng_.GenerateBinomial(h, 1, config_.p)) {
                PushEdge(i, j);
            }
        }
    }
}

} // namespace kagen
