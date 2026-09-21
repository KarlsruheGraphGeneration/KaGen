#include "kagen/generators/generator.h"

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/edgeweight_generators/default_generator.h"
#include "kagen/edgeweight_generators/edge_weight_generator.h"
#include "kagen/edgeweight_generators/euclidean_distance_generator.h"
#include "kagen/edgeweight_generators/hashing_based_generator.h"
#include "kagen/edgeweight_generators/uniform_random_generator.h"
#include "kagen/edgeweight_generators/voiding_generator.h"
#include "kagen/kagen.h"
#include "kagen/tools/converter.h"
#include "kagen/tools/postprocessor.h"
#include "kagen/vertexweight_generators/default_generator.h"
#include "kagen/vertexweight_generators/uniform_random_generator.h"
#include "kagen/vertexweight_generators/vertex_weight_generator.h"
#include "kagen/vertexweight_generators/voiding_generator.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <optional>
#include <sstream>

#ifdef KAGEN_XXHASH_FOUND
    #include "kagen/tools/random_permutation.h"
#endif

namespace kagen {
Generator::~Generator() = default;

Generator* Generator::Generate(const GraphRepresentation representation) {
    Reset();
    desired_representation_ = representation;

    switch (desired_representation_) {
        case GraphRepresentation::EDGE_LIST:
            GenerateEdgeList();
            break;

        case GraphRepresentation::CSR:
            GenerateCSR();
            break;
    }

    return this;
}

Generator* Generator::Finalize(MPI_Comm comm) {
    switch (desired_representation_) {
        case GraphRepresentation::EDGE_LIST:
            FinalizeEdgeList(comm);
            break;

        case GraphRepresentation::CSR:
            FinalizeCSR(comm);
            break;
    }

    graph_.representation = desired_representation_;

    return this;
}

std::unique_ptr<kagen::EdgeWeightGenerator>
CreateEdgeWeightGenerator(const EdgeWeightConfig weight_config, MPI_Comm comm, const VertexRange vertex_range) {
    switch (weight_config.generator_type) {
        case EdgeWeightGeneratorType::DEFAULT:
            return std::make_unique<DefaultEdgeWeightGenerator>(weight_config);
        case EdgeWeightGeneratorType::VOIDING:
            return std::make_unique<VoidingEdgeWeightGenerator>(weight_config);
        case EdgeWeightGeneratorType::HASHING_BASED:
            return std::make_unique<HashingBasedEdgeWeightGenerator>(weight_config, vertex_range);
        case EdgeWeightGeneratorType::EUCLIDEAN_DISTANCE:
            return std::make_unique<EuclideanDistanceEdgeWeightGenerator>(weight_config);
        case EdgeWeightGeneratorType::UNIFORM_RANDOM:
            return std::make_unique<UniformRandomEdgeWeightGenerator>(weight_config, comm, vertex_range);
    }

    throw std::runtime_error("invalid weight generator type");
}

void Generator::GenerateEdgeWeights(EdgeWeightConfig weight_config, MPI_Comm comm) {
    std::unique_ptr<kagen::EdgeWeightGenerator> edge_weight_generator =
        CreateEdgeWeightGenerator(weight_config, comm, graph_.vertex_range);

    switch (desired_representation_) {
        case GraphRepresentation::EDGE_LIST:
            edge_weight_generator->GenerateEdgeWeights(graph_.edges, graph_.edge_weights);
            break;
        case GraphRepresentation::CSR:
            edge_weight_generator->GenerateEdgeWeights(graph_.xadj, graph_.adjncy, graph_.edge_weights);
            break;
    }
}

namespace {
template <typename Permutator>
auto ApplyPermutationAndComputeSendBuffersEdgeList(
    const Graph& graph, const std::vector<VertexRange>& recv_ranges, Permutator&& permute) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Edgelist edges = graph.edges;
    for (auto& [src, dst]: edges) {
        src = permute(src);
        dst = permute(dst);
    }
    std::unordered_map<PEID, std::vector<SInt>>  send_buffers;
    std::unordered_map<PEID, std::vector<SSInt>> edge_weights_send_buffers;
    std::unordered_map<PEID, std::vector<SSInt>> vertex_weights_send_buffers;

    bool has_vertex_weights = !graph.vertex_weights.empty();
    bool has_edge_weights   = graph.NumberOfLocalEdges() == 0 || !graph.edge_weights.empty();

    if (has_vertex_weights) {
        for (std::size_t i = 0; i < graph.NumberOfLocalVertices(); ++i) {
            const SInt          global_id          = graph.vertex_range.first + i;
            const SInt          permuted_global_id = permute(global_id);
            const SSInt         weight             = graph.vertex_weights[i];
            const PEID          target_pe          = FindPEInRange(permuted_global_id, recv_ranges);
            std::vector<SSInt>& weights_send_buf   = vertex_weights_send_buffers[target_pe];
            weights_send_buf.push_back(permuted_global_id);
            weights_send_buf.push_back(weight);
        }
    }

    for (std::size_t i = 0; i < edges.size(); ++i) {
        const auto& [src, dst]       = edges[i];
        const PEID         target_pe = FindPEInRange(src, recv_ranges);
        std::vector<SInt>& send_buf  = send_buffers[target_pe];
        send_buf.push_back(src);
        send_buf.push_back(dst);
        // [Permuted_Src_Id, Degree, EdgeWeights]
        if (has_edge_weights) {
            std::vector<SSInt>& weights_send_buf = edge_weights_send_buffers[target_pe];
            weights_send_buf.push_back(graph.edge_weights[i]);
        }
    }
    return std::make_tuple(
        std::move(send_buffers), std::move(edge_weights_send_buffers), std::move(vertex_weights_send_buffers));
}

template <typename Permutator>
auto ApplyPermutationAndComputeSendBuffersCSR(
    const Graph& graph, const std::vector<VertexRange>& recv_ranges, Permutator&& permute) {
    AdjncyArray permuted_adjncy = graph.adjncy;

    for (auto& edge: permuted_adjncy) {
        edge = permute(edge);
    }

    std::unordered_map<PEID, std::vector<SInt>>  send_buffers;
    std::unordered_map<PEID, std::vector<SSInt>> edge_weights_send_buffers;
    std::unordered_map<PEID, std::vector<SSInt>> vertex_weights_send_buffers;

    bool has_edge_weights   = graph.NumberOfLocalEdges() == 0 || !graph.edge_weights.empty();
    bool has_vertex_weights = !graph.vertex_weights.empty();

    // xadj/adjncy are indexed by the physically present row space, which is *not* vertex_range: on a PE holding
    // a replica of a split vertex (left_partial_vertex), that vertex has a row here even though it is credited
    // to -- and counted in the vertex_range of -- the lower-rank neighbor. Using vertex_range.first would label
    // every row with the ID of the next vertex, and run one row past the end of the range; on the last PE that
    // last row would ask for permute(n), which is outside the permutation's domain.
    const VertexRange row_range = graph.PhysicalVertexRange();

    for (std::size_t i = 0; i + 1 < graph.xadj.size(); ++i) {
        const SInt         degree             = graph.xadj[i + 1] - graph.xadj[i];
        const SInt         global_id          = row_range.first + i;
        const SInt         permuted_global_id = permute(global_id);
        const PEID         target_pe          = FindPEInRange(permuted_global_id, recv_ranges);
        std::vector<SInt>& send_buf           = send_buffers[target_pe];
        auto               edge_begin_offset  = graph.xadj[i];
        auto               edge_end_offset    = graph.xadj[i + 1];
        // [Permuted_Src_Id, Degree, [Permuted_Dst_Ids]] with #Permuted_Dst_Ids = Degree
        send_buf.push_back(permuted_global_id);
        send_buf.push_back(degree);
        send_buf.insert(
            send_buf.end(), permuted_adjncy.begin() + edge_begin_offset, permuted_adjncy.begin() + edge_end_offset);
        // [Permuted_Src_Id, Degree, EdgeWeights]
        if (has_edge_weights) {
            std::vector<SSInt>& weights_send_buf = edge_weights_send_buffers[target_pe];
            weights_send_buf.push_back(permuted_global_id);
            weights_send_buf.push_back(degree);
            weights_send_buf.insert(
                weights_send_buf.end(), graph.edge_weights.begin() + edge_begin_offset,
                graph.edge_weights.begin() + edge_end_offset);
        }
        // [Permuted_Src_Id, VertexWeight]
        // vertex_weights is indexed by vertex_range, which excludes a left-partial replica row: that vertex's
        // weight is held -- and sent -- by its canonical owner, the lower-rank neighbor, so skip it here rather
        // than reading past the end of the (shorter) weight array.
        if (has_vertex_weights && global_id >= graph.vertex_range.first && global_id < graph.vertex_range.second) {
            std::vector<SSInt>& vertex_weights_send_buf = vertex_weights_send_buffers[target_pe];
            vertex_weights_send_buf.push_back(permuted_global_id);
            vertex_weights_send_buf.push_back(graph.vertex_weights[global_id - graph.vertex_range.first]);
        }
    }
    return std::make_tuple(
        std::move(send_buffers), std::move(edge_weights_send_buffers), std::move(vertex_weights_send_buffers));
}

template <typename Permutator>
auto ApplyPermutationAndComputeSendBuffers(
    const Graph& graph, const std::vector<VertexRange>& recv_ranges, Permutator&& permutator) {
    switch (graph.representation) {
        case GraphRepresentation::EDGE_LIST:
            return ApplyPermutationAndComputeSendBuffersEdgeList(
                graph, recv_ranges, std::forward<Permutator>(permutator));
        case GraphRepresentation::CSR:
            return ApplyPermutationAndComputeSendBuffersCSR(graph, recv_ranges, std::forward<Permutator>(permutator));
        default:
            throw std::runtime_error("Unexpected graph representation type.");
    }
}

[[maybe_unused]] inline auto ConstructPermutedGraphCSR(
    VertexRange recv_range, const std::vector<SInt>& recv_edges, const std::vector<SSInt>& recv_edge_weights,
    const std::vector<SSInt>& recv_vertex_weights) {
    std::size_t       num_local_vertices = recv_range.second - recv_range.first;
    std::vector<SInt> degree_count(num_local_vertices, 0);

    // Scan received data for degrees. A vertex whose own edges were split across several PEs before the
    // permutation (see Graph::left_partial_vertex/right_partial_vertex) arrives as one message per holder --
    // all of them addressed to this PE, since they all carry the same permuted ID -- so the shares have to be
    // summed here. Assigning would drop every share but the last, and then size adjncy too small for the
    // copies below.
    for (std::size_t cur_pos = 0; cur_pos < recv_edges.size();) {
        const SInt src_id = recv_edges[cur_pos];
        const SInt degree = recv_edges[cur_pos + 1];
        degree_count[src_id - recv_range.first] += degree;
        // skip edges
        cur_pos += 1 + degree + 1;
    }
    XadjArray xadj(num_local_vertices + 1, 0);
    // compute xadj array for received graph
    std::exclusive_scan(degree_count.begin(), degree_count.end(), xadj.begin(), SInt{0});
    if (num_local_vertices > 0) {
        xadj.back() = degree_count.back() + xadj[num_local_vertices - 1];
    }

    // compute adjncy for received graph
    const std::size_t num_local_edges = xadj.back();
    XadjArray         adjncy(num_local_edges);
    // Per-row write cursor rather than xadj[local_src_id] directly, so that the several shares of a formerly
    // split vertex are appended one after the other instead of overwriting each other at the row start. They
    // are appended in source rank order -- the order ExchangeMessageBuffers concatenates messages in -- which
    // is the order the shares appeared in along the graph's edges.
    std::vector<SInt> write_pos(xadj.begin(), xadj.end() - 1);
    for (std::size_t cur_pos = 0; cur_pos < recv_edges.size();) {
        const SInt global_src_id = recv_edges[cur_pos];
        const SInt degree        = recv_edges[cur_pos + 1];
        const SInt local_src_id  = global_src_id - recv_range.first;
        std::copy_n(recv_edges.begin() + cur_pos + 2, degree, adjncy.begin() + write_pos[local_src_id]);
        write_pos[local_src_id] += degree;
        // forward to next received src vertex
        cur_pos += 1 + degree + 1;
    }
    // compute edge weights for received graph
    EdgeWeights edge_weights(recv_edge_weights.empty() ? 0 : num_local_edges);
    // The weight messages carry the same (vertex, share) sequence as the edge messages above and are exchanged
    // the same way, so an identical cursor walk keeps every share's weights aligned with its edges.
    std::vector<SInt> weight_write_pos(xadj.begin(), xadj.end() - 1);
    for (std::size_t cur_pos = 0; cur_pos < recv_edge_weights.size();) {
        const SInt global_src_id = static_cast<SInt>(recv_edge_weights[cur_pos]);
        const SInt degree        = static_cast<SInt>(recv_edge_weights[cur_pos + 1]);
        const SInt local_src_id  = global_src_id - recv_range.first;
        std::copy_n(
            recv_edge_weights.begin() + cur_pos + 2, degree, edge_weights.begin() + weight_write_pos[local_src_id]);
        weight_write_pos[local_src_id] += degree;
        // forward to next received src vertex
        cur_pos += 1 + degree + 1;
    }
    // Sized by the local vertex count, not by the number of received entries: the two differ if some local
    // vertex received no weight (e.g. an empty local range), and the indexed writes below assume the former.
    std::vector<SSInt> vertex_weights(recv_vertex_weights.empty() ? 0 : num_local_vertices, 0);
    for (std::size_t i = 0; i < recv_vertex_weights.size(); i += 2) {
        const auto global_src_id     = recv_vertex_weights[i];
        const auto local_src_id      = global_src_id - recv_range.first;
        vertex_weights[local_src_id] = recv_vertex_weights[i + 1];
    }
    return std::make_tuple(std::move(xadj), std::move(adjncy), std::move(edge_weights), std::move(vertex_weights));
}

[[maybe_unused]] inline auto ConstructPermutedGraphEdgeList(
    VertexRange recv_range, const std::vector<SInt>& recv_edges, const std::vector<SSInt>& recv_edge_weights,
    const std::vector<SSInt>& recv_vertex_weights) {
    std::size_t       num_local_vertices = recv_range.second - recv_range.first;
    std::vector<SInt> degree(num_local_vertices, 0);
    int               rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for (std::size_t i = 0; i < recv_edges.size(); i += 2) {
        const auto src = recv_edges[i];
        ++degree[src - recv_range.first];
    }
    std::vector<SInt> write_idx(num_local_vertices);
    std::exclusive_scan(degree.begin(), degree.end(), write_idx.begin(), SInt{0});
    Edgelist    edgelist(recv_edges.size() / 2); // edges are sent flat - not as pairs
    EdgeWeights edge_weights(recv_edge_weights.size());
    for (std::size_t i = 0; i < recv_edges.size(); i += 2) {
        const auto src          = recv_edges[i];
        const auto dst          = recv_edges[i + 1];
        const auto local_src_id = src - recv_range.first;
        const auto idx          = write_idx[local_src_id];
        edgelist[idx]           = std::make_pair(src, dst);
        if (!edge_weights.empty()) {
            edge_weights[idx] = recv_edge_weights[i / 2];
        }
        ++write_idx[local_src_id];
    }
    std::vector<SSInt> vertex_weights(recv_vertex_weights.empty() ? 0 : num_local_vertices, 0);
    for (std::size_t i = 0; i < recv_vertex_weights.size(); i += 2) {
        const auto global_id     = recv_vertex_weights[i];
        const auto local_id      = global_id - recv_range.first;
        vertex_weights[local_id] = recv_vertex_weights[i + 1];
    }
    return std::make_tuple(std::move(edgelist), std::move(edge_weights), std::move(vertex_weights));
}

// The balance a permuted graph has to be restored to. For generated graphs that is config.redistribution; for
// file graphs the equivalent knob is the input distribution (--distribution), which FileGraphGenerator applies
// in place of config.redistribution. ROOT and EXPLICIT are not balance requests -- they say where the *input*
// should be placed, and permuting necessarily reassigns vertices to PEs anyway -- so they are left alone.
GraphRedistribution EffectiveRedistribution(const PGeneratorConfig& config) {
    if (config.generator == GeneratorType::FILE) {
        switch (config.input_graph.distribution) {
            case GraphDistribution::BALANCE_EDGES:
                return GraphRedistribution::BALANCE_EDGES;
            case GraphDistribution::BALANCE_EDGES_TRUE:
                return GraphRedistribution::BALANCE_EDGES_TRUE;
            default:
                break;
        }
    }
    return config.redistribution;
}

// Relabels every vertex with its permuted ID and then re-establishes the requested edge balance in the new ID
// space. Relabeling itself needs no communication, so the whole operation is a single all-to-all -- the one the
// redistribution primitive performs -- and vertex_range plus the split metadata come back from that primitive
// rather than being recomputed here. For BALANCE_EDGES_TRUE that metadata describes the *new* split vertices:
// permuting does not change the degree sequence, only which IDs the balance boundaries fall on, so a hub that
// was split before is (some other hub) split after.
template <typename Permutator>
void PermuteAndRebalance(
    Graph& graph, const SInt n, const GraphRedistribution redistribution, Permutator&& permute, MPI_Comm comm) {
    const bool csr = graph.representation == GraphRepresentation::CSR;

    Edgelist edges;
    if (csr) {
        // PhysicalVertexRange(), not vertex_range: xadj is indexed by the physically present row space, which
        // includes a left-partial split vertex that vertex_range excludes.
        edges = BuildEdgeListFromCSR(graph.PhysicalVertexRange(), graph.xadj, graph.adjncy);
        graph.FreeCSR();
    } else {
        edges = std::move(graph.edges);
        graph.FreeEdgelist();
    }
    for (auto& [src, dst]: edges) {
        src = permute(src);
        dst = permute(dst);
    }

    // remap_round_robin=false: the IDs have just been permuted, they must not be remapped a second time.
    Edgelist redistributed;
    switch (redistribution) {
        case GraphRedistribution::BALANCE_EDGES:
            graph.vertex_range = RedistributeEdgesBalanced(edges, redistributed, n, /*remap_round_robin=*/false, comm);
            graph.has_split_vertices = false;
            graph.left_partial_vertex.reset();
            graph.right_partial_vertex.reset();
            break;

        case GraphRedistribution::BALANCE_EDGES_TRUE: {
            const EdgeBalancedDistribution distribution =
                RedistributeEdgesTrueBalance(edges, redistributed, n, /*remap_round_robin=*/false, comm);
            graph.vertex_range         = distribution.vertex_range;
            graph.has_split_vertices   = distribution.has_split_vertices;
            graph.left_partial_vertex  = distribution.left_partial_vertex;
            graph.right_partial_vertex = distribution.right_partial_vertex;
            break;
        }

        case GraphRedistribution::BALANCE_VERTICES:
            throw std::runtime_error("PermuteAndRebalance called for a vertex-balanced distribution");
    }

    if (csr) {
        // Same row space as every other split-aware CSR builder (see EdgeListOnlyGenerator::FinalizeCSR): a PE
        // holding a replica of its first vertex needs a row for it, and that vertex sits just below the
        // gap-free vertex_range. Unweighted throughout -- weighted graphs never take this path.
        EdgeWeights no_edge_weights;
        std::tie(graph.xadj, graph.adjncy) =
            BuildCSRFromEdgeList(graph.PhysicalVertexRange(), redistributed, no_edge_weights);
    } else {
        graph.edges = std::move(redistributed);
    }
}
} // namespace

void Generator::PermuteVertices([[maybe_unused]] const PGeneratorConfig& config, [[maybe_unused]] MPI_Comm comm) {
#ifdef KAGEN_XXHASH_FOUND
    int size = -1;
    int rank = -1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // config.n is the generator's own parameter and not every generator sets it: the file generator takes its
    // vertex count from the input file and leaves config.n at 0, which would underflow the permutation domain
    // below. vertex_range is a gap-free partition of [0, n) for both representations, so its largest end is n.
    SInt n = config.n;
    if (n == 0) {
        n = graph_.vertex_range.second;
        MPI_Allreduce(MPI_IN_PLACE, &n, 1, KAGEN_MPI_SINT, MPI_MAX, comm);
    }
    if (n == 0) {
        return; // empty graph -- nothing to permute
    }

    // Coordinates are stored per local vertex and are not carried through the exchange below, so the permuted
    // graph would silently end up with its old PE's coordinates. Rejected outright rather than shipping wrong
    // data. Collective, like the weight check below: a PE with an empty local graph has no coordinates even
    // when the graph as a whole does, and everything past this point is collective.
    bool has_coordinates = !graph_.coordinates.first.empty() || !graph_.coordinates.second.empty();
    MPI_Allreduce(MPI_IN_PLACE, &has_coordinates, 1, MPI_C_BOOL, MPI_LOR, comm);
    if (has_coordinates) {
        throw ConfigurationError(
            "vertex permutation (--permute) is not supported together with coordinate output (--coordinates): "
            "the permutation reassigns vertices to PEs, but coordinates are not moved along with them");
    }

    auto permutator = random_permutation::FeistelPseudoRandomPermutation::buildPermutation(n - 1, 0);
    auto permute    = [&permutator](SInt v) {
        return permutator.f(v);
    };

    // Permuting relabels vertices; it does not preserve how many *edges* a PE holds, since each PE's new vertex
    // set is an essentially random subset of the old one. So for the edge-balanced modes the balance has to be
    // re-established in the permuted ID space -- otherwise --permute silently downgrades
    // --redistribution=balance-edges[-strict] to the plain vertex-balanced distribution computed below, leaving
    // the per-PE edge count to the whims of the permutation (a hub vertex lands, with its entire degree, on
    // whichever PE happens to own its permuted ID).
    const GraphRedistribution redistribution = EffectiveRedistribution(config);
    if (redistribution != GraphRedistribution::BALANCE_VERTICES) {
        // The RedistributeEdges*() primitives carry no weights alongside the edge list they move -- the same
        // restriction FinalizeGraphFragment() documents for edge-balanced input distributions -- and
        // rebalancing changes both which PE holds which edges and which vertices it owns, so weights cannot be
        // carried along unchanged. Checked collectively: a PE with an empty local graph sees empty weight
        // arrays even when the graph as a whole is weighted, and PermuteAndRebalance() below is collective.
        bool has_weights = !graph_.vertex_weights.empty() || !graph_.edge_weights.empty();
        MPI_Allreduce(MPI_IN_PLACE, &has_weights, 1, MPI_C_BOOL, MPI_LOR, comm);
        if (has_weights) {
            std::stringstream msg;
            msg << "vertex permutation (--permute) together with '" << redistribution
                << "' is not supported for weighted graphs: re-establishing the edge balance after permuting "
                   "moves edges and vertices independently of their weights; use "
                   "--redistribution=balance-vertices, --drop-edge-weights and/or --drop-vertex-weights";
            throw ConfigurationError(msg.str());
        }

        PermuteAndRebalance(graph_, n, redistribution, permute, comm);
        return;
    }

    // all PE get n / size vertices
    // the first n modulo size PEs obtain one additional vertices.
    const SInt vertices_per_pe             = n / size;
    const PEID num_pe_with_additional_node = n % size;
    const bool has_pe_additional_node      = rank < num_pe_with_additional_node;
    const SInt begin_vertices              = std::min(num_pe_with_additional_node, rank) + rank * vertices_per_pe;
    const SInt end_vertices                = begin_vertices + vertices_per_pe + has_pe_additional_node;

    VertexRange              recv_range{begin_vertices, end_vertices};
    std::vector<VertexRange> recv_ranges = AllgatherVertexRange(recv_range, comm);

    auto [send_buffers, edge_weight_send_buffers, vertex_weight_send_buffers] =
        ApplyPermutationAndComputeSendBuffers(graph_, recv_ranges, permute);
    auto recv_edges          = ExchangeMessageBuffers(std::move(send_buffers), KAGEN_MPI_SINT, comm);
    auto recv_edge_weights   = ExchangeMessageBuffers(std::move(edge_weight_send_buffers), KAGEN_MPI_SSINT, comm);
    auto recv_vertex_weights = ExchangeMessageBuffers(std::move(vertex_weight_send_buffers), KAGEN_MPI_SSINT, comm);

    switch (desired_representation_) {
        case GraphRepresentation::EDGE_LIST: {
            auto [permuted_edgelist, permuted_edge_weights, permuted_vertex_weights] =
                ConstructPermutedGraphEdgeList(recv_range, recv_edges, recv_edge_weights, recv_vertex_weights);
            graph_.edges          = std::move(permuted_edgelist);
            graph_.edge_weights   = std::move(permuted_edge_weights);
            graph_.vertex_weights = std::move(permuted_vertex_weights);
            break;
        }
        case GraphRepresentation::CSR: {
            auto [permuted_xadj, permuted_adjncy, permuted_edge_weights, permuted_vertex_weights] =
                ConstructPermutedGraphCSR(recv_range, recv_edges, recv_edge_weights, recv_vertex_weights);

            graph_.xadj           = std::move(permuted_xadj);
            graph_.adjncy         = std::move(permuted_adjncy);
            graph_.edge_weights   = std::move(permuted_edge_weights);
            graph_.vertex_weights = std::move(permuted_vertex_weights);
            break;
        }
    }
    SetVertexRange(recv_range);

    // The exchange above routes every edge of vertex v to the single PE owning permute(v), and recv_ranges is a
    // gap-free partition of [0, n), so a vertex whose edges used to live on several PEs is whole again here --
    // the permutation resolves every split by construction. Clear the now-stale metadata: its offsets and
    // counts refer to the pre-exchange local edge list, and a left-over left_partial_vertex would make
    // PhysicalVertexRange() claim a row space one vertex larger than the one that actually exists.
    SetHasSplitVertices(false);
    SetPartialVertices(std::nullopt, std::nullopt);
#endif // KAGEN_XXHASH_FOUND
}

std::unique_ptr<kagen::VertexWeightGenerator>
CreateVertexWeightGenerator(const VertexWeightConfig weight_config, MPI_Comm comm) {
    switch (weight_config.generator_type) {
        case VertexWeightGeneratorType::DEFAULT:
            return std::make_unique<DefaultVertexWeightGenerator>(weight_config);
        case VertexWeightGeneratorType::VOIDING:
            return std::make_unique<VoidingVertexWeightGenerator>(weight_config);
        case VertexWeightGeneratorType::UNIFORM_RANDOM:
            return std::make_unique<UniformRandomVertexWeightGenerator>(weight_config, comm);
    }

    throw std::runtime_error("invalid weight generator type");
}

void Generator::GenerateVertexWeights(VertexWeightConfig weight_config, MPI_Comm comm) {
    std::unique_ptr<kagen::VertexWeightGenerator> vertex_weight_generator =
        CreateVertexWeightGenerator(weight_config, comm);

    switch (desired_representation_) {
        case GraphRepresentation::EDGE_LIST:
            vertex_weight_generator->GenerateVertexWeights(graph_.vertex_range, graph_.edges, graph_.vertex_weights);
            break;
        case GraphRepresentation::CSR:
            vertex_weight_generator->GenerateVertexWeights(
                graph_.vertex_range, graph_.xadj, graph_.adjncy, graph_.vertex_weights);
            break;
    }
}

void Generator::FinalizeEdgeList(MPI_Comm) {}

void Generator::FinalizeCSR(MPI_Comm) {}

void CSROnlyGenerator::GenerateEdgeList() {
    GenerateCSR();
}

void CSROnlyGenerator::FinalizeEdgeList(MPI_Comm comm) {
    if (graph_.xadj.empty()) {
        return;
    }

    // Otherwise, we have generated the graph in CSR representation, but
    // actually want edge list representation -> transform graph
    FinalizeCSR(comm);
    graph_.edges = BuildEdgeListFromCSR(graph_.PhysicalVertexRange(), graph_.xadj, graph_.adjncy);
    {
        XadjArray tmp;
        std::swap(graph_.xadj, tmp);
    }
    {
        AdjncyArray tmp;
        std::swap(graph_.adjncy, tmp);
    }
}

void EdgeListOnlyGenerator::GenerateCSR() {
    GenerateEdgeList();
}

void EdgeListOnlyGenerator::FinalizeCSR(MPI_Comm comm) {
    if (!graph_.xadj.empty()) {
        return;
    }

    // Otherwise, we have generated the graph in edge list representation, but
    // actually want CSR format --> transform graph
    FinalizeEdgeList(comm);

    // BuildCSRFromEdgeList gives every vertex in the given range a row (indexing by `from - range.first`),
    // isolated vertices included as empty rows. Building from vertex_range (the gap-free ownership range) would
    // underflow on a PE holding a *replica* of its first vertex (left_partial_vertex set): that vertex is
    // credited to the lower-rank canonical PE and so lies just below vertex_range, yet its edges are physically
    // here. PhysicalVertexRange() is exactly vertex_range extended down by one to cover that replica row, giving
    // the same physically-present row-space layout the strict CSR file reader produces (see
    // FinalizeGraphFragment). graph_.vertex_range itself is left untouched -- it keeps meaning the gap-free
    // ownership range, same as for the edge-list representation.
    const VertexRange csr_range          = graph_.PhysicalVertexRange();
    std::tie(graph_.xadj, graph_.adjncy) = BuildCSRFromEdgeList(csr_range, graph_.edges, graph_.edge_weights);
    {
        Edgelist tmp;
        std::swap(graph_.edges, tmp);
    }
}

SInt Generator::GetNumberOfEdges() const {
    return std::max(graph_.adjncy.size(), graph_.edges.size());
}

bool Generator::HasSplitVertices() const {
    return graph_.has_split_vertices;
}

Graph Generator::Take() {
    return std::move(graph_);
}

Edgelist Generator::TakeNonlocalEdges() {
    return std::move(nonlocal_edges_);
}

void Generator::SetVertexRange(const VertexRange vertex_range) {
    graph_.vertex_range = vertex_range;
}

void Generator::SetHasSplitVertices(const bool has_split_vertices) {
    graph_.has_split_vertices = has_split_vertices;
}

void Generator::SetPartialVertices(
    std::optional<SplitVertexInfo> left_partial_vertex, std::optional<SplitVertexInfo> right_partial_vertex) {
    graph_.left_partial_vertex  = left_partial_vertex;
    graph_.right_partial_vertex = right_partial_vertex;
}

void Generator::FilterDuplicateEdges() {
    std::sort(graph_.edges.begin(), graph_.edges.end());
    auto it = std::unique(graph_.edges.begin(), graph_.edges.end());
    graph_.edges.erase(it, graph_.edges.end());
}

void Generator::Reset() {
    graph_.Clear();
}

GeneratorFactory::~GeneratorFactory() = default;

PGeneratorConfig
GeneratorFactory::NormalizeParameters(PGeneratorConfig config, PEID, const PEID size, const bool output) const {
    if (config.k == 0) {
        config.k = static_cast<SInt>(size);
        if (output) {
            std::cout << "Setting number of chunks to " << config.k << std::endl;
        }
    }
    return config;
}

namespace {
bool IsPowerOfTwo(const SInt value) {
    return (value & (value - 1)) == 0;
}

bool IsSquare(const SInt value) {
    const SInt root = std::round(std::sqrt(value));
    return root * root == value;
}

bool IsCubic(const SInt value) {
    const SInt root = std::round(std::cbrt(value));
    return root * root * root == value;
}
} // namespace

void GeneratorFactory::EnsureSquarePowerOfTwoChunkSize(
    PGeneratorConfig& config, const PEID size, const bool output) const {
    if (config.k == 0) {
        if (IsSquare(size) && IsPowerOfTwo(size)) {
            config.k = static_cast<SInt>(size);
        } else {
            const SInt l = std::ceil(std::log2(size));
            config.k     = 1 << l;
            if (!IsSquare(config.k)) {
                config.k *= 2;
            }
            while (std::ceil(1.0 * config.k / size) > (1.0 + config.max_vertex_imbalance) * config.k / size) {
                config.k <<= 2;
            }
        }
        if (output) {
            std::cout << "Setting number of chunks to " << config.k << std::endl;
        }
    } else if (config.k < static_cast<SInt>(size) || !IsSquare(config.k) || !IsPowerOfTwo(config.k)) {
        throw ConfigurationError("number of chunks must be square power of two and larger than number of PEs");
    }
}

void GeneratorFactory::EnsureCubicPowerOfTwoChunkSize(
    PGeneratorConfig& config, const PEID size, const bool output) const {
    if (config.k == 0) {
        if (IsCubic(size) && IsPowerOfTwo(size)) {
            config.k = static_cast<SInt>(size);
        } else {
            const SInt l = std::ceil(std::log2(size));
            config.k     = 1 << l;
            if (!IsCubic(config.k)) {
                config.k *= 2;
            }
            if (!IsCubic(config.k)) {
                config.k *= 2;
            }

            while (std::ceil(1.0 * config.k / size) > (1.0 + config.max_vertex_imbalance) * config.k / size) {
                config.k <<= 3;
            }
        }
        if (output) {
            std::cout << "Setting number of chunks to " << config.k << std::endl;
        }
    } else if (config.k < static_cast<SInt>(size) || !IsCubic(config.k)) {
        throw ConfigurationError("number of chunks must be cubic and larger than the number of PEs");
    }
}

void GeneratorFactory::EnsureOneChunkPerPE(PGeneratorConfig& config, const PEID size) const {
    if (config.k != static_cast<SInt>(size)) {
        throw ConfigurationError("number of chunks must match the number of PEs");
    }
}
} // namespace kagen
