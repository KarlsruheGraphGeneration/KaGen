#include "kagen/io.h"

#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/io/coordinates.h"
#include "kagen/io/dot.h"
#include "kagen/io/edgelist.h"
#include "kagen/io/freight-netl.h"
#include "kagen/io/hmetis.h"
#include "kagen/io/metis.h"
#include "kagen/io/netd_are.h"
#include "kagen/io/parhip.h"
#include "kagen/kagen.h"
#include "kagen/tools/postprocessor.h"
#include "kagen/tools/utils.h"

#include <mpi.h>

#include <algorithm>
#include <memory>
#include <sstream>

namespace kagen {
const std::unordered_map<FileFormat, std::unique_ptr<FileFormatFactory>>& GetGraphFormatFactories() {
    static std::unordered_map<FileFormat, std::unique_ptr<FileFormatFactory>> factories;
    if (factories.empty()) {
        factories[FileFormat::NOOP]                        = std::make_unique<NoopFactory>();
        factories[FileFormat::EDGE_LIST]                   = std::make_unique<EdgelistFactory>();
        factories[FileFormat::EDGE_LIST_UNDIRECTED]        = std::make_unique<UndirectedEdgelistFactory>();
        factories[FileFormat::BINARY_EDGE_LIST]            = std::make_unique<BinaryEdgelistFactory>();
        factories[FileFormat::BINARY_EDGE_LIST_UNDIRECTED] = std::make_unique<UndirectedBinaryEdgelistFactory>();
        factories[FileFormat::PLAIN_EDGE_LIST]             = std::make_unique<PlainEdgelistFactory>();
        factories[FileFormat::XTRAPULP]                    = std::make_unique<XtrapulpFactory>();
        factories[FileFormat::METIS]                       = std::make_unique<MetisFactory>();
        factories[FileFormat::HMETIS]                      = std::make_unique<HmetisFactory>();
        factories[FileFormat::HMETIS_DIRECTED]             = std::make_unique<DirectedHmetisFactory>();
        factories[FileFormat::DOT]                         = std::make_unique<DotFactory>();
        factories[FileFormat::DOT_DIRECTED]                = std::make_unique<DirectedDotFactory>();
        factories[FileFormat::COORDINATES]                 = std::make_unique<CoordinatesFactory>();
        factories[FileFormat::PARHIP]                      = std::make_unique<ParhipFactory>();

        // Experimental formats
        factories[FileFormat::FREIGHT_NETL_EP]           = std::make_unique<FreightNetlEpFactory>();
        factories[FileFormat::FREIGHT_NETL]              = std::make_unique<FreightNetlFactory>();
        factories[FileFormat::HMETIS_EP]                 = std::make_unique<HmetisEpFactory>();
        factories[FileFormat::WEIGHTED_BINARY_EDGE_LIST] = std::make_unique<WeightedBinaryEdgelistFactory>();
        factories[FileFormat::NETD_ARE]                  = std::make_unique<NetDAreFactory>();
    }
    return factories;
}

const std::unique_ptr<FileFormatFactory>& GetGraphFormatFactory(const FileFormat format) {
    const auto& factories = GetGraphFormatFactories();
    const auto  it        = factories.find(format);
    if (it != factories.end()) {
        return (*it).second;
    }

    std::stringstream error_msg;
    error_msg << "there is no file format with name " << format;
    throw IOError(error_msg.str());
}

namespace {
std::string GetExtension(const std::string& filename) {
    const auto last_dot_pos = filename.find_last_of('.');
    if (last_dot_pos != std::string::npos) {
        return filename.substr(last_dot_pos + 1);
    }
    return filename;
}
} // namespace

std::unique_ptr<GraphReader>
CreateGraphReader(const std::string& filename, const InputGraphConfig& config, const PEID rank, const PEID size) {
    const std::string extension = GetExtension(filename);
    const auto&       factories = GetGraphFormatFactories();

    // Each file format can register itself under multiple extensions, which are sorted by priority.
    // Thus, after each loop over all formats, increase priority and try again, until there are no more candidates left.
    std::size_t candidates = factories.size();
    for (std::size_t priority = 0; candidates > 0; ++priority) {
        candidates = 0;

        for (const auto& [format, factory]: factories) {
            const auto& extensions = factory->DefaultExtensions();
            if (extensions.size() <= priority) {
                continue;
            }

            ++candidates;

            if (extensions[priority] == extension) {
                auto reader = factory->CreateReader(config, rank, size);
                if (reader != nullptr) {
                    return reader;
                }
            }
        }
    }

    std::stringstream error_msg;
    error_msg << "no file format found for filename " << filename << " with file extension " << extension;
    throw IOError(error_msg.str());
}

std::unique_ptr<GraphReader>
CreateGraphReader(const FileFormat format, const InputGraphConfig& config, const PEID rank, const PEID size) {
    if (format == FileFormat::EXTENSION) {
        return CreateGraphReader(config.filename, config, rank, size);
    }

    const auto& factories = GetGraphFormatFactories();
    const auto  it        = factories.find(format);
    if (it != factories.end()) {
        auto reader = (*it).second->CreateReader(config, rank, size);
        if (reader != nullptr) {
            return reader;
        }
    }

    std::stringstream error_msg;
    error_msg << "file format " << format << " not available for reading";
    throw IOError(error_msg.str());
}

namespace {

std::vector<SInt> ReadExplicitVertexDistribution(const std::string& filename, const bool is_prefix_sum) {
    MappedFileToker toker(filename);
    toker.SkipSpaces();

    std::vector<SInt> number_of_vertices;
    while (toker.ValidPosition()) {
        toker.SkipSpaces();
        number_of_vertices.push_back(toker.ScanUnsigned());
        toker.SkipLine();
    }

    if (!is_prefix_sum) {
        number_of_vertices.push_back(0);
        std::exclusive_scan(number_of_vertices.begin(), number_of_vertices.end(), number_of_vertices.begin(), SInt{0});
    }

    return number_of_vertices;
}

// Collectively checks that the union of all PEs' local edge slices is globally sorted by tail vertex: each PE's
// edges are non-decreasing by tail, and the last tail of the nearest non-empty PE to the left is <= this PE's
// first tail (equality is allowed -- that is a split boundary vertex). Only then does a contiguous global edge
// range map to a coherent, splittable vertex partition, which the direct strict-read path requires for
// edge-list readers (whose on-disk order is otherwise arbitrary). Empty PEs contribute the sentinel n.
bool CheckGloballyTailSorted(const Edgelist& edges, const SInt n, MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const bool local_sorted =
        std::is_sorted(edges.begin(), edges.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

    const bool has_local  = !edges.empty();
    const SInt first_tail = has_local ? edges.front().first : n;
    const SInt last_tail  = has_local ? edges.back().first : n;

    std::vector<SInt> all_last_tails(static_cast<std::size_t>(size));
    MPI_Allgather(&last_tail, 1, KAGEN_MPI_SINT, all_last_tails.data(), 1, KAGEN_MPI_SINT, comm);

    bool boundary_ok = true;
    if (has_local) {
        for (PEID pe = rank - 1; pe >= 0; --pe) {
            if (all_last_tails[static_cast<std::size_t>(pe)] != n) { // nearest non-empty predecessor
                boundary_ok = all_last_tails[static_cast<std::size_t>(pe)] <= first_tail;
                break;
            }
        }
    }

    bool ok = local_sorted && boundary_ok;
    MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_C_BOOL, MPI_LAND, comm);
    return ok;
}

// Heals a directly-read strict edge slice (edge list, globally tail-sorted) to whole-vertex ("vertex-atomic")
// ownership for --distribution=balance-edges: each PE hands its leading partial boundary vertex's edges to its
// left neighbor (the lower-rank PE, which owns the whole vertex) and receives its right neighbor's, via a single
// neighbor exchange. Edge weights travel with their edges. Returns false -- signalling the caller to fall back to
// RedistributeEdgesBalanced -- when a single exchange cannot resolve ownership: any empty PE, or a vertex whose
// degree spans more than two PEs.
bool HealToVertexAtomic(Graph& graph, const SInt n, MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    Edgelist&    edges      = graph.edges;
    EdgeWeights& weights    = graph.edge_weights;
    const bool   weighted   = !weights.empty();
    const bool   has_local  = !edges.empty();
    const SInt   first_tail = has_local ? edges.front().first : n;
    const SInt   last_tail  = has_local ? edges.back().first : n;

    const BoundaryOwnership bo = ComputeBoundaryOwnership(first_tail, last_tail, has_local, n, comm);

    bool needs_fallback = !has_local || (bo.is_left_partial && bo.is_right_partial && first_tail == last_tail);
    MPI_Allreduce(MPI_IN_PLACE, &needs_fallback, 1, MPI_C_BOOL, MPI_LOR, comm);
    if (needs_fallback) {
        return false;
    }

    // Split off the leading run (all edges whose tail == first_tail) to send left, if this PE shares first_tail
    // with its left neighbor.
    std::vector<SInt>  send_edges;
    std::vector<SSInt> send_weights;
    SInt               keep_from = 0;
    if (bo.is_left_partial) {
        auto tail_greater = [](SInt v, const std::pair<SInt, SInt>& e) {
            return v < e.first;
        };
        const auto split_it = std::upper_bound(edges.begin(), edges.end(), first_tail, tail_greater);
        keep_from           = static_cast<SInt>(std::distance(edges.begin(), split_it));
        send_edges.reserve(static_cast<std::size_t>(keep_from) * 2);
        for (SInt i = 0; i < keep_from; ++i) {
            send_edges.push_back(edges[i].first);
            send_edges.push_back(edges[i].second);
            if (weighted) {
                send_weights.push_back(weights[i]);
            }
        }
    }

    const PEID left_dest = bo.is_left_partial ? rank - 1 : MPI_PROC_NULL;
    const PEID right_src = bo.is_right_partial ? rank + 1 : MPI_PROC_NULL;

    SInt send_count = static_cast<SInt>(send_edges.size());
    SInt recv_count = 0;
    MPI_Sendrecv(
        &send_count, 1, KAGEN_MPI_SINT, left_dest, 0, &recv_count, 1, KAGEN_MPI_SINT, right_src, 0, comm,
        MPI_STATUS_IGNORE);
    std::vector<SInt> recv_edges(static_cast<std::size_t>(recv_count));
    MPI_Sendrecv(
        send_edges.data(), send_count, KAGEN_MPI_SINT, left_dest, 1, recv_edges.data(), recv_count, KAGEN_MPI_SINT,
        right_src, 1, comm, MPI_STATUS_IGNORE);

    std::vector<SSInt> recv_weights;
    if (weighted) {
        SInt send_w = static_cast<SInt>(send_weights.size());
        SInt recv_w = 0;
        MPI_Sendrecv(
            &send_w, 1, KAGEN_MPI_SINT, left_dest, 2, &recv_w, 1, KAGEN_MPI_SINT, right_src, 2, comm,
            MPI_STATUS_IGNORE);
        recv_weights.resize(static_cast<std::size_t>(recv_w));
        MPI_Sendrecv(
            send_weights.data(), send_w, KAGEN_MPI_SSINT, left_dest, 3, recv_weights.data(), recv_w, KAGEN_MPI_SSINT,
            right_src, 3, comm, MPI_STATUS_IGNORE);
    }

    // Rebuild: drop the sent leading run, keep the rest, append the received run (all for last_tail, so tail
    // order is preserved).
    Edgelist    new_edges;
    EdgeWeights new_weights;
    new_edges.reserve((edges.size() - static_cast<std::size_t>(keep_from)) + recv_edges.size() / 2);
    for (SInt i = keep_from; i < static_cast<SInt>(edges.size()); ++i) {
        new_edges.push_back(edges[i]);
        if (weighted) {
            new_weights.push_back(weights[i]);
        }
    }
    for (std::size_t i = 0; i < recv_edges.size(); i += 2) {
        new_edges.emplace_back(recv_edges[i], recv_edges[i + 1]);
    }
    if (weighted) {
        new_weights.insert(new_weights.end(), recv_weights.begin(), recv_weights.end());
    }

    edges   = std::move(new_edges);
    weights = std::move(new_weights);
    return true;
}

} // namespace

GraphFragment ReadGraphFragment(
    GraphReader& reader, const GraphRepresentation representation, const InputGraphConfig& config, const PEID rank,
    const PEID size) {
    const auto [n, m] = reader.ReadSize();

    // Direct edge-balanced read: each PE reads its own contiguous edge slice [from_edge, to_edge) straight from
    // disk -- the same closed-form partition the redistribution primitives target -- without the vertex-balanced
    // read + all-to-all. No edge movement, no second in-memory copy of the input.
    //  - BALANCE_EDGES_TRUE takes it whenever the reader supports a strict edge-range read (ParHIP/METIS always;
    //    a globally tail-sorted edge list after a cheap check).
    //  - BALANCE_EDGES takes it only for readers that would otherwise redistribute (the sortedness-check readers,
    //    i.e. edge lists); ParHIP/METIS already serve BALANCE_EDGES directly via FindNodeByEdge below. The strict
    //    slice is then healed to whole-vertex ownership in FinalizeGraphFragment.
    const bool direct_strict = config.distribution == GraphDistribution::BALANCE_EDGES_TRUE;
    const bool direct_balance_edges =
        config.distribution == GraphDistribution::BALANCE_EDGES && reader.StrictEdgeRangeRequiresSortednessCheck();
    if ((direct_strict || direct_balance_edges) && reader.CanReadStrictEdgeRange()
        && !reader.HasDeficit(ReaderDeficits::UNKNOWN_NUM_EDGES)) {
        const auto [from_edge, to_edge] = ComputeRange(m, size, rank);

        GraphFragment fragment;
        fragment.graph                   = reader.ReadStrictEdgeRange(from_edge, to_edge, representation);
        fragment.deficits                = reader.Deficits();
        fragment.strict_edge_direct      = true;
        fragment.strict_needs_sort_check = reader.StrictEdgeRangeRequiresSortednessCheck();
        fragment.heal_to_vertex_atomic   = direct_balance_edges;
        fragment.n                       = reader.HasDeficit(ReaderDeficits::UNKNOWN_NUM_VERTICES) ? 0 : n;
        return fragment;
    }

    const bool reader_requires_redistribution = reader.HasDeficit(ReaderDeficits::REQUIRES_REDISTRIBUTION);
    const bool needs_postprocessing =
        reader_requires_redistribution || config.distribution == GraphDistribution::BALANCE_EDGES_TRUE;

    // If we need postprocessing, always generate an edge list because postprocessing is not implemented for CSR
    GraphRepresentation actual_representation = needs_postprocessing ? GraphRepresentation::EDGE_LIST : representation;

    SInt from    = 0;
    SInt to_node = std::numeric_limits<SInt>::max();
    SInt to_edge = std::numeric_limits<SInt>::max();

    switch (config.distribution) {
        case GraphDistribution::ROOT: {
            if (rank == 0) {
                from    = 0;
                to_node = n;
            } else {
                from    = n;
                to_node = n;
            }
            break;
        }

        case GraphDistribution::BALANCE_VERTICES:
        case GraphDistribution::BALANCE_EDGES_TRUE: {
            // True edge-balancing needs the real (deduplicated, exactly-balanced) edge positions, which
            // FindNodeByEdge below cannot provide -- it only supports the vertex-atomic approximation used by
            // BALANCE_EDGES. Read an ordinary vertex-balanced partition here (works for every reader, including
            // METIS/ParHIP, since it's just a plain vertex-range read) and let FinalizeGraphFragment's
            // postprocessing pass do the real redistribution.
            std::tie(from, to_node) = ComputeRange(n, size, rank);
            break;
        }

        case GraphDistribution::BALANCE_EDGES: {
            if (reader_requires_redistribution) {
                // FindNodeByEdge is not meaningful for these readers (e.g. plain edgelist, where it's an
                // unimplemented stub); read an arbitrary partition instead and let FinalizeGraphFragment's
                // postprocessing pass do the real edge-balanced redistribution.
                std::tie(from, to_node) = ComputeRange(n, size, rank);
            } else {
                const auto edge_range = ComputeRange(m, size, rank);
                from                  = reader.FindNodeByEdge(edge_range.first);
                to_edge               = edge_range.second;
            }
            break;
        }

        case GraphDistribution::EXPLICIT: {
            auto distribution = ReadExplicitVertexDistribution(
                config.explicit_distribution_filename, config.explicit_distribution_is_prefix_sum);
            from    = distribution[rank];
            to_node = distribution[rank + 1];
            break;
        }
    }

    return {
        reader.Read(from, to_node, to_edge, actual_representation),
        reader.Deficits(),
    };
}

Graph FinalizeGraphFragment(GraphFragment fragment, const InputGraphConfig& config, const bool output, MPI_Comm comm) {
    // Direct strict edge-balanced read: the fragment already holds this PE's exact edge slice. All that remains
    // is to compute the gap-free vertex_range and split-vertex metadata from one MPI_Allgather of the per-PE
    // boundary tails -- no edge movement, no in-memory duplication. Edge-list inputs are first checked for
    // global tail-sortedness; if that fails, fall back to the redistribution path (reusing the edges already in
    // memory, no re-read).
    if (fragment.strict_edge_direct) {
        Graph& graph = fragment.graph;
        // Trust the reader's exact vertex count when it has one (even if it is 0, i.e. an empty graph); only
        // edge-list readers (UNKNOWN_NUM_VERTICES) must derive it from the gathered edges.
        const SInt n = (fragment.deficits & ReaderDeficits::UNKNOWN_NUM_VERTICES)
                           ? FindNumberOfVerticesInEdgelist(graph.edges, comm)
                           : fragment.n;

        const bool sorted_ok = !fragment.strict_needs_sort_check || CheckGloballyTailSorted(graph.edges, n, comm);

        // --distribution=balance-edges: heal the strict slice to whole-vertex ownership. If that is not possible
        // (empty PE or a vertex spanning >2 PEs), fall back to the redistribution path below.
        if (sorted_ok && fragment.heal_to_vertex_atomic && !HealToVertexAtomic(graph, n, comm)) {
            fragment.strict_edge_direct = false;
        } else if (sorted_ok) {
            if (output) {
                std::cout << "reading edge-balanced ... " << std::flush;
            }

            if (graph.representation == GraphRepresentation::EDGE_LIST) {
                const EdgeBalancedDistribution dist = ComputeEdgeBalancedBoundaries(graph.edges, n, comm);
                graph.vertex_range                  = dist.vertex_range;
                graph.has_split_vertices            = dist.has_split_vertices;
                graph.left_partial_vertex           = dist.left_partial_vertex;
                graph.right_partial_vertex          = dist.right_partial_vertex;
            } else {
                // CSR: the reader gives us the physically-present row space [first_tail, last_tail + 1] that xadj
                // indexes -- a split boundary vertex genuinely has a (partial) row on both PEs sharing it, so
                // this row space can overlap the neighbors' by one vertex at each split. Use it (via first/last
                // tail) to compute left_partial_vertex/right_partial_vertex/has_split_vertices below, then
                // replace graph.vertex_range with the gap-free ownership range (same meaning as for EDGE_LIST);
                // PhysicalVertexRange() reconstructs the row space from vertex_range + left_partial_vertex for
                // anything that needs to index xadj/adjncy by global vertex id. The first and last rows are
                // guaranteed non-empty, so their tails are the boundary tails.
                const bool has_local  = !graph.adjncy.empty();
                const SInt first_tail = has_local ? graph.vertex_range.first : n;
                const SInt last_tail  = has_local ? graph.vertex_range.second - 1 : n;

                const BoundaryOwnership bo = ComputeBoundaryOwnership(first_tail, last_tail, has_local, n, comm);
                if (has_local) {
                    const SInt num_rows = graph.vertex_range.second - graph.vertex_range.first;
                    if (bo.is_left_partial) {
                        graph.left_partial_vertex = SplitVertexInfo{first_tail, 0, graph.xadj[1] - graph.xadj[0]};
                    }
                    if (bo.is_right_partial) {
                        const SInt offset          = graph.xadj[num_rows - 1];
                        graph.right_partial_vertex = SplitVertexInfo{last_tail, offset, graph.xadj[num_rows] - offset};
                    }
                }
                graph.has_split_vertices = bo.has_split_vertices;
                graph.vertex_range       = bo.vertex_range;
            }

            return std::move(fragment.graph);
        }

        // Not globally tail-sorted: fall through to the redistribution path with the edges already in memory.
        fragment.strict_edge_direct = false;
    }

    const bool needs_postprocessing = (fragment.deficits & ReaderDeficits::REQUIRES_REDISTRIBUTION)
                                      || config.distribution == GraphDistribution::BALANCE_EDGES_TRUE;
    if (needs_postprocessing) {
        if (fragment.graph.representation == GraphRepresentation::CSR) {
            throw std::invalid_argument("not implemented");
        }

        if (config.distribution == GraphDistribution::BALANCE_EDGES
            || config.distribution == GraphDistribution::BALANCE_EDGES_TRUE) {
            // A PE with an empty local partition (possible for a small graph split across many PEs) would see
            // empty weight arrays locally even though the graph as a whole is weighted, so this check must be
            // collective -- otherwise some PEs would throw while others proceed into the (also collective)
            // redistribution below, deadlocking.
            bool has_weights = !fragment.graph.edge_weights.empty() || !fragment.graph.vertex_weights.empty();
            MPI_Allreduce(MPI_IN_PLACE, &has_weights, 1, MPI_C_BOOL, MPI_LOR, comm);
            if (has_weights) {
                // None of the redistribution primitives below carry edge_weights/vertex_weights alongside the
                // Edgelist they redistribute -- doing so would desync the weights from their edges/vertices
                // rather than just reordering them, since redistribution changes both which PE holds which
                // edges and which vertex range each PE owns.
                throw ConfigurationError(
                    "edge-balanced redistribution (--distribution=balance-edges or balance-edges-strict) is not "
                    "supported together with vertex- or edge-weighted input; use --distribution=balance-vertices, "
                    "--drop-edge-weights, and/or --drop-vertex-weights");
            }
        }

        const PEID size = GetCommSize(comm);
        const PEID rank = GetCommRank(comm);

        if (output) {
            std::cout << "redistributing edges ... " << std::flush;
        }

        const SInt n = [&] {
            SInt n = 0;
            if (fragment.deficits & ReaderDeficits::UNKNOWN_NUM_VERTICES) {
                n = FindNumberOfVerticesInEdgelist(fragment.graph.edges, comm);
            } else {
                n = fragment.graph.vertex_range.second;
                MPI_Bcast(&n, 1, KAGEN_MPI_SINT, size - 1, comm);
            }
            return n;
        }();

        // remap_round_robin=false throughout: file-graph vertex IDs are the user's own on-disk IDs, so they
        // must not be remapped.
        switch (config.distribution) {
            case GraphDistribution::ROOT:
            case GraphDistribution::EXPLICIT:
            case GraphDistribution::BALANCE_VERTICES: {
                std::tie(fragment.graph.vertex_range.first, fragment.graph.vertex_range.second) =
                    ComputeRange(n, size, rank);
                RedistributeEdgesByVertexRange(fragment.graph.edges, fragment.graph.vertex_range, comm);
                break;
            }
            case GraphDistribution::BALANCE_EDGES: {
                Edgelist source = std::move(fragment.graph.edges);
                fragment.graph.vertex_range =
                    RedistributeEdgesBalanced(source, fragment.graph.edges, n, /*remap_round_robin=*/false, comm);
                break;
            }
            case GraphDistribution::BALANCE_EDGES_TRUE: {
                Edgelist                       source = std::move(fragment.graph.edges);
                const EdgeBalancedDistribution distribution =
                    RedistributeEdgesTrueBalance(source, fragment.graph.edges, n, /*remap_round_robin=*/false, comm);
                fragment.graph.vertex_range         = distribution.vertex_range;
                fragment.graph.has_split_vertices   = distribution.has_split_vertices;
                fragment.graph.left_partial_vertex  = distribution.left_partial_vertex;
                fragment.graph.right_partial_vertex = distribution.right_partial_vertex;
                break;
            }
        }
    }

    return std::move(fragment.graph);
}

void WriteGraph(GraphWriter& writer, const OutputGraphConfig& config, const bool output, MPI_Comm comm) {
    const PEID size = GetCommSize(comm);
    const PEID rank = GetCommRank(comm);

    const std::string filename = config.distributed ? config.filename + "." + std::to_string(rank) : config.filename;

    // Overwrite file if it already exists
    if (rank == 0) {
        std::ofstream out(filename);
    }

    if (config.distributed) {
        // Distributed output: each PE writes its part of the graph to
        // its own file This allows parallel writes to parallel file
        // systems

        if (output) {
            std::cout << "Writing graph to [" << filename << ".0";
            if (size > 2) {
                std::cout << ", ...";
            }
            if (size > 1) {
                std::cout << ", " << filename << "." << size - 1;
            }
            std::cout << "] ... " << std::flush;
        }

        bool continue_with_next_pass = true;
        for (int pass = 0; continue_with_next_pass; ++pass) {
            continue_with_next_pass = writer.Write(pass, filename);
        }

        if (output) {
            std::cout << "OK" << std::endl;
        }
    } else {
        // Sequential output (default): all PEs write the the same file,
        // sequentially

        if (output) {
            std::cout << "Writing graph to " << filename << " ..." << std::endl;
        }

        bool continue_with_next_pass = true;
        for (int pass = 0; continue_with_next_pass; ++pass) {
            for (PEID pe = 0; pe < size; ++pe) {
                if (output) {
                    std::cout << "  Writing subgraph of PE " << pe + 1 << " / " << size << " (pass " << pass << ") ... "
                              << std::flush;
                }
                if (rank == pe) {
                    continue_with_next_pass = writer.Write(pass, filename);
                }
                MPI_Barrier(comm);
                if (output) {
                    std::cout << "OK" << std::endl;
                }
            }
        }
    }
}
} // namespace kagen
