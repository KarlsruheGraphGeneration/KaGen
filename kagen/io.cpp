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

} // namespace

GraphFragment ReadGraphFragment(
    GraphReader& reader, const GraphRepresentation representation, const InputGraphConfig& config, const PEID rank,
    const PEID size) {
    const auto [n, m] = reader.ReadSize();

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
                fragment.graph.vertex_range       = distribution.vertex_range;
                fragment.graph.has_split_vertices = distribution.has_split_vertices;
                fragment.graph.partial_vertices   = distribution.partial_vertices;
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
