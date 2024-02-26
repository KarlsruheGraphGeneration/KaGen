#include "kagen/io.h"

#include "kagen/context.h"
#include "kagen/io/coordinates.h"
#include "kagen/io/dot.h"
#include "kagen/io/edgelist.h"
#include "kagen/io/freight-netl.h"
#include "kagen/io/hmetis.h"
#include "kagen/io/metis.h"
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
        factories[FileFormat::FREIGHT_NETL_EP] = std::make_unique<FreightNetlEpFactory>();
        factories[FileFormat::FREIGHT_NETL]    = std::make_unique<FreightNetlFactory>();
        factories[FileFormat::HMETIS_EP]       = std::make_unique<HmetisEpFactory>();
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

GraphFragment ReadGraphFragment(
    GraphReader& reader, const GraphRepresentation representation, const InputGraphConfig& config, const PEID rank,
    const PEID size) {
    const auto [n, m] = [&] {
        if (reader.HasDeficit(ReaderDeficits::UNKNOWN_NUM_VERTICES)
            && reader.HasDeficit(ReaderDeficits::UNKNOWN_NUM_EDGES)) {
            return std::pair<SInt, SInt>{size, size};
        }

        auto [n, m] = reader.ReadSize();
        if (reader.HasDeficit(ReaderDeficits::UNKNOWN_NUM_VERTICES)) {
            n = size;
        }
        if (reader.HasDeficit(ReaderDeficits::UNKNOWN_NUM_EDGES)) {
            m = size;
        }
        return std::pair<SInt, SInt>{n, m};
    }();

    if (reader.HasDeficit(ReaderDeficits::REQUIRES_REDISTRIBUTION)
        && config.distribution == GraphDistribution::BALANCE_EDGES) {
        throw std::invalid_argument("not implemented");
    }

    // If we need postprocessing, always generate an edge list because postprocessing is not implemented for CSR
    GraphRepresentation actual_representation =
        reader.HasDeficit(ReaderDeficits::REQUIRES_REDISTRIBUTION) ? GraphRepresentation::EDGE_LIST : representation;

    SInt from    = 0;
    SInt to_node = std::numeric_limits<SInt>::max();
    SInt to_edge = std::numeric_limits<SInt>::max();

    switch (config.distribution) {
        case GraphDistribution::ROOT:
            if (rank == 0) {
                from    = 0;
                to_node = n;
            } else {
                from    = n;
                to_node = n;
            }
            break;
        case GraphDistribution::BALANCE_VERTICES:
            std::tie(from, to_node) = ComputeRange(n, size, rank);
            break;

        case GraphDistribution::BALANCE_EDGES: {
            const auto edge_range = ComputeRange(m, size, rank);
            from                  = reader.FindNodeByEdge(edge_range.first);
            to_edge               = edge_range.second;
            break;
        }
    }

    return {
        reader.Read(from, to_node, to_edge, actual_representation),
        reader.Deficits(),
    };
}

Graph FinalizeGraphFragment(GraphFragment fragment, const bool output, MPI_Comm comm) {
    if (fragment.deficits & ReaderDeficits::REQUIRES_REDISTRIBUTION) {
        if (fragment.graph.representation == GraphRepresentation::CSR) {
            throw std::invalid_argument("not implemented");
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

        std::tie(fragment.graph.vertex_range.first, fragment.graph.vertex_range.second) = ComputeRange(n, size, rank);
        RedistributeEdgesByVertexRange(fragment.graph.edges, fragment.graph.vertex_range, comm);
    }

    return std::move(fragment.graph);
}

void WriteGraph(GraphWriter& writer, const OutputGraphConfig& config, const bool output, MPI_Comm comm) {
    const PEID size = GetCommSize(comm);
    const PEID rank = GetCommRank(comm);

    const std::string filename = config.distributed ? config.filename + "." + std::to_string(rank) : config.filename;

    // Overwrite file if it already exists
    { std::ofstream out(filename); }

    if (config.distributed) {
        // Distributed output: each PE writes its part of the graph to its own file
        // This allows parallel writes to parallel file systems

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
        // Sequential output (default): all PEs write the the same file, sequentially

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
