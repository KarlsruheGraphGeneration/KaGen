#include "kagen/generators/generator.h"
#include "kagen/context.h"
#include "kagen/definitions.h"
#include "kagen/tools/converter.h"
#include "kagen/generators/edgeweights/edge_weight_generator.h"
#include "kagen/generators/edgeweights/hashed/hashed.h"
#include "kagen/generators/edgeweights/random/random.h"
#include "kagen/generators/edgeweights/constant/constant.h"
#include "kagen/generators/edgeweights/distance/distance.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>

namespace kagen {
    Generator::~Generator() = default;

    void Generator::Generate(const GraphRepresentation representation) {
        Reset();

        representation_ = representation;

        switch (representation_) {
            case GraphRepresentation::EDGE_LIST:
                GenerateEdgeList();
                break;

            case GraphRepresentation::CSR:
                GenerateCSR();
                break;
        }
    }

    void Generator::Finalize(MPI_Comm comm) {
        switch (representation_) {
            case GraphRepresentation::EDGE_LIST:
                FinalizeEdgeList(comm);
                break;

            case GraphRepresentation::CSR:
                FinalizeCSR(comm);
                break;
        }
    }

    void Generator::FinalizeEdgeList(MPI_Comm comm) {}

    void Generator::FinalizeCSR(MPI_Comm) {}

    void CSROnlyGenerator::GenerateEdgeList() {
        GenerateCSR();
    }

    void CSROnlyGenerator::FinalizeEdgeList(MPI_Comm comm) {
        if (xadj_.empty()) {
            return;
        }

        // Otherwise, we have generated the graph in CSR representation, but
        // actually want edge list representation -> transform graph
        FinalizeCSR(comm);
        edges_ = BuildEdgeListFromCSR(vertex_range_, xadj_, adjncy_);
        {
            XadjArray tmp;
            std::swap(xadj_, tmp);
        }
        {
            AdjncyArray tmp;
            std::swap(adjncy_, tmp);
        }
    }

    void EdgeListOnlyGenerator::GenerateCSR() {
        GenerateEdgeList();
    }

    void EdgeListOnlyGenerator::FinalizeCSR(MPI_Comm comm) {
        if (!xadj_.empty()) {
            return;
        }

        // Otherwise, we have generated the graph in edge list representation, but
        // actually want CSR format --> transform graph
        FinalizeEdgeList(comm);
        std::tie(xadj_, adjncy_) = BuildCSRFromEdgeList(vertex_range_, edges_, edge_weights_);
        {
            Edgelist tmp;
            std::swap(edges_, tmp);
        }
    }

    std::unique_ptr<kagen::EdgeWeightGeneratorFactory> CreateEdgeWeightGeneratorFactory(const EdgeWeightType type) {
        switch (type) {
            case EdgeWeightType::HASHED:
                return std::make_unique<HashedEdgeWeightGeneratorFactory>();
            case EdgeWeightType::DISTANCE:
                return std::make_unique<DistanceEdgeWeightGeneratorFactory>();
            case EdgeWeightType::CONSTANT:
                return std::make_unique<ConstantEdgeWeightGeneratorFactory>();
            case EdgeWeightType::RANDOM:
                return std::make_unique<RandomEdgeWeightGeneratorFactory>();
        }

        throw std::runtime_error("invalid graph generator type");
    }

    void Generator::GenerateEdgeWeights(EdgeWeightType type, MPI_Comm comm) {
        std::unique_ptr<kagen::EdgeWeightGeneratorFactory> factory = CreateEdgeWeightGeneratorFactory(type);
        std::unique_ptr<kagen::EdgeWeightGenerator> generator = factory->Create();
        if (xadj_.empty()) {
            edge_weights_ = generator->GenerateEdgeWeights(edges_, coordinates_);
            std::cout << "Edgelist: edge weight size " << edge_weights_.size() << " edge weight first element "
                      << edge_weights_[0] << std::endl;
        } else {
            edge_weights_ = generator->GenerateEdgeWeights(xadj_, adjncy_, coordinates_);
            std::cout << "CSR: edge weight size " << edge_weights_.size() << " edge weight first element "
                      << edge_weights_[0] << std::endl;
        }

//        if(type == EdgeWeightType::RANDOM) {
//            ValidateEdgeWeight(comm);
//        }
    }

    void Generator::ValidateEdgeWeight(MPI_Comm comm) {
//        PEID size;
//        MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//        const int num_local_vertices = vertex_range_.second - vertex_range_.first;
//        std::vector<int> recvcounts(size);
//        MPI_Allgather(&num_local_vertices, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
//        std::vector<int> displs(size);
//        std::exclusive_scan(recvcounts.begin(), recvcounts.end(), displs.begin(), 0);
//
//        global_graph.vertex_weights.resize(recvcounts.back() + displs.back());
//        MPI_Allgatherv(
//                local_graph.vertex_weights.data(), num_local_vertices, KAGEN_MPI_SSINT,
//                global_graph.vertex_weights.data(),
//                recvcounts.data(), displs.data(), KAGEN_MPI_SSINT, MPI_COMM_WORLD);
//
//
//        int has_edge_weights = !local_graph.edge_weights.empty();
//        MPI_Allreduce(MPI_IN_PLACE, &has_edge_weights, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
//        if (has_edge_weights) {
//            const int num_local_edges = std::max<int>(local_graph.edges.size(), local_graph.adjncy.size());
//            std::vector<int> recvcounts(size);
//            MPI_Allgather(&num_local_edges, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
//            std::vector<int> displs(size);
//            std::exclusive_scan(recvcounts.begin(), recvcounts.end(), displs.begin(), 0);
//
//            global_graph.edge_weights.resize(recvcounts.back() + displs.back());
//            MPI_Allgatherv(
//                    local_graph.edge_weights.data(), num_local_edges, KAGEN_MPI_SSINT, global_graph.edge_weights.data(),
//                    recvcounts.data(), displs.data(), KAGEN_MPI_SSINT, MPI_COMM_WORLD);
//        }
    }

    SInt Generator::GetNumberOfEdges() const {
        return std::max(adjncy_.size(), edges_.size());
    }

    Graph Generator::Take() {
        return {
                vertex_range_,
                representation_,
                std::move(edges_),
                std::move(xadj_),
                std::move(adjncy_),
                std::move(vertex_weights_),
                std::move(edge_weights_),
                std::move(coordinates_)};
    }

    void Generator::SetVertexRange(const VertexRange vertex_range) {
        vertex_range_ = vertex_range;
    }

    void Generator::FilterDuplicateEdges() {
        std::sort(edges_.begin(), edges_.end());
        auto it = std::unique(edges_.begin(), edges_.end());
        edges_.erase(it, edges_.end());
    }

    void Generator::Reset() {
        edges_.clear();
        xadj_.clear();
        adjncy_.clear();
        vertex_weights_.clear();
        edge_weights_.clear();
        coordinates_.first.clear();
        coordinates_.second.clear();
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
            PGeneratorConfig &config, const PEID size, const bool output) const {
        if (config.k == 0) {
            if (IsSquare(size) && IsPowerOfTwo(size)) {
                config.k = static_cast<SInt>(size);
            } else {
                const SInt l = std::ceil(std::log2(size));
                config.k = 1 << l;
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
            PGeneratorConfig &config, const PEID size, const bool output) const {
        if (config.k == 0) {
            if (IsCubic(size) && IsPowerOfTwo(size)) {
                config.k = static_cast<SInt>(size);
            } else {
                const SInt l = std::ceil(std::log2(size));
                config.k = 1 << l;
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

    void GeneratorFactory::EnsureOneChunkPerPE(PGeneratorConfig &config, const PEID size) const {
        if (config.k != static_cast<SInt>(size)) {
            throw ConfigurationError("number of chunks must match the number of PEs");
        }
    }
}
// namespace kagen
