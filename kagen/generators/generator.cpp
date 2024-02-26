#include "kagen/generators/generator.h"

#include "kagen/context.h"
#include "kagen/kagen.h"
#include "kagen/tools/converter.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>

namespace kagen {
Generator::~Generator() = default;

void Generator::Generate(const GraphRepresentation representation) {
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
}

void Generator::Finalize(MPI_Comm comm) {
    switch (desired_representation_) {
        case GraphRepresentation::EDGE_LIST:
            FinalizeEdgeList(comm);
            break;

        case GraphRepresentation::CSR:
            FinalizeCSR(comm);
            break;
    }

    graph_.representation = desired_representation_;
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
    graph_.edges = BuildEdgeListFromCSR(graph_.vertex_range, graph_.xadj, graph_.adjncy);
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
    std::tie(graph_.xadj, graph_.adjncy) = BuildCSRFromEdgeList(graph_.vertex_range, graph_.edges, graph_.edge_weights);
    {
        Edgelist tmp;
        std::swap(graph_.edges, tmp);
    }
}

SInt Generator::GetNumberOfEdges() const {
    return std::max(graph_.adjncy.size(), graph_.edges.size());
}

Graph Generator::Take() {
    return std::move(graph_);
}

void Generator::SetVertexRange(const VertexRange vertex_range) {
    graph_.vertex_range = vertex_range;
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
