#include "kagen/context.h"
#include "kagen/generators/file/file_graph.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <mpi.h>

#include "io/parhip.h"
#include "tests/gather.h"
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <utility>

using namespace kagen;

using WeightRange   = std::pair<SInt, SInt>;
using GeneratorFunc = std::function<Graph(KaGen&, SInt, SInt)>;
MATCHER_P(EqualAdjacenyStructure, graph, "") {
    bool same_edge_list = graph.edges == arg.edges;
    bool same_xadj      = graph.xadj == arg.xadj;
    bool same_adjcny    = graph.adjncy == arg.adjncy;
    return same_edge_list && same_xadj && same_adjcny;
}

MATCHER_P(EqualWeights, graph, "") {
    bool same_edge_weights   = graph.edge_weights == arg.edge_weights;
    bool same_vertex_weights = graph.vertex_weights == arg.vertex_weights;
    return same_edge_weights && same_vertex_weights;
}

struct ParhipReadWriteTestFixture
    : public ::testing::TestWithParam<std::tuple<std::string, GeneratorFunc, GraphDistribution, GraphRepresentation>> {
};

INSTANTIATE_TEST_SUITE_P(
    GenericGeneratorTest, ParhipReadWriteTestFixture,
    ::testing::Values(
        std::make_tuple(
            "GNM", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateUndirectedGNM(n, m); }),
            GraphDistribution::BALANCE_VERTICES, GraphRepresentation::CSR),
        std::make_tuple(
            "RMAT", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRMAT(n, m, 0.56, 0.19, 0.19); }),
            GraphDistribution::BALANCE_VERTICES, GraphRepresentation::CSR),
        std::make_tuple(
            "GNM", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateUndirectedGNM(n, m); }),
            GraphDistribution::BALANCE_EDGES, GraphRepresentation::CSR),
        std::make_tuple(
            "RMAT", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRMAT(n, m, 0.56, 0.19, 0.19); }),
            GraphDistribution::BALANCE_EDGES, GraphRepresentation::CSR),
        std::make_tuple(
            "GNM", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateUndirectedGNM(n, m); }),
            GraphDistribution::ROOT, GraphRepresentation::CSR),
        std::make_tuple(
            "RMAT", GeneratorFunc([](KaGen& gen, SInt n, SInt m) { return gen.GenerateRMAT(n, m, 0.56, 0.19, 0.19); }),
            GraphDistribution::ROOT, GraphRepresentation::CSR)),
    [](const ::testing::TestParamInfo<ParhipReadWriteTestFixture::ParamType>& info) -> std::string {
        const std::string gen  = std::get<0>(info.param);
        const int         dist = static_cast<std::underlying_type_t<GraphDistribution>>(std::get<2>(info.param));
        const int         rep  = static_cast<std::underlying_type_t<GraphRepresentation>>(std::get<3>(info.param));
        return gen + "_" + std::to_string(dist) + "_" + std::to_string(rep);
    });

namespace {
template <typename Write>
void write_graph(Write&& write, int rank, int size) {
    bool continue_write = false;
    int  round          = 0;
    do {
        for (int i = 0; i < size; ++i) {
            if (i == rank) {
                continue_write = write(round);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        ++round;
    } while (continue_write);
}
std::filesystem::path get_temp_dir() {
    if (const char* runner_env = std::getenv("RUNNER_TEMP"); runner_env && *runner_env) {
        return std::filesystem::path(runner_env);
    }
    return std::filesystem::path(::testing::TempDir());
}

std::filesystem::path get_file_path(const std::string& instance_id) {
    return get_temp_dir() / (instance_id + ".parhip");
}

std::string get_instance_id(std::string const& test_name) {
    auto pos = test_name.rfind('/');
    return (pos == std::string::npos) ? test_name : test_name.substr(pos + 1);
}

SInt read_version(const std::string& filename) {
    std::ifstream in(filename, std::ios_base::binary);
    EXPECT_TRUE(in.good()) << "cannot open " << filename;
    parhip::ParhipID version = 0;
    in.read(reinterpret_cast<char*>(&version), sizeof(parhip::ParhipID));
    return static_cast<SInt>(version);
}
} // namespace
//

TEST_P(ParhipReadWriteTestFixture, default_write_read_in_parhip_format) {
    const ::testing::TestInfo* test_info   = ::testing::UnitTest::GetInstance()->current_test_info();
    const std::string          instance_id = get_instance_id(test_info->name());
    GeneratorFunc              generate    = std::get<1>(GetParam());
    const SInt                 n           = 1000;
    const SInt                 m           = 16 * n;
    const WeightRange          weight_range{1, 100};
    MPI_Comm                   comm = MPI_COMM_WORLD;
    int                        rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // setup
    kagen::KaGen generator(comm);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(
        kagen::EdgeWeightGeneratorType::UNIFORM_RANDOM, weight_range.first, weight_range.second);

    Graph generated_graph = generate(generator, n, m);

    // writer setup
    GraphInfo         info(generated_graph, comm);
    OutputGraphConfig config;
    config.filename = get_file_path(instance_id);
    kagen::ParhipWriter writer(config, generated_graph, info, rank, size);
    auto                write = [&](int round) {
        return writer.Write(round, config.filename);
    };
    write_graph(write, rank, size);

    // read graph
    const auto read_graph =
        generator.GenerateFromOptionString("type=file;filename=" + config.filename + ";edgeweights_generator=default");

    const auto total_generated_graph = kagen::testing::GatherGraph(generated_graph);
    const auto total_read_graph      = kagen::testing::GatherGraph(read_graph);
    EXPECT_THAT(total_read_graph, EqualAdjacenyStructure(total_generated_graph));
    EXPECT_THAT(total_read_graph, EqualWeights(total_generated_graph));
}

TEST_P(ParhipReadWriteTestFixture, write_from_csr_read_in_parhip_format) {
    const ::testing::TestInfo* test_info   = ::testing::UnitTest::GetInstance()->current_test_info();
    const std::string          instance_id = get_instance_id(test_info->name());
    GeneratorFunc              generate    = std::get<1>(GetParam());
    const SInt                 n           = 1000;
    const SInt                 m           = 16 * n;
    const WeightRange          weight_range{1, 100};
    MPI_Comm                   comm = MPI_COMM_WORLD;
    int                        rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // setup
    kagen::KaGen generator(comm);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(
        kagen::EdgeWeightGeneratorType::UNIFORM_RANDOM, weight_range.first, weight_range.second);

    Graph generated_graph = generate(generator, n, m);

    // writer setup
    GraphInfo         info(generated_graph, comm);
    OutputGraphConfig config;
    config.filename = get_file_path(instance_id);
    kagen::ParhipWriter writer(config, generated_graph, info, rank, size);
    auto                write = [&](int round) {
        return writer.WriteFromCSR<kagen::SInt, kagen::SSInt, kagen::SSInt>(
            round, config.filename, generated_graph.xadj, generated_graph.adjncy, nullptr,
            &generated_graph.edge_weights);
    };
    write_graph(write, rank, size);

    // read graph
    const auto read_graph =
        generator.GenerateFromOptionString("type=file;filename=" + config.filename + ";edgeweights_generator=default");

    const auto total_generated_graph = kagen::testing::GatherGraph(generated_graph);
    const auto total_read_graph      = kagen::testing::GatherGraph(read_graph);
    EXPECT_THAT(total_read_graph, EqualAdjacenyStructure(total_generated_graph));
    EXPECT_THAT(total_read_graph, EqualWeights(total_generated_graph));
}

// The same graph, once written with 32 bit and once with 64 bit vertex IDs (edge IDs are always 64 bit), must be read
// back identically -- independent of how the reader distributes the vertices among the PEs.
struct ParhipVertexIDWidthTestFixture : public ::testing::TestWithParam<std::string> {};

INSTANTIATE_TEST_SUITE_P(
    ParhipReaderTest, ParhipVertexIDWidthTestFixture,
    ::testing::Values(std::string("balance-vertices"), std::string("balance-edges")),
    [](const ::testing::TestParamInfo<ParhipVertexIDWidthTestFixture::ParamType>& info) -> std::string {
        std::string name = info.param;
        std::replace(name.begin(), name.end(), '-', '_');
        return name;
    });

TEST_P(ParhipVertexIDWidthTestFixture, read_graph_with_32bit_and_64bit_vertex_ids) {
    const ::testing::TestInfo* test_info    = ::testing::UnitTest::GetInstance()->current_test_info();
    const std::string          instance_id  = get_instance_id(test_info->name());
    const std::string          distribution = GetParam();
    const SInt                 n            = 1000;
    const SInt                 m            = 16 * n;
    const WeightRange          weight_range{1, 100};
    MPI_Comm                   comm = MPI_COMM_WORLD;
    int                        rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // setup
    kagen::KaGen generator(comm);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(
        kagen::EdgeWeightGeneratorType::UNIFORM_RANDOM, weight_range.first, weight_range.second);

    Graph generated_graph = generator.GenerateUndirectedGNM(n, m);

    // write the same graph twice, once with 32 bit and once with 64 bit vertex IDs
    const GraphInfo info(generated_graph, comm);

    auto write_with_vtx_width = [&](const int vtx_width) {
        OutputGraphConfig config;
        config.filename  = get_file_path(instance_id + "_" + std::to_string(vtx_width) + "bit");
        config.vtx_width = vtx_width;
        kagen::ParhipWriter writer(config, generated_graph, info, rank, size);
        write_graph([&](int round) { return writer.Write(round, config.filename); }, rank, size);
        return config.filename;
    };

    const std::string filename_32bit = write_with_vtx_width(32);
    const std::string filename_64bit = write_with_vtx_width(64);

    // the two files must only differ in their vertex ID width
    const SInt version_32bit = read_version(filename_32bit);
    const SInt version_64bit = read_version(filename_64bit);
    EXPECT_TRUE(parhip::Has32BitVertexIDs(version_32bit));
    EXPECT_FALSE(parhip::Has32BitVertexIDs(version_64bit));
    EXPECT_FALSE(parhip::Has32BitEdgeIDs(version_32bit));
    EXPECT_FALSE(parhip::Has32BitEdgeIDs(version_64bit));
    EXPECT_EQ(std::filesystem::file_size(filename_64bit) - std::filesystem::file_size(filename_32bit), 4 * 2 * m);

    // read both graphs back with the same distribution
    auto read = [&](const std::string& filename) {
        return generator.GenerateFromOptionString(
            "type=file;filename=" + filename + ";distribution=" + distribution + ";edgeweights_generator=default");
    };
    const auto read_graph_32bit = read(filename_32bit);
    const auto read_graph_64bit = read(filename_64bit);

    // both files must be distributed in the same way ...
    EXPECT_EQ(read_graph_32bit.vertex_range, read_graph_64bit.vertex_range);

    // ... and yield the graph that was written
    const auto total_generated_graph  = kagen::testing::GatherGraph(generated_graph);
    const auto total_read_graph_32bit = kagen::testing::GatherGraph(read_graph_32bit);
    const auto total_read_graph_64bit = kagen::testing::GatherGraph(read_graph_64bit);

    EXPECT_THAT(total_read_graph_64bit, EqualAdjacenyStructure(total_generated_graph));
    EXPECT_THAT(total_read_graph_64bit, EqualWeights(total_generated_graph));
    EXPECT_THAT(total_read_graph_32bit, EqualAdjacenyStructure(total_generated_graph));
    EXPECT_THAT(total_read_graph_32bit, EqualWeights(total_generated_graph));
}

TEST_P(ParhipReadWriteTestFixture, write_from_csr_read_in_parhip_format_32bit_edges) {
    const ::testing::TestInfo* test_info   = ::testing::UnitTest::GetInstance()->current_test_info();
    const std::string          instance_id = get_instance_id(test_info->name());
    GeneratorFunc              generate    = std::get<1>(GetParam());
    const SInt                 n           = 1000;
    const SInt                 m           = 16 * n;
    const WeightRange          weight_range{1, 100};
    MPI_Comm                   comm = MPI_COMM_WORLD;
    int                        rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // setup
    kagen::KaGen generator(comm);
    generator.UseCSRRepresentation();
    generator.ConfigureEdgeWeightGeneration(
        kagen::EdgeWeightGeneratorType::UNIFORM_RANDOM, weight_range.first, weight_range.second);
    Graph generated_graph = generate(generator, n, m);

    // writer setup
    GraphInfo         info(generated_graph, comm);
    OutputGraphConfig config;
    config.filename     = get_file_path(instance_id);
    config.adjwgt_width = 32;
    kagen::ParhipWriter       writer(config, generated_graph, info, rank, size);
    std::vector<std::int32_t> edge_weight_32_bit(generated_graph.edge_weights.size());
    std::fill(edge_weight_32_bit.begin(), edge_weight_32_bit.end(), static_cast<std::int32_t>(size));

    auto write = [&](int round) {
        return writer.WriteFromCSR<kagen::SInt, kagen::SSInt, std::int32_t>(
            round, config.filename, generated_graph.xadj, generated_graph.adjncy, nullptr, &edge_weight_32_bit);
    };
    write_graph(write, rank, size);

    // read graph
    generator.UseCSRRepresentation();
    const auto read_graph =
        generator.GenerateFromOptionString("type=file;filename=" + config.filename + ";edgeweights_generator=default");
    const auto total_generated_graph = kagen::testing::GatherGraph(generated_graph);
    const auto total_read_graph      = kagen::testing::GatherGraph(read_graph);
    EXPECT_THAT(total_read_graph, EqualAdjacenyStructure(total_generated_graph));
    EXPECT_THAT(total_read_graph.edge_weights, ::testing::Each(size));
    EXPECT_TRUE(total_read_graph.vertex_weights.empty());
}
