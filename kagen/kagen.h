#pragma once

#ifdef __cplusplus
    #include <algorithm>
    #include <cstdint>
    #include <memory>
    #include <string>
    #include <tuple>
    #include <type_traits>
    #include <unordered_map>
    #include <utility>
    #include <vector>
#endif

#include <mpi.h>
#include <stdbool.h>
#include <stddef.h>

#define KAGEN_VERSION_MAJOR 1
#define KAGEN_VERSION_MINOR 0
#define KAGEN_VERSION_PATCH 2

#ifdef __cplusplus
namespace kagen {
std::string BuildDescription();
}
#endif

//
// KaGen data types
//

#define KAGEN_MPI_UINT    MPI_UNSIGNED
#define KAGEN_MPI_PEID    MPI_INT
#define KAGEN_MPI_HPFLOAT MPI_LONG_DOUBLE
#define KAGEN_MPI_LPFLOAT MPI_DOUBLE
#define KAGEN_MPI_SINT    MPI_UNSIGNED_LONG_LONG
#define KAGEN_MPI_SSINT   MPI_LONG_LONG

// C++ interface
#ifdef __cplusplus
namespace kagen {
using SInt          = unsigned long long;
using SSInt         = long long;
using Edgelist      = std::vector<std::pair<SInt, SInt>>;
using Edgelist32    = std::vector<std::pair<int, int>>;
using VertexRange   = std::pair<SInt, SInt>;
using PEID          = int;
using HPFloat       = long double;
using LPFloat       = double;
using Coordinates2D = std::vector<std::tuple<HPFloat, HPFloat>>;
using Coordinates3D = std::vector<std::tuple<HPFloat, HPFloat, HPFloat>>;
using Coordinates   = std::pair<Coordinates2D, Coordinates3D>;
using VertexWeights = std::vector<SSInt>;
using EdgeWeights   = std::vector<SSInt>;
using XadjArray     = std::vector<SInt>;
using AdjncyArray   = std::vector<SInt>;
} // namespace kagen
#endif

// C interface
typedef unsigned long long kagen_index;
typedef long long          kagen_weight;

//
// Config enums (C++ interface only)
//

#ifdef __cplusplus
namespace kagen {
enum class FileFormat {
    NOOP,
    EXTENSION,
    EDGE_LIST,
    EDGE_LIST_UNDIRECTED,
    BINARY_EDGE_LIST,
    BINARY_EDGE_LIST_UNDIRECTED,
    PLAIN_EDGE_LIST,
    METIS,
    HMETIS,
    HMETIS_DIRECTED,
    HMETIS_EP,
    DOT,
    DOT_DIRECTED,
    COORDINATES,
    PARHIP,
    XTRAPULP,
    FREIGHT_NETL,
    FREIGHT_NETL_EP,
};

std::unordered_map<std::string, FileFormat> GetOutputFormatMap();

std::unordered_map<std::string, FileFormat> GetInputFormatMap();

std::ostream& operator<<(std::ostream& out, FileFormat format);

enum class GeneratorType {
    GNM_DIRECTED,
    GNM_UNDIRECTED,
    GNP_DIRECTED,
    GNP_UNDIRECTED,
    RGG_2D,
    RGG_3D,
    RDG_2D,
    RDG_3D,
    GRID_2D,
    GRID_3D,
    PATH_DIRECTED,
    BA,
    KRONECKER,
    RHG,
    RMAT,
    IMAGE_MESH,
    FILE,
};

std::unordered_map<std::string, GeneratorType> GetGeneratorTypeMap();

std::ostream& operator<<(std::ostream& out, GeneratorType generator_type);

enum class StatisticsLevel : std::uint8_t {
    NONE     = 0,
    BASIC    = 1,
    ADVANCED = 2,
};

bool operator<=(StatisticsLevel a, StatisticsLevel b);

std::unordered_map<std::string, StatisticsLevel> GetStatisticsLevelMap();

std::ostream& operator<<(std::ostream& out, StatisticsLevel statistics_level);

enum class ImageMeshWeightModel : std::uint8_t {
    L2         = 0,
    INV_L2     = 1,
    RATIO      = 2,
    INV_RATIO  = 3,
    SIMILARITY = 4,
};

std::unordered_map<std::string, ImageMeshWeightModel> GetImageMeshWeightModelMap();

std::ostream& operator<<(std::ostream& out, ImageMeshWeightModel weight_model);

enum class GraphRepresentation {
    EDGE_LIST,
    CSR,
};

std::unordered_map<std::string, GraphRepresentation> GetGraphRepresentationMap();

std::ostream& operator<<(std::ostream& out, GraphRepresentation representation);

enum class GraphDistribution {
    ROOT,
    BALANCE_VERTICES,
    BALANCE_EDGES,
};

std::unordered_map<std::string, GraphDistribution> GetGraphDistributionMap();

std::ostream& operator<<(std::ostream& out, GraphDistribution distribution);
} // namespace kagen
#endif

//
// C++ interface
//

#ifdef __cplusplus
namespace kagen {
struct Graph {
    VertexRange         vertex_range;
    GraphRepresentation representation;

    // Edge list representation only
    Edgelist edges;

    // CSR representation only
    XadjArray   xadj;
    AdjncyArray adjncy;

    VertexWeights vertex_weights;
    EdgeWeights   edge_weights;
    Coordinates   coordinates;

    SInt NumberOfLocalVertices() const;

    SInt NumberOfLocalEdges() const;

    void SortEdgelist();

    void FreeEdgelist();
    void FreeCSR();

    void Clear();

    template <typename T = SInt>
    std::vector<std::pair<T, T>> TakeEdges() {
        return TakeVector<std::pair<T, T>>(edges);
    }

    template <typename T = SInt>
    std::vector<T> TakeXadj() {
        return TakeVector<T>(xadj);
    }

    template <typename T = SInt>
    std::vector<T> TakeAdjncy() {
        return TakeVector<T>(adjncy);
    }

    template <typename T = SSInt>
    std::vector<T> TakeVertexWeights() {
        return TakeVector<T>(vertex_weights);
    }

    template <typename T = SSInt>
    std::vector<T> TakeEdgeWeights() {
        return TakeVector<T>(edge_weights);
    }

    std::tuple<VertexRange, Edgelist, XadjArray, AdjncyArray, VertexWeights, EdgeWeights, Coordinates> tuple() && {
        return std::make_tuple(
            vertex_range, std::move(edges), std::move(xadj), std::move(adjncy), std::move(vertex_weights),
            std::move(edge_weights), std::move(coordinates));
    }

private:
    template <typename To, typename From, std::enable_if_t<std::is_same_v<typename From::value_type, To>, bool> = true>
    std::vector<To> TakeVector(From& from) {
        return std::move(from);
    }

    template <typename To, typename From, std::enable_if_t<!std::is_same_v<typename From::value_type, To>, bool> = true>
    std::vector<To> TakeVector(From& from) {
        std::vector<To> copy(from.size());
        std::copy(from.begin(), from.end(), copy.begin());
        std::vector<typename From::value_type> free;
        std::swap(from, free);
        return copy;
    }
};

class KaGen {
public:
    KaGen(MPI_Comm comm);

    KaGen(const KaGen&) = delete;
    KaGen(KaGen&&) noexcept;

    KaGen& operator=(const KaGen&) = delete;
    KaGen& operator=(KaGen&&) noexcept;

    ~KaGen();

    /*!
     * Sets the seed for the random number generator (Default: 0).
     * Must be the same on all PEs.
     * @param seed Seed for random number generation.
     */
    void SetSeed(int seed);

    /*!
     * If enabled, KaGen will verify that the generated graph is simple and undirected.
     * This feature is only for debugging; unless explicitly configured otherwise, KaGen
     * should *never* fail this verification.
     */
    void EnableUndirectedGraphVerification();

    void EnableBasicStatistics();

    void EnableAdvancedStatistics();

    /*!
     * If enabled, KaGen will print information to stdout and stderr (but only on rank 0).
     *
     * @param header If set to true, KaGen will also print a banner and a summary of the
     * configuration parameters.
     */
    void EnableOutput(bool header);

    /*!
     * If set to true, KaGen will use higher-precision floating point numbers (80 bit instead of 64 bit on x86 systems)
     * to generate graphs. Currently, this only affects the random hyperbolic graph generator. Per default, KaGen will
     * decide automatically which precision to use.
     *
     * Note that "higher-precision" refers to "long double". On systems where "long double" has the same precision as
     * "double", this option does nothing.
     *
     * @param state If true, always use higher-precision floating point numbers; if false, never use them.
     */
    void UseHPFloats(bool state);

    /*!
     * Controls the number of chunks KaGen uses for graph generation. If not set explicitly, KaGen will choose the
     * number of chunks automatically. Usually, you do not have to use this option.
     *
     * @param k Number of chunks to be used.
     */
    void SetNumberOfChunks(SInt k);

    /*!
     * Represents the generated graph as a list of edges (from, to).
     * This representation requires 2 * |E| memory and is the default representation.
     */
    void UseEdgeListRepresentation();

    /*!
     * Represents the generated graph in compressed sparse row format.
     * This representation requires |V| + |E| memory. However, not all generators support this representation;
     * for generators that do not support it, KaGen first generates the graph as a list of edges, sorts it and
     * then builds the CSR data structures. This requires |V|+3|E| memory and |E|*log(|E|) time.
     */
    void UseCSRRepresentation();

    /*!
     * Generates a graph with options given by a string of options in `key=value` or `flag` format:
     * `key1=value1;flag1;...`.
     *
     * Use one of the following flags to select a graph model:
     * - gnm_undirected
     * - gnm_directed
     * - gnp_undirected
     * - gnp_directed
     * - rgg2d
     * - rgg3d
     * - grid2d
     * - grid3d
     * - rdg2d
     * - rdg3d
     * - rhg
     * - ba
     * - kronecker
     * - rmat
     *
     * Use the following keys to specify generator properties:
     * - n=<SInt>             -- number of nodes
     * - N=<SInt>             -- number of nodes as a power of 2
     * - m=<SInt>             -- number of edges
     * - M=<SInt>             -- number of edges as a power of 2
     * - k=<SInt>             -- number of chunks
     * - prob=<HPFloat>       -- edge probability (varius generators)
     * - radius=<HPFloat>     -- edge radius (RGG2D/3D)
     * - gamma=<HPFloat>      -- power law exponent (RHG)
     * - avg_degree=<HPFloat> -- average degree (RHG)
     * - min_degree=<SInt>    -- minimum degree (BA)
     * - grid_x=<SInt>        -- grid width (GRID2D/3D)
     * - grid_y=<SInt>        -- grid height (GRID2D/3D)
     * - grid_z=<SInt>        -- grid depth (GRID3D)
     * - rmat_a=<HPFloat>     -- RMat probability for block A (RMAT)
     * - rmat_b=<HPFloat>     -- RMat probability for block B (RMAT)
     * - rmat_c=<HPFloat>     -- RMat probability for block C (RMAT)
     * - periodic[=0|1]       -- periodic boundary condition (various generators)
     * - coordinates          -- generate coordinates (only for geometric graph generators)
     *
     * Depending on the selected generator type, some options are mandatory, some are optional and some are ignored.
     * The following example generates a RGG2D graph with 100 nodes and 200 edges: `rgg2d;n=100;m=200`.
     *
     * @type options Options string with key=value pairs.
     * @return The generated graph.
     */
    Graph GenerateFromOptionString(const std::string& options);

    Graph GenerateDirectedGNM(SInt n, SInt m, bool self_loops = false);

    Graph GenerateUndirectedGNM(SInt n, SInt m, bool self_loops = false);

    Graph GenerateDirectedGNP(SInt n, LPFloat p, bool self_loops = false);

    Graph GenerateUndirectedGNP(SInt n, LPFloat p, bool self_loops = false);

    Graph GenerateRGG2D(SInt n, LPFloat r, bool coordinates = false);

    Graph GenerateRGG2D_NM(SInt n, SInt m, bool coordinates = false);

    Graph GenerateRGG2D_MR(SInt m, LPFloat r, bool coordinates = false);

    Graph GenerateRGG3D(SInt n, LPFloat r, bool coordinates = false);

    Graph GenerateRGG3D_NM(SInt n, SInt m, bool coordinates = false);

    Graph GenerateRGG3D_MR(SInt m, LPFloat r, bool coordinates = false);

    Graph GenerateRDG2D(SInt n, bool periodic, bool coordinates = false);

    Graph GenerateRDG2D_M(SInt m, bool periodic, bool coordinates = false);

    Graph GenerateRDG3D(SInt n, bool coordinates = false);

    Graph GenerateRDG3D_M(SInt m, bool coordinates = false);

    Graph GenerateBA(SInt n, SInt d, bool directed = false, bool self_loops = false);

    Graph GenerateBA_NM(SInt n, SInt m, bool directed = false, bool self_loops = false);

    Graph GenerateBA_MD(SInt m, SInt d, bool directed = false, bool self_loops = false);

    Graph GenerateRHG(LPFloat gamma, SInt n, LPFloat d, bool coordinates = false);

    Graph GenerateRHG_NM(LPFloat gamma, SInt n, SInt m, bool coordinates = false);

    Graph GenerateRHG_MD(LPFloat gamma, SInt m, LPFloat d, bool coordinates = false);

    Graph GenerateGrid2D(SInt grid_x, SInt grid_y, LPFloat p, bool periodic = false, bool coordinates = false);

    Graph GenerateGrid2D_N(SInt n, LPFloat p, bool periodic = false, bool coordinates = false);

    Graph GenerateGrid2D_NM(SInt n, SInt m, bool periodic = false, bool coordinates = false);

    Graph
    GenerateGrid3D(SInt grid_x, SInt grid_y, SInt grid_z, LPFloat p, bool periodic = false, bool coordinates = false);

    Graph GenerateGrid3D_N(SInt n, LPFloat p, bool periodic = false, bool coordinates = false);

    Graph GenerateGrid3D_NM(SInt n, SInt m, bool periodic = false, bool coordinates = false);

    Graph GenerateDirectedPath(unsigned long long n, bool permute = false, bool periodic = false);

    Graph GenerateKronecker(SInt n, SInt m, bool directed = false, bool self_loops = false);

    Graph GenerateRMAT(SInt n, SInt m, LPFloat a, LPFloat b, LPFloat c, bool directed = false, bool self_loops = false);

    Graph ReadFromFile(std::string const& filename, const FileFormat format, const GraphDistribution distribution);

private:
    void SetDefaults();

    MPI_Comm                                 comm_;
    std::unique_ptr<struct PGeneratorConfig> config_;
    GraphRepresentation                      representation_;
};

/*!
 * Returns an array A of size |comm| + 1 s.t. PE i contains nodes in the range [A[i], A[i + 1]).
 *
 * @param graph Graph generated by any of the Generate* functions.
 * @param idx_mpi_type MPI type corresponding to the template parameter IDX.
 * @param comm MPI communicator that was used to generate the graph.
 *
 * @tparam IDX Data type to be used for the entries in A. Must be large enough to represent the global number of nodes
 * in the graph.
 *
 * @return Vertex distribution as described above.
 */
template <typename IDX>
std::vector<IDX> BuildVertexDistribution(const Graph& graph, MPI_Datatype idx_mpi_type, MPI_Comm comm) {
    PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::vector<IDX> distribution(size + 1);
    distribution[rank + 1] = graph.vertex_range.second;
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, distribution.data() + 1, 1, idx_mpi_type, comm);

    return distribution;
}
} // namespace kagen
#endif

//
// C interface
//

typedef struct kagen_obj   kagen_obj;
typedef struct kagen_graph kagen_graph;

typedef struct {
    kagen_index source;
    kagen_index target;
} kagen_edge;

#ifdef __cplusplus
extern "C" {
#endif

kagen_obj* kagen_create(MPI_Comm comm);
void       kagen_free(kagen_obj* gen);

void          kagen_graph_vertex_range(kagen_graph* result, kagen_index* begin, kagen_index* end);
kagen_edge*   kagen_graph_edge_list(kagen_graph* result, size_t* nedges);
kagen_index*  kagen_graph_csr_xadj(kagen_graph* result, size_t* nvertices);
kagen_index*  kagen_graph_csr_adjncy(kagen_graph* result, size_t* nedges);
kagen_weight* kagen_graph_vertex_weights(kagen_graph* result, size_t* size);
kagen_weight* kagen_graph_edge_weights(kagen_graph* result, size_t* size);
void          kagen_graph_free(kagen_graph* result);

void kagen_set_seed(kagen_obj* gen, int seed);
void kagen_enable_undirected_graph_verification(kagen_obj* gen);
void kagen_enable_basic_statistics(kagen_obj* gen);
void kagen_enable_advanced_statistics(kagen_obj* gen);
void kagen_enable_output(kagen_obj* gen, bool header);
void kagen_use_hp_floats(kagen_obj* gen, bool state);
void kagen_set_numer_of_chunks(kagen_obj* gen, unsigned long long k);
void kagen_use_edge_list_representation(kagen_obj* gen);
void kagen_use_csr_representation(kagen_obj* gen);

kagen_graph* kagen_generate_from_option_string(kagen_obj* gen, const char* options);

kagen_graph* kagen_generate_directed_gnm(kagen_obj* gen, unsigned long long n, unsigned long long m, bool self_loops);
kagen_graph* kagen_generate_undirected_gnm(kagen_obj* gen, unsigned long long n, unsigned long long m, bool self_loops);
kagen_graph* kagen_generate_directed_gnp(kagen_obj* gen, unsigned long long n, double p, bool self_loops);
kagen_graph* kagen_generate_undirected_gnp(kagen_obj* gen, unsigned long long n, double p, bool self_loops);

kagen_graph* kagen_generate_rgg2d(kagen_obj* gen, unsigned long long n, double r);
kagen_graph* kagen_generate_rgg2d_nm(kagen_obj* gen, unsigned long long n, unsigned long long m);
kagen_graph* kagen_generate_rgg2d_mr(kagen_obj* gen, unsigned long long m, double r);

kagen_graph* kagen_generate_rgg3d(kagen_obj* gen, unsigned long long n, double r);
kagen_graph* kagen_generate_rgg3d_nm(kagen_obj* gen, unsigned long long n, unsigned long long m);
kagen_graph* kagen_generate_rgg3d_mr(kagen_obj* gen, unsigned long long m, double r);

kagen_graph* kagen_generate_rdg2d(kagen_obj* gen, unsigned long long n, bool periodic);
kagen_graph* kagen_generate_rdg2d_m(kagen_obj* gen, unsigned long long m, bool periodic);
kagen_graph* kagen_generate_rdg3d(kagen_obj* gen, unsigned long long n);
kagen_graph* kagen_generate_rdg3d_m(kagen_obj* gen, unsigned long long m);

kagen_graph*
kagen_generate_ba(kagen_obj* gen, unsigned long long n, unsigned long long d, bool directed, bool self_loops);
kagen_graph*
kagen_generate_ba_nm(kagen_obj* gen, unsigned long long n, unsigned long long m, bool directed, bool self_loops);
kagen_graph*
kagen_generate_ba_md(kagen_obj* gen, unsigned long long m, unsigned long long d, bool directed, bool self_loops);

kagen_graph* kagen_generate_rhg(kagen_obj* gen, double gamma, unsigned long long n, double d);
kagen_graph* kagen_generate_rhg_nm(kagen_obj* gen, double gamma, unsigned long long n, unsigned long long m);
kagen_graph* kagen_generate_rhg_md(kagen_obj* gen, double gamma, unsigned long long m, double d);

kagen_graph*
kagen_generate_grid2d(kagen_obj* gen, unsigned long long grid_x, unsigned long long grid_y, double p, bool periodic);
kagen_graph* kagen_generate_grid2d_n(kagen_obj* gen, unsigned long long n, double p, bool periodic);
kagen_graph* kagen_generate_grid3d(
    kagen_obj* gen, unsigned long long grid_x, unsigned long long grid_y, unsigned long long grid_z, double p,
    bool periodic);
kagen_graph* kagen_generate_grid3d_n(kagen_obj* gen, unsigned long long n, double p, bool periodic);

kagen_graph* kagen_generate_directed_path(kagen_obj* gen, unsigned long long n, bool permuted, bool periodic);

kagen_graph*
kagen_generate_kronecker(kagen_obj* gen, unsigned long long n, unsigned long long m, bool directed, bool self_loops);

kagen_graph* kagen_generate_rmat(
    kagen_obj* gen, unsigned long long n, unsigned long long m, double a, double b, double c, bool directed,
    bool self_loops);

void kagen_build_vertex_distribution(kagen_graph* result, kagen_index* dist, MPI_Comm comm);

#ifdef __cplusplus
}
#endif
