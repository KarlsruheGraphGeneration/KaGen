#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <kagen/kagen.h>

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;

// ============================================================================
// PyGraph: wraps kagen::Graph, exposes data as numpy arrays
// ============================================================================

struct PyGraph {
    kagen::Graph graph;

    explicit PyGraph(kagen::Graph&& g) : graph(std::move(g)) {}

    py::tuple vertex_range() const {
        return py::make_tuple(
            static_cast<uint64_t>(graph.vertex_range.first),
            static_cast<uint64_t>(graph.vertex_range.second));
    }

    uint64_t num_vertices() const {
        return static_cast<uint64_t>(graph.vertex_range.second - graph.vertex_range.first);
    }

    uint64_t num_edges() const {
        if (graph.representation == kagen::GraphRepresentation::EDGE_LIST) {
            return static_cast<uint64_t>(graph.edges.size());
        } else {
            return static_cast<uint64_t>(graph.adjncy.size());
        }
    }

    // Edge list as (N, 2) uint64 array
    py::array_t<uint64_t> edges() const {
        const auto n = graph.edges.size();
        if (n == 0) {
            std::vector<py::ssize_t> shape = {0, 2};
            return py::array_t<uint64_t>(shape);
        }
        // pair<SInt,SInt> is layout-compatible with uint64_t[2]
        std::vector<py::ssize_t> shape = {static_cast<py::ssize_t>(n), 2};
        auto arr = py::array_t<uint64_t>(shape);
        auto buf = arr.mutable_unchecked<2>();
        for (size_t i = 0; i < n; ++i) {
            buf(i, 0) = static_cast<uint64_t>(graph.edges[i].first);
            buf(i, 1) = static_cast<uint64_t>(graph.edges[i].second);
        }
        return arr;
    }

    // CSR xadj array
    py::array_t<uint64_t> xadj() const {
        const auto n = graph.xadj.size();
        auto arr = py::array_t<uint64_t>(static_cast<py::ssize_t>(n));
        auto buf = arr.mutable_unchecked<1>();
        for (size_t i = 0; i < n; ++i) {
            buf(i) = static_cast<uint64_t>(graph.xadj[i]);
        }
        return arr;
    }

    // CSR adjncy array
    py::array_t<uint64_t> adjncy() const {
        const auto n = graph.adjncy.size();
        auto arr = py::array_t<uint64_t>(static_cast<py::ssize_t>(n));
        auto buf = arr.mutable_unchecked<1>();
        for (size_t i = 0; i < n; ++i) {
            buf(i) = static_cast<uint64_t>(graph.adjncy[i]);
        }
        return arr;
    }

    py::array_t<int64_t> vertex_weights() const {
        const auto n = graph.vertex_weights.size();
        auto arr = py::array_t<int64_t>(static_cast<py::ssize_t>(n));
        auto buf = arr.mutable_unchecked<1>();
        for (size_t i = 0; i < n; ++i) {
            buf(i) = static_cast<int64_t>(graph.vertex_weights[i]);
        }
        return arr;
    }

    py::array_t<int64_t> edge_weights() const {
        const auto n = graph.edge_weights.size();
        auto arr = py::array_t<int64_t>(static_cast<py::ssize_t>(n));
        auto buf = arr.mutable_unchecked<1>();
        for (size_t i = 0; i < n; ++i) {
            buf(i) = static_cast<int64_t>(graph.edge_weights[i]);
        }
        return arr;
    }

    // 2D coordinates as (N, 2) float64 array (downcast from long double)
    py::array_t<double> coordinates_2d() const {
        const auto& coords = graph.coordinates.first;
        const auto n = coords.size();
        if (n == 0) {
            std::vector<py::ssize_t> shape = {0, 2};
            return py::array_t<double>(shape);
        }
        std::vector<py::ssize_t> shape2d = {static_cast<py::ssize_t>(n), 2};
        auto arr = py::array_t<double>(shape2d);
        auto buf = arr.mutable_unchecked<2>();
        for (size_t i = 0; i < n; ++i) {
            buf(i, 0) = static_cast<double>(std::get<0>(coords[i]));
            buf(i, 1) = static_cast<double>(std::get<1>(coords[i]));
        }
        return arr;
    }

    // 3D coordinates as (N, 3) float64 array (downcast from long double)
    py::array_t<double> coordinates_3d() const {
        const auto& coords = graph.coordinates.second;
        const auto n = coords.size();
        if (n == 0) {
            std::vector<py::ssize_t> shape = {0, 3};
            return py::array_t<double>(shape);
        }
        std::vector<py::ssize_t> shape3d = {static_cast<py::ssize_t>(n), 3};
        auto arr = py::array_t<double>(shape3d);
        auto buf = arr.mutable_unchecked<2>();
        for (size_t i = 0; i < n; ++i) {
            buf(i, 0) = static_cast<double>(std::get<0>(coords[i]));
            buf(i, 1) = static_cast<double>(std::get<1>(coords[i]));
            buf(i, 2) = static_cast<double>(std::get<2>(coords[i]));
        }
        return arr;
    }

    void sort_edgelist() {
        graph.SortEdgelist();
    }
};

// ============================================================================
// Helper: parse edge/vertex weight generator type from string
// ============================================================================

static kagen::EdgeWeightGeneratorType parse_edge_weight_type(const std::string& s) {
    auto map = kagen::GetEdgeWeightGeneratorTypeMap();
    auto it = map.find(s);
    if (it == map.end()) {
        throw py::value_error("Unknown edge weight generator type: '" + s + "'");
    }
    return it->second;
}

static kagen::VertexWeightGeneratorType parse_vertex_weight_type(const std::string& s) {
    auto map = kagen::GetVertexWeightGeneratorTypeMap();
    auto it = map.find(s);
    if (it == map.end()) {
        throw py::value_error("Unknown vertex weight generator type: '" + s + "'");
    }
    return it->second;
}

// ============================================================================
// PyKaGen: wraps kagen::KaGen
// ============================================================================

struct PyKaGen {
    kagen::KaGen gen;

    PyKaGen() : gen(MPI_COMM_WORLD) {}

    void set_seed(int seed) { gen.SetSeed(seed); }
    void enable_undirected_graph_verification() { gen.EnableUndirectedGraphVerification(); }
    void enable_basic_statistics() { gen.EnableBasicStatistics(); }
    void enable_advanced_statistics() { gen.EnableAdvancedStatistics(); }
    void enable_vertex_permutation() { gen.EnableVertexPermutation(); }

    void configure_edge_weight_generation(const std::string& generator, uint64_t begin, uint64_t end) {
        gen.ConfigureEdgeWeightGeneration(parse_edge_weight_type(generator), begin, end);
    }

    void configure_vertex_weight_generation(const std::string& generator, uint64_t begin, uint64_t end) {
        gen.ConfigureVertexWeightGeneration(parse_vertex_weight_type(generator), begin, end);
    }

    void enable_output(bool header) { gen.EnableOutput(header); }
    void use_hp_floats(bool state) { gen.UseHPFloats(state); }
    void set_num_threads(int threads) { gen.SetNumberOfThreads(threads); }
    void set_number_of_chunks(uint64_t k) { gen.SetNumberOfChunks(k); }
    void use_edge_list_representation() { gen.UseEdgeListRepresentation(); }
    void use_csr_representation() { gen.UseCSRRepresentation(); }

    // ---- Generators ----

    PyGraph generate_from_option_string(const std::string& options) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateFromOptionString(options));
    }

    PyGraph generate_directed_gnm(uint64_t n, uint64_t m, bool self_loops) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateDirectedGNM(n, m, self_loops));
    }

    PyGraph generate_undirected_gnm(uint64_t n, uint64_t m, bool self_loops) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateUndirectedGNM(n, m, self_loops));
    }

    PyGraph generate_directed_gnp(uint64_t n, double p, bool self_loops) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateDirectedGNP(n, p, self_loops));
    }

    PyGraph generate_undirected_gnp(uint64_t n, double p, bool self_loops) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateUndirectedGNP(n, p, self_loops));
    }

    PyGraph generate_rgg2d(uint64_t n, double r, bool coordinates) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateRGG2D(n, r, coordinates));
    }

    PyGraph generate_rgg2d_nm(uint64_t n, uint64_t m, bool coordinates) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateRGG2D_NM(n, m, coordinates));
    }

    PyGraph generate_rgg2d_mr(uint64_t m, double r, bool coordinates) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateRGG2D_MR(m, r, coordinates));
    }

    PyGraph generate_rgg3d(uint64_t n, double r, bool coordinates) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateRGG3D(n, r, coordinates));
    }

    PyGraph generate_rgg3d_nm(uint64_t n, uint64_t m, bool coordinates) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateRGG3D_NM(n, m, coordinates));
    }

    PyGraph generate_rgg3d_mr(uint64_t m, double r, bool coordinates) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateRGG3D_MR(m, r, coordinates));
    }

    PyGraph generate_rdg2d(uint64_t n, bool periodic, bool coordinates) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateRDG2D(n, periodic, coordinates));
    }

    PyGraph generate_rdg2d_m(uint64_t m, bool periodic, bool coordinates) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateRDG2D_M(m, periodic, coordinates));
    }

    PyGraph generate_rdg3d(uint64_t n, bool coordinates) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateRDG3D(n, coordinates));
    }

    PyGraph generate_rdg3d_m(uint64_t m, bool coordinates) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateRDG3D_M(m, coordinates));
    }

    PyGraph generate_ba(uint64_t n, uint64_t d, bool directed, bool self_loops) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateBA(n, d, directed, self_loops));
    }

    PyGraph generate_ba_nm(uint64_t n, uint64_t m, bool directed, bool self_loops) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateBA_NM(n, m, directed, self_loops));
    }

    PyGraph generate_ba_md(uint64_t m, uint64_t d, bool directed, bool self_loops) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateBA_MD(m, d, directed, self_loops));
    }

    PyGraph generate_rhg(double gamma, uint64_t n, double d, bool coordinates) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateRHG(gamma, n, d, coordinates));
    }

    PyGraph generate_rhg_nm(double gamma, uint64_t n, uint64_t m, bool coordinates) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateRHG_NM(gamma, n, m, coordinates));
    }

    PyGraph generate_rhg_md(double gamma, uint64_t m, double d, bool coordinates) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateRHG_MD(gamma, m, d, coordinates));
    }

    PyGraph generate_grid2d(uint64_t x, uint64_t y, double p, bool periodic, bool coordinates) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateGrid2D(x, y, p, periodic, coordinates));
    }

    PyGraph generate_grid2d_n(uint64_t n, double p, bool periodic, bool coordinates) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateGrid2D_N(n, p, periodic, coordinates));
    }

    PyGraph generate_grid2d_nm(uint64_t n, uint64_t m, bool periodic, bool coordinates) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateGrid2D_NM(n, m, periodic, coordinates));
    }

    PyGraph generate_grid3d(uint64_t x, uint64_t y, uint64_t z, double p, bool periodic, bool coordinates) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateGrid3D(x, y, z, p, periodic, coordinates));
    }

    PyGraph generate_grid3d_n(uint64_t n, double p, bool periodic, bool coordinates) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateGrid3D_N(n, p, periodic, coordinates));
    }

    PyGraph generate_grid3d_nm(uint64_t n, uint64_t m, bool periodic, bool coordinates) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateGrid3D_NM(n, m, periodic, coordinates));
    }

    PyGraph generate_directed_path(uint64_t n, bool permute, bool periodic) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateDirectedPath(n, permute, periodic));
    }

    PyGraph generate_kronecker(uint64_t n, uint64_t m, bool directed, bool self_loops) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateKronecker(n, m, directed, self_loops));
    }

    PyGraph generate_rmat(uint64_t n, uint64_t m, double a, double b, double c, bool directed, bool self_loops) {
        py::gil_scoped_release release;
        return PyGraph(gen.GenerateRMAT(n, m, a, b, c, directed, self_loops));
    }
};

// ============================================================================
// Module definition
// ============================================================================

PYBIND11_MODULE(_kagen, m) {
    m.doc() = "KaGen: Communication-free Graph Generators";
    m.attr("__version__") = std::to_string(KAGEN_VERSION_MAJOR) + "." +
                            std::to_string(KAGEN_VERSION_MINOR) + "." +
                            std::to_string(KAGEN_VERSION_PATCH);

#ifdef KAGEN_CGAL_FOUND
    m.attr("has_cgal") = true;
#else
    m.attr("has_cgal") = false;
#endif

    // ---- Graph class ----
    py::class_<PyGraph>(m, "Graph")
        .def("vertex_range", &PyGraph::vertex_range,
             "Returns (first_vertex, first_invalid_vertex) tuple.")
        .def("num_vertices", &PyGraph::num_vertices,
             "Number of local vertices.")
        .def("num_edges", &PyGraph::num_edges,
             "Number of local edges.")
        .def("edges", &PyGraph::edges,
             "Edge list as numpy array of shape (num_edges, 2), dtype=uint64.")
        .def("xadj", &PyGraph::xadj,
             "CSR row offsets as numpy array, dtype=uint64.")
        .def("adjncy", &PyGraph::adjncy,
             "CSR column indices as numpy array, dtype=uint64.")
        .def("vertex_weights", &PyGraph::vertex_weights,
             "Vertex weights as numpy array, dtype=int64.")
        .def("edge_weights", &PyGraph::edge_weights,
             "Edge weights as numpy array, dtype=int64.")
        .def("coordinates_2d", &PyGraph::coordinates_2d,
             "2D coordinates as numpy array of shape (N, 2), dtype=float64.")
        .def("coordinates_3d", &PyGraph::coordinates_3d,
             "3D coordinates as numpy array of shape (N, 3), dtype=float64.")
        .def("sort_edgelist", &PyGraph::sort_edgelist,
             "Sort the edge list in-place.");

    // ---- KaGen class ----
    py::class_<PyKaGen>(m, "KaGen")
        .def(py::init<>())

        // Configuration
        .def("set_seed", &PyKaGen::set_seed, py::arg("seed"),
             "Set the PRNG seed (must be same on all PEs).")
        .def("enable_undirected_graph_verification", &PyKaGen::enable_undirected_graph_verification)
        .def("enable_basic_statistics", &PyKaGen::enable_basic_statistics)
        .def("enable_advanced_statistics", &PyKaGen::enable_advanced_statistics)
        .def("enable_vertex_permutation", &PyKaGen::enable_vertex_permutation)
        .def("configure_edge_weight_generation", &PyKaGen::configure_edge_weight_generation,
             py::arg("generator"), py::arg("weight_range_begin"), py::arg("weight_range_end"),
             "Configure edge weight generation. Generator: 'hashing_based', 'uniform_random', etc.")
        .def("configure_vertex_weight_generation", &PyKaGen::configure_vertex_weight_generation,
             py::arg("generator"), py::arg("weight_range_begin"), py::arg("weight_range_end"))
        .def("enable_output", &PyKaGen::enable_output, py::arg("header") = false)
        .def("use_hp_floats", &PyKaGen::use_hp_floats, py::arg("state"))
        .def("set_num_threads", &PyKaGen::set_num_threads, py::arg("threads"),
             "Set the number of threads per MPI process for intra-node parallelism.")
        .def("set_number_of_chunks", &PyKaGen::set_number_of_chunks, py::arg("k"))
        .def("use_edge_list_representation", &PyKaGen::use_edge_list_representation)
        .def("use_csr_representation", &PyKaGen::use_csr_representation)

        // Option string
        .def("generate_from_option_string", &PyKaGen::generate_from_option_string,
             py::arg("options"))

        // Erdos-Renyi
        .def("generate_directed_gnm", &PyKaGen::generate_directed_gnm,
             py::arg("n"), py::arg("m"), py::arg("self_loops") = false)
        .def("generate_undirected_gnm", &PyKaGen::generate_undirected_gnm,
             py::arg("n"), py::arg("m"), py::arg("self_loops") = false)
        .def("generate_directed_gnp", &PyKaGen::generate_directed_gnp,
             py::arg("n"), py::arg("p"), py::arg("self_loops") = false)
        .def("generate_undirected_gnp", &PyKaGen::generate_undirected_gnp,
             py::arg("n"), py::arg("p"), py::arg("self_loops") = false)

        // RGG
        .def("generate_rgg2d", &PyKaGen::generate_rgg2d,
             py::arg("n"), py::arg("r"), py::arg("coordinates") = false)
        .def("generate_rgg2d_nm", &PyKaGen::generate_rgg2d_nm,
             py::arg("n"), py::arg("m"), py::arg("coordinates") = false)
        .def("generate_rgg2d_mr", &PyKaGen::generate_rgg2d_mr,
             py::arg("m"), py::arg("r"), py::arg("coordinates") = false)
        .def("generate_rgg3d", &PyKaGen::generate_rgg3d,
             py::arg("n"), py::arg("r"), py::arg("coordinates") = false)
        .def("generate_rgg3d_nm", &PyKaGen::generate_rgg3d_nm,
             py::arg("n"), py::arg("m"), py::arg("coordinates") = false)
        .def("generate_rgg3d_mr", &PyKaGen::generate_rgg3d_mr,
             py::arg("m"), py::arg("r"), py::arg("coordinates") = false)

        // RDG (requires CGAL)
        .def("generate_rdg2d", &PyKaGen::generate_rdg2d,
             py::arg("n"), py::arg("periodic") = false, py::arg("coordinates") = false)
        .def("generate_rdg2d_m", &PyKaGen::generate_rdg2d_m,
             py::arg("m"), py::arg("periodic") = false, py::arg("coordinates") = false)
        .def("generate_rdg3d", &PyKaGen::generate_rdg3d,
             py::arg("n"), py::arg("coordinates") = false)
        .def("generate_rdg3d_m", &PyKaGen::generate_rdg3d_m,
             py::arg("m"), py::arg("coordinates") = false)

        // BA
        .def("generate_ba", &PyKaGen::generate_ba,
             py::arg("n"), py::arg("d"), py::arg("directed") = false, py::arg("self_loops") = false)
        .def("generate_ba_nm", &PyKaGen::generate_ba_nm,
             py::arg("n"), py::arg("m"), py::arg("directed") = false, py::arg("self_loops") = false)
        .def("generate_ba_md", &PyKaGen::generate_ba_md,
             py::arg("m"), py::arg("d"), py::arg("directed") = false, py::arg("self_loops") = false)

        // RHG
        .def("generate_rhg", &PyKaGen::generate_rhg,
             py::arg("gamma"), py::arg("n"), py::arg("d"), py::arg("coordinates") = false)
        .def("generate_rhg_nm", &PyKaGen::generate_rhg_nm,
             py::arg("gamma"), py::arg("n"), py::arg("m"), py::arg("coordinates") = false)
        .def("generate_rhg_md", &PyKaGen::generate_rhg_md,
             py::arg("gamma"), py::arg("m"), py::arg("d"), py::arg("coordinates") = false)

        // Grid
        .def("generate_grid2d", &PyKaGen::generate_grid2d,
             py::arg("grid_x"), py::arg("grid_y"), py::arg("p"),
             py::arg("periodic") = false, py::arg("coordinates") = false)
        .def("generate_grid2d_n", &PyKaGen::generate_grid2d_n,
             py::arg("n"), py::arg("p"), py::arg("periodic") = false, py::arg("coordinates") = false)
        .def("generate_grid2d_nm", &PyKaGen::generate_grid2d_nm,
             py::arg("n"), py::arg("m"), py::arg("periodic") = false, py::arg("coordinates") = false)
        .def("generate_grid3d", &PyKaGen::generate_grid3d,
             py::arg("grid_x"), py::arg("grid_y"), py::arg("grid_z"), py::arg("p"),
             py::arg("periodic") = false, py::arg("coordinates") = false)
        .def("generate_grid3d_n", &PyKaGen::generate_grid3d_n,
             py::arg("n"), py::arg("p"), py::arg("periodic") = false, py::arg("coordinates") = false)
        .def("generate_grid3d_nm", &PyKaGen::generate_grid3d_nm,
             py::arg("n"), py::arg("m"), py::arg("periodic") = false, py::arg("coordinates") = false)

        // Path
        .def("generate_directed_path", &PyKaGen::generate_directed_path,
             py::arg("n"), py::arg("permute") = false, py::arg("periodic") = false)

        // Kronecker / RMAT
        .def("generate_kronecker", &PyKaGen::generate_kronecker,
             py::arg("n"), py::arg("m"), py::arg("directed") = false, py::arg("self_loops") = false)
        .def("generate_rmat", &PyKaGen::generate_rmat,
             py::arg("n"), py::arg("m"), py::arg("a"), py::arg("b"), py::arg("c"),
             py::arg("directed") = false, py::arg("self_loops") = false);

    // Initialize pseudo-MPI once at module load
    MPI_Init(nullptr, nullptr);
}
