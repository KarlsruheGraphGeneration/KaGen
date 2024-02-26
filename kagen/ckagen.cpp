#include <mpi.h>

#include "kagen.h"

struct kagen_obj {
    kagen::KaGen* gen_ptr;
};

struct kagen_graph {
    kagen::Graph* result_ptr;
};

kagen_obj* kagen_create(MPI_Comm comm) {
    kagen_obj* gen_ptr = new kagen_obj;
    gen_ptr->gen_ptr   = new kagen::KaGen(comm);
    return gen_ptr;
}

void kagen_free(kagen_obj* gen) {
    if (gen == nullptr) {
        return;
    }
    delete gen->gen_ptr;
    delete gen;
}

void kagen_graph_vertex_range(kagen_graph* result, kagen_index* begin, kagen_index* end) {
    *begin = result->result_ptr->vertex_range.first;
    *end   = result->result_ptr->vertex_range.second;
}

kagen_edge* kagen_graph_edge_list(kagen_graph* result, size_t* nedges) {
    if (nedges != nullptr) {
        *nedges = result->result_ptr->edges.size();
    }
    return reinterpret_cast<kagen_edge*>(result->result_ptr->edges.data());
}

kagen_index* kagen_graph_csr_xadj(kagen_graph* result, size_t* nvertices) {
    if (nvertices != nullptr) {
        *nvertices = result->result_ptr->xadj.size() - 1;
    }
    return result->result_ptr->xadj.data();
}

kagen_index* kagen_graph_csr_adjncy(kagen_graph* result, size_t* nedges) {
    if (nedges != nullptr) {
        *nedges = result->result_ptr->adjncy.size();
    }
    return result->result_ptr->adjncy.data();
}

kagen_weight* kagen_graph_vertex_weights(kagen_graph* result, size_t* size) {
    if (size != nullptr) {
        *size = result->result_ptr->vertex_weights.size();
    }
    return result->result_ptr->vertex_weights.data();
}

kagen_weight* kagen_graph_edge_weights(kagen_graph* result, size_t* size) {
    if (size != nullptr) {
        *size = result->result_ptr->edge_weights.size();
    }
    return result->result_ptr->edge_weights.data();
}

void kagen_graph_free(kagen_graph* result) {
    if (result == nullptr) {
        return;
    }
    delete result->result_ptr;
    delete result;
}

void kagen_set_seed(kagen_obj* gen, int seed) {
    gen->gen_ptr->SetSeed(seed);
}

void kagen_enable_undirected_graph_verification(kagen_obj* gen) {
    gen->gen_ptr->EnableUndirectedGraphVerification();
}

void kagen_enable_basic_statistics(kagen_obj* gen) {
    gen->gen_ptr->EnableBasicStatistics();
}

void kagen_enable_advanced_statistics(kagen_obj* gen) {
    gen->gen_ptr->EnableAdvancedStatistics();
}

void kagen_enable_output(kagen_obj* gen, bool header) {
    gen->gen_ptr->EnableOutput(header);
}

void kagen_use_hp_floats(kagen_obj* gen, bool state) {
    gen->gen_ptr->UseHPFloats(state);
}
void kagen_set_numer_of_chunks(kagen_obj* gen, unsigned long long k) {
    gen->gen_ptr->SetNumberOfChunks(k);
}

void kagen_use_edge_list_representation(kagen_obj* gen) {
    gen->gen_ptr->UseEdgeListRepresentation();
}

void kagen_use_csr_representation(kagen_obj* gen) {
    gen->gen_ptr->UseCSRRepresentation();
}

kagen_graph* kagen_generate_from_option_string(kagen_obj* gen, const char* options) {
    auto result = gen->gen_ptr->GenerateFromOptionString(options);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph* kagen_generate_directed_gnm(kagen_obj* gen, unsigned long long n, unsigned long long m, bool self_loops) {
    auto result = gen->gen_ptr->GenerateDirectedGNM(n, m, self_loops);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph*
kagen_generate_undirected_gnm(kagen_obj* gen, unsigned long long n, unsigned long long m, bool self_loops) {
    auto result = gen->gen_ptr->GenerateUndirectedGNM(n, m, self_loops);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph* kagen_generate_directed_gnp(kagen_obj* gen, unsigned long long n, double p, bool self_loops) {
    auto result = gen->gen_ptr->GenerateDirectedGNP(n, p, self_loops);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph* kagen_generate_undirected_gnp(kagen_obj* gen, unsigned long long n, double p, bool self_loops) {
    auto result = gen->gen_ptr->GenerateUndirectedGNP(n, p, self_loops);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph* kagen_generate_rgg2d(kagen_obj* gen, unsigned long long n, double r) {
    auto result = gen->gen_ptr->GenerateRGG2D(n, r);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph* kagen_generate_rgg2d_nm(kagen_obj* gen, unsigned long long n, unsigned long long m) {
    auto result = gen->gen_ptr->GenerateRGG2D_NM(n, m);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph* kagen_generate_rgg2d_mr(kagen_obj* gen, unsigned long long m, double r) {
    auto result = gen->gen_ptr->GenerateRGG2D_MR(m, r);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph* kagen_generate_rgg3d(kagen_obj* gen, unsigned long long n, double r) {
    auto result = gen->gen_ptr->GenerateRGG3D(n, r);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph* kagen_generate_rgg3d_nm(kagen_obj* gen, unsigned long long n, unsigned long long m) {
    auto result = gen->gen_ptr->GenerateRGG3D_NM(n, m);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph* kagen_generate_rgg3d_mr(kagen_obj* gen, unsigned long long m, double r) {
    auto result = gen->gen_ptr->GenerateRGG3D_MR(m, r);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph* kagen_generate_rdg2d(kagen_obj* gen, unsigned long long n, bool periodic) {
    auto result = gen->gen_ptr->GenerateRDG2D(n, periodic);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph* kagen_generate_rdg2d_m(kagen_obj* gen, unsigned long long m, bool periodic) {
    auto result = gen->gen_ptr->GenerateRDG2D_M(m, periodic);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph* kagen_generate_rdg3d(kagen_obj* gen, unsigned long long n) {
    auto result = gen->gen_ptr->GenerateRDG3D(n);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph* kagen_generate_rdg3d_m(kagen_obj* gen, unsigned long long m) {
    auto result = gen->gen_ptr->GenerateRDG3D_M(m);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph*
kagen_generate_ba(kagen_obj* gen, unsigned long long n, unsigned long long d, bool directed, bool self_loops) {
    auto result = gen->gen_ptr->GenerateBA(n, d, directed, self_loops);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}
kagen_graph*
kagen_generate_ba_nm(kagen_obj* gen, unsigned long long n, unsigned long long m, bool directed, bool self_loops) {
    auto result = gen->gen_ptr->GenerateBA_NM(n, m, directed, self_loops);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph*
kagen_generate_ba_md(kagen_obj* gen, unsigned long long m, unsigned long long d, bool directed, bool self_loops) {
    auto result = gen->gen_ptr->GenerateBA_MD(m, d, directed, self_loops);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph* kagen_generate_rhg(kagen_obj* gen, double gamma, unsigned long long n, double d) {
    auto result = gen->gen_ptr->GenerateRHG(gamma, n, d);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph* kagen_generate_rhg_nm(kagen_obj* gen, double gamma, unsigned long long n, unsigned long long m) {
    auto result = gen->gen_ptr->GenerateRHG_NM(gamma, n, m);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph* kagen_generate_rhg_md(kagen_obj* gen, double gamma, unsigned long long m, double d) {
    auto result = gen->gen_ptr->GenerateRHG_MD(gamma, m, d);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph*
kagen_generate_grid2d(kagen_obj* gen, unsigned long long grid_x, unsigned long long grid_y, double p, bool periodic) {
    auto result = gen->gen_ptr->GenerateGrid2D(grid_x, grid_y, p, periodic);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph* kagen_generate_grid2d_n(kagen_obj* gen, unsigned long long n, double p, bool periodic) {
    auto result = gen->gen_ptr->GenerateGrid2D_N(n, p, periodic);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}
kagen_graph* kagen_generate_grid3d(
    kagen_obj* gen, unsigned long long grid_x, unsigned long long grid_y, unsigned long long grid_z, double p,
    bool periodic) {
    auto result = gen->gen_ptr->GenerateGrid3D(grid_x, grid_y, grid_z, p, periodic);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}
kagen_graph* kagen_generate_grid3d_n(kagen_obj* gen, unsigned long long n, double p, bool periodic) {
    auto result = gen->gen_ptr->GenerateGrid3D_N(n, p, periodic);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph* kagen_generate_directed_path(kagen_obj* gen, unsigned long long n, bool permute, bool periodic) {
    auto         result         = gen->gen_ptr->GenerateDirectedPath(n, permute, periodic);
    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph*
kagen_generate_kronecker(kagen_obj* gen, unsigned long long n, unsigned long long m, bool directed, bool self_loops) {
    auto result = gen->gen_ptr->GenerateKronecker(n, m, directed, self_loops);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

kagen_graph* kagen_generate_rmat(
    kagen_obj* gen, unsigned long long n, unsigned long long m, double a, double b, double c, bool directed,
    bool self_loops) {
    auto result = gen->gen_ptr->GenerateRMAT(n, m, a, b, c, directed, self_loops);

    kagen_graph* result_wrapper = new kagen_graph;
    result_wrapper->result_ptr  = new kagen::Graph(std::move(result));
    return result_wrapper;
}

void kagen_build_vertex_distribution(kagen_graph* result, kagen_index* dist, MPI_Comm comm) {
    if (dist == nullptr) {
        return;
    }

    kagen::PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    dist[0]        = 0;
    dist[rank + 1] = result->result_ptr->vertex_range.second;
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, dist + 1, 1, MPI_UNSIGNED_LONG_LONG, comm);
}
