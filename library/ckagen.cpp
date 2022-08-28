#include <algorithm>

#include <mpi.h>

#include "ckagen.h"
#include "kagen.h"

struct kagen_gen_t {
    kagen::KaGen* gen_ptr;
};

struct kagen_result_t {
    kagen::KaGenResult* result_ptr;
};

kagen_gen* kagen_create(MPI_Comm comm) {
    kagen_gen_t* gen_ptr = new kagen_gen_t;
    gen_ptr->gen_ptr     = new kagen::KaGen(comm);
    return gen_ptr;
}

void kagen_free(kagen_gen* gen) {
    if (gen == NULL) {
        return;
    }
    delete gen->gen_ptr;
    delete gen;
}

void kagen_result_vertex_range(kagen_result* result, unsigned long long* begin, unsigned long long* end) {
    *begin = result->result_ptr->vertex_range.first;
    *end   = result->result_ptr->vertex_range.second;
}

kagen_edge* kagen_result_edge_list(kagen_result* result, size_t* nedges) {
    *nedges = result->result_ptr->edges.size();
    return reinterpret_cast<kagen_edge*>(result->result_ptr->edges.data());
}

void kagen_result_free(kagen_result* result) {
    if (result == NULL) {
        return;
    }
    delete result->result_ptr;
    delete result;
}

void kagen_set_seed(kagen_gen* gen, int seed) {
    gen->gen_ptr->SetSeed(seed);
}

void kagen_enable_undirected_graph_verification(kagen_gen* gen) {
    gen->gen_ptr->EnableUndirectedGraphVerification();
}

void kagen_enable_basic_statistics(kagen_gen* gen) {
    gen->gen_ptr->EnableBasicStatistics();
}

void kagen_enable_advanced_statistics(kagen_gen* gen) {
    gen->gen_ptr->EnableAdvancedStatistics();
}

void kagen_enable_output(kagen_gen* gen, bool header) {
    gen->gen_ptr->EnableOutput(header);
}

void kagen_use_hpf_floats(kagen_gen* gen, bool state) {
    gen->gen_ptr->UseHPFloats(state);
}
void kagen_set_numer_of_chunks(kagen_gen* gen, unsigned long long k) {
    gen->gen_ptr->SetNumberOfChunks(k);
}

kagen_result* kagen_generate_directed_gnm(kagen_gen* gen, unsigned long long n, unsigned long long m, bool self_loops) {
    auto result = gen->gen_ptr->GenerateDirectedGNM(n, m, self_loops);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result*
kagen_generate_undirected_gnm(kagen_gen* gen, unsigned long long n, unsigned long long m, bool self_loops) {
    auto result = gen->gen_ptr->GenerateUndirectedGNM(n, m, self_loops);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result* kagen_generate_directed_gnp(kagen_gen* gen, unsigned long long n, double p, bool self_loops) {
    auto result = gen->gen_ptr->GenerateDirectedGNP(n, p, self_loops);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result* kagen_generate_undirected_gnp(kagen_gen* gen, unsigned long long n, double p, bool self_loops) {
    auto result = gen->gen_ptr->GenerateUndirectedGNP(n, p, self_loops);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result* kagen_generate_rgg2d(kagen_gen* gen, unsigned long long n, double r) {
    auto result = gen->gen_ptr->GenerateRGG2D(n, r);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result* kagen_generate_rgg2d_nm(kagen_gen* gen, unsigned long long n, unsigned long long m) {
    auto result = gen->gen_ptr->GenerateRGG2D_NM(n, m);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result* kagen_generate_rgg2d_mr(kagen_gen* gen, unsigned long long m, double r) {
    auto result = gen->gen_ptr->GenerateRGG2D_MR(m, r);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result* kagen_generate_rgg3d(kagen_gen* gen, unsigned long long n, double r) {
    auto result = gen->gen_ptr->GenerateRGG3D(n, r);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result* kagen_generate_rgg3d_nm(kagen_gen* gen, unsigned long long n, unsigned long long m) {
    auto result = gen->gen_ptr->GenerateRGG3D_NM(n, m);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result* kagen_generate_rgg3d_mr(kagen_gen* gen, unsigned long long m, double r) {
    auto result = gen->gen_ptr->GenerateRGG3D_MR(m, r);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result* kagen_generate_rdg2d(kagen_gen* gen, unsigned long long n, bool periodic) {
    auto result = gen->gen_ptr->GenerateRDG2D(n, periodic);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result* kagen_generate_rdg2d_m(kagen_gen* gen, unsigned long long m, bool periodic) {
    auto result = gen->gen_ptr->GenerateRDG2D_M(m, periodic);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result* kagen_generate_rdg3d(kagen_gen* gen, unsigned long long n) {
    auto result = gen->gen_ptr->GenerateRDG3D(n);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result* kagen_generate_rdg3d_m(kagen_gen* gen, unsigned long long m) {
    auto result = gen->gen_ptr->GenerateRDG3D_M(m);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result*
kagen_generate_ba(kagen_gen* gen, unsigned long long n, unsigned long long d, bool directed, bool self_loops) {
    auto result = gen->gen_ptr->GenerateBA(n, d, directed, self_loops);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}
kagen_result*
kagen_generate_ba_nm(kagen_gen* gen, unsigned long long n, unsigned long long m, bool directed, bool self_loops) {
    auto result = gen->gen_ptr->GenerateBA_NM(n, m, directed, self_loops);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result*
kagen_generate_ba_md(kagen_gen* gen, unsigned long long m, unsigned long long d, bool directed, bool self_loops) {
    auto result = gen->gen_ptr->GenerateBA_MD(m, d, directed, self_loops);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result* kagen_generate_rhg(kagen_gen* gen, double gamma, unsigned long long n, double d) {
    auto result = gen->gen_ptr->GenerateRHG(gamma, n, d);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result* kagen_generate_rhg_nm(kagen_gen* gen, double gamma, unsigned long long n, unsigned long long m) {
    auto result = gen->gen_ptr->GenerateRHG_NM(gamma, n, m);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result* kagen_generate_rhg_md(kagen_gen* gen, double gamma, unsigned long long m, double d) {
    auto result = gen->gen_ptr->GenerateRHG_MD(gamma, m, d);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result*
kagen_generate_grid2d(kagen_gen* gen, unsigned long long grid_x, unsigned long long grid_y, double p, bool periodic) {
    auto result = gen->gen_ptr->GenerateGrid2D(grid_x, grid_y, p, periodic);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result* kagen_generate_grid2d_n(kagen_gen* gen, unsigned long long n, double p, bool periodic) {
    auto result = gen->gen_ptr->GenerateGrid2D_N(n, p, periodic);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}
kagen_result* kagen_generate_grid3d(
    kagen_gen* gen, unsigned long long grid_x, unsigned long long grid_y, unsigned long long grid_z, double p,
    bool periodic) {
    auto result = gen->gen_ptr->GenerateGrid3D(grid_x, grid_y, grid_z, p, periodic);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}
kagen_result* kagen_generate_grid3d_n(kagen_gen* gen, unsigned long long n, double p, bool periodic) {
    auto result = gen->gen_ptr->GenerateGrid3D_N(n, p, periodic);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result*
kagen_generate_kronecker(kagen_gen* gen, unsigned long long n, unsigned long long m, bool directed, bool self_loops) {
    auto result = gen->gen_ptr->GenerateKronecker(n, m, directed, self_loops);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result* kagen_generate_rmat(
    kagen_gen* gen, unsigned long long n, unsigned long long m, double a, double b, double c, bool directed,
    bool self_loops) {
    auto result = gen->gen_ptr->GenerateRMAT(n, m, a, b, c, directed, self_loops);

    kagen_result* result_wrapper = new kagen_result;
    result_wrapper->result_ptr   = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

void kagen_build_vertex_distribution(kagen_result* result, unsigned long long* dist, MPI_Comm comm) {
    if (dist == NULL) {
        return;
    }

    kagen::PEID rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    dist[0]        = 0;
    dist[rank + 1] = result->result_ptr->vertex_range.second;
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, dist + 1, 1, MPI_UNSIGNED_LONG_LONG, comm);
}

void kagen_build_csr(kagen_result* result, unsigned long long* xadj, unsigned long long* adjncy) {
    if (xadj == NULL || adjncy == NULL) {
        return;
    }
    if (!std::is_sorted(result->result_ptr->edges.begin(), result->result_ptr->edges.end())) {
        std::sort(result->result_ptr->edges.begin(), result->result_ptr->edges.end());
    }

    const unsigned long long num_local_nodes =
        result->result_ptr->vertex_range.second - result->result_ptr->vertex_range.first;

    unsigned long long cur_vertex = 0;
    unsigned long long cur_edge   = 0;

    xadj[0] = 0;
    for (const auto& [from, to]: result->result_ptr->edges) {
        while (from - result->result_ptr->vertex_range.first > cur_vertex) {
            xadj[++cur_vertex] = cur_edge;
        }
        adjncy[cur_edge++] = to;
    }
    while (cur_vertex < num_local_nodes) {
        xadj[++cur_vertex] = cur_edge;
    }
}
