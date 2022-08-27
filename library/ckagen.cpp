#include "ckagen.h"
#include <algorithm>
#include <mpi.h>
#include "kagen.h"

struct kagen_gen {
    void* gen_ptr;
};

struct kagen_result {
    void* result_ptr;
};

kagen_gen_t* kagen_create(MPI_Comm comm) {
    kagen_gen_t*  gen_ptr        = (kagen_gen_t*)malloc(sizeof(kagen_gen_t));
    kagen::KaGen* kagen_instance = new kagen::KaGen(comm);
    gen_ptr->gen_ptr             = kagen_instance;
    return gen_ptr;
}

void kagen_free(kagen_gen_t* gen) {
    if (gen == NULL) {
        return;
    }
    delete static_cast<kagen::KaGen*>(gen->gen_ptr);
    free(gen);
}

void kagen_result_vertex_range(kagen_result_t* result, unsigned long long* begin, unsigned long long* end) {
    kagen::KaGenResult* result_instance = static_cast<kagen::KaGenResult*>(result->result_ptr);
    *begin                              = result_instance->vertex_range.first;
    *end                                = result_instance->vertex_range.second;
}

kagen_edge_t* kagen_result_edge_list(kagen_result_t* result, size_t* nedges) {
    kagen::KaGenResult* result_instance = static_cast<kagen::KaGenResult*>(result->result_ptr);
    *nedges                             = result_instance->edges.size();
    return (kagen_edge_t*)(result_instance->edges.data());
}

void kagen_result_free(kagen_result_t* result) {
    if (result == NULL) {
        return;
    }
    delete static_cast<kagen::KaGenResult*>(result->result_ptr);
    free(result);
}

void kagen_set_seed(kagen_gen_t* gen, int seed) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    gen_instance->SetSeed(seed);
}

void kagen_enable_undirected_graph_verification(kagen_gen_t* gen) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    gen_instance->EnableUndirectedGraphVerification();
}

void kagen_enable_basic_statistics(kagen_gen_t* gen) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    gen_instance->EnableBasicStatistics();
}

void kagen_enable_advanced_statistics(kagen_gen_t* gen) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    gen_instance->EnableAdvancedStatistics();
}

void kagen_enable_output(kagen_gen_t* gen, bool header) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    gen_instance->EnableOutput(header);
}

void kagen_use_hpf_floats(kagen_gen_t* gen, bool state) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    gen_instance->UseHPFloats(state);
}
void kagen_set_numer_of_chunks(kagen_gen_t* gen, unsigned long long k) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    gen_instance->SetNumberOfChunks(k);
}

kagen_result_t*
kagen_generate_directed_gnm(kagen_gen_t* gen, unsigned long long n, unsigned long long m, bool self_loops) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateDirectedGNM(n, m, self_loops);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t*
kagen_generate_undirected_gnm(kagen_gen_t* gen, unsigned long long n, unsigned long long m, bool self_loops) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateUndirectedGNM(n, m, self_loops);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t* kagen_generate_directed_gnp(kagen_gen_t* gen, unsigned long long n, double p, bool self_loops) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateDirectedGNP(n, p, self_loops);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t* kagen_generate_undirected_gnp(kagen_gen_t* gen, unsigned long long n, double p, bool self_loops) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateUndirectedGNP(n, p, self_loops);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t* kagen_generate_rgg2d(kagen_gen_t* gen, unsigned long long n, double r) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateRGG2D(n, r);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t* kagen_generate_rgg2d_nm(kagen_gen_t* gen, unsigned long long n, unsigned long long m) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateRGG2D_NM(n, m);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t* kagen_generate_rgg2d_mr(kagen_gen_t* gen, unsigned long long m, double r) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateRGG2D_MR(m, r);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t* kagen_generate_rgg3d(kagen_gen_t* gen, unsigned long long n, double r) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateRGG3D(n, r);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t* kagen_generate_rgg3d_nm(kagen_gen_t* gen, unsigned long long n, unsigned long long m) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateRGG3D_NM(n, m);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t* kagen_generate_rgg3d_mr(kagen_gen_t* gen, unsigned long long m, double r) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateRGG3D_MR(m, r);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t* kagen_generate_rdg2d(kagen_gen_t* gen, unsigned long long n, bool periodic) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateRDG2D(n, periodic);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t* kagen_generate_rdg2d_m(kagen_gen_t* gen, unsigned long long m, bool periodic) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateRDG2D_M(m, periodic);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t* kagen_generate_rdg3d(kagen_gen_t* gen, unsigned long long n) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateRDG3D(n);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t* kagen_generate_rdg3d_m(kagen_gen_t* gen, unsigned long long m) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateRDG3D_M(m);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t*
kagen_generate_ba(kagen_gen_t* gen, unsigned long long n, unsigned long long d, bool directed, bool self_loops) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateBA(n, d, directed, self_loops);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}
kagen_result_t*
kagen_generate_ba_nm(kagen_gen_t* gen, unsigned long long n, unsigned long long m, bool directed, bool self_loops) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateBA_NM(n, m, directed, self_loops);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t*
kagen_generate_ba_md(kagen_gen_t* gen, unsigned long long m, unsigned long long d, bool directed, bool self_loops) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateBA_MD(m, d, directed, self_loops);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t* kagen_generate_rhg(kagen_gen_t* gen, double gamma, unsigned long long n, double d) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateRHG(gamma, n, d);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t* kagen_generate_rhg_nm(kagen_gen_t* gen, double gamma, unsigned long long n, unsigned long long m) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateRHG_NM(gamma, n, m);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t* kagen_generate_rhg_md(kagen_gen_t* gen, double gamma, unsigned long long m, double d) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateRHG_MD(gamma, m, d);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t*
kagen_generate_grid2d(kagen_gen_t* gen, unsigned long long grid_x, unsigned long long grid_y, double p, bool periodic) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateGrid2D(grid_x, grid_y, p, periodic);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t* kagen_generate_grid2d_n(kagen_gen_t* gen, unsigned long long n, double p, bool periodic) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateGrid2D_N(n, p, periodic);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}
kagen_result_t* kagen_generate_grid3d(
    kagen_gen_t* gen, unsigned long long grid_x, unsigned long long grid_y, unsigned long long grid_z, double p,
    bool periodic) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateGrid3D(grid_x, grid_y, grid_z, p, periodic);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}
kagen_result_t* kagen_generate_grid3d_n(kagen_gen_t* gen, unsigned long long n, double p, bool periodic) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateGrid3D_N(n, p, periodic);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t*
kagen_generate_kronecker(kagen_gen_t* gen, unsigned long long n, unsigned long long m, bool directed, bool self_loops) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateKronecker(n, m, directed, self_loops);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

kagen_result_t* kagen_generate_rmat(
    kagen_gen_t* gen, unsigned long long n, unsigned long long m, double a, double b, double c, bool directed,
    bool self_loops) {
    kagen::KaGen* gen_instance = static_cast<kagen::KaGen*>(gen->gen_ptr);
    auto          result       = gen_instance->GenerateRMAT(n, m, a, b, c, directed, self_loops);

    kagen_result_t* result_wrapper = (kagen_result_t*)malloc(sizeof(kagen_result_t));
    result_wrapper->result_ptr     = new kagen::KaGenResult(std::move(result));
    return result_wrapper;
}

void kagen_build_vertex_distribution(kagen_result_t* result, unsigned long long* dist, MPI_Comm comm) {
  if (dist == NULL) {
    return;
  }

  kagen::PEID rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  dist[0] = 0;
  dist[rank + 1] = static_cast<kagen::KaGenResult*>(result->result_ptr)->vertex_range.second;
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, dist + 1, 1, MPI_UNSIGNED_LONG_LONG, comm);
}

void kagen_build_csr(kagen_result_t* result, unsigned long long* xadj, unsigned long long* adjncy) {
    if (xadj == NULL || adjncy == NULL) {
        return;
    }
    kagen::KaGenResult* result_instance = static_cast<kagen::KaGenResult*>(result->result_ptr);
    if (!std::is_sorted(result_instance->edges.begin(), result_instance->edges.end())) {
        std::sort(result_instance->edges.begin(), result_instance->edges.end());
    }

    const unsigned long long num_local_nodes =
        result_instance->vertex_range.second - result_instance->vertex_range.first;

    unsigned long long cur_vertex = 0;
    unsigned long long cur_edge   = 0;

    xadj[0] = 0;
    for (const auto& [from, to]: result_instance->edges) {
        while (from - result_instance->vertex_range.first > cur_vertex) {
            xadj[++cur_vertex] = cur_edge;
        }
        adjncy[cur_edge++] = to;
    }
    while (cur_vertex < num_local_nodes) {
        xadj[++cur_vertex] = cur_edge;
    }
}
