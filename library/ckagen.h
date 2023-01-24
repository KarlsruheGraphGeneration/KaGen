#pragma once

#include <mpi.h>
#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef unsigned long long kagen_index_t;
typedef long long          kagen_weight_t;

typedef struct kagen_t        kagen_t;
typedef struct kagen_result_t kagen_result_t;

struct kagen_edge_t {
    kagen_index_t source;
    kagen_index_t target;
};
typedef struct kagen_edge_t kagen_edge_t;

kagen_t* kagen_create(MPI_Comm comm);
void     kagen_free(kagen_t* gen);

void            kagen_result_vertex_range(kagen_result_t* result, kagen_index_t* begin, kagen_index_t* end);
kagen_edge_t*   kagen_result_edge_list(kagen_result_t* result, size_t* nedges);
kagen_index_t*  kagen_result_csr_xadj(kagen_result_t* result, size_t* nvertices);
kagen_index_t*  kagen_result_csr_adjncy(kagen_result_t* result, size_t* nedges);
kagen_weight_t* kagen_result_vertex_weights(kagen_result_t* result, size_t* size);
kagen_weight_t* kagen_result_edge_weights(kagen_result_t* result, size_t* size);
void            kagen_result_free(kagen_result_t* result);

void kagen_set_seed(kagen_t* gen, int seed);
void kagen_enable_undirected_graph_verification(kagen_t* gen);
void kagen_enable_basic_statistics(kagen_t* gen);
void kagen_enable_advanced_statistics(kagen_t* gen);
void kagen_enable_output(kagen_t* gen, bool header);
void kagen_use_hp_floats(kagen_t* gen, bool state);
void kagen_set_numer_of_chunks(kagen_t* gen, unsigned long long k);
void kagen_use_edge_list_representation(kagen_t* gen);
void kagen_use_csr_representation(kagen_t* gen);

kagen_result_t* kagen_generate_from_option_string(kagen_t* gen, const char* options);

kagen_result_t* kagen_generate_directed_gnm(kagen_t* gen, unsigned long long n, unsigned long long m, bool self_loops);
kagen_result_t*
kagen_generate_undirected_gnm(kagen_t* gen, unsigned long long n, unsigned long long m, bool self_loops);
kagen_result_t* kagen_generate_directed_gnp(kagen_t* gen, unsigned long long n, double p, bool self_loops);
kagen_result_t* kagen_generate_undirected_gnp(kagen_t* gen, unsigned long long n, double p, bool self_loops);

kagen_result_t* kagen_generate_rgg2d(kagen_t* gen, unsigned long long n, double r);
kagen_result_t* kagen_generate_rgg2d_nm(kagen_t* gen, unsigned long long n, unsigned long long m);
kagen_result_t* kagen_generate_rgg2d_mr(kagen_t* gen, unsigned long long m, double r);

kagen_result_t* kagen_generate_rgg3d(kagen_t* gen, unsigned long long n, double r);
kagen_result_t* kagen_generate_rgg3d_nm(kagen_t* gen, unsigned long long n, unsigned long long m);
kagen_result_t* kagen_generate_rgg3d_mr(kagen_t* gen, unsigned long long m, double r);

kagen_result_t* kagen_generate_rdg2d(kagen_t* gen, unsigned long long n, bool periodic);
kagen_result_t* kagen_generate_rdg2d_m(kagen_t* gen, unsigned long long m, bool periodic);
kagen_result_t* kagen_generate_rdg3d(kagen_t* gen, unsigned long long n);
kagen_result_t* kagen_generate_rdg3d_m(kagen_t* gen, unsigned long long m);

kagen_result_t*
kagen_generate_ba(kagen_t* gen, unsigned long long n, unsigned long long d, bool directed, bool self_loops);
kagen_result_t*
kagen_generate_ba_nm(kagen_t* gen, unsigned long long n, unsigned long long m, bool directed, bool self_loops);
kagen_result_t*
kagen_generate_ba_md(kagen_t* gen, unsigned long long m, unsigned long long d, bool directed, bool self_loops);

kagen_result_t* kagen_generate_rhg(kagen_t* gen, double gamma, unsigned long long n, double d);
kagen_result_t* kagen_generate_rhg_nm(kagen_t* gen, double gamma, unsigned long long n, unsigned long long m);
kagen_result_t* kagen_generate_rhg_md(kagen_t* gen, double gamma, unsigned long long m, double d);

kagen_result_t*
kagen_generate_grid2d(kagen_t* gen, unsigned long long grid_x, unsigned long long grid_y, double p, bool periodic);
kagen_result_t* kagen_generate_grid2d_n(kagen_t* gen, unsigned long long n, double p, bool periodic);
kagen_result_t* kagen_generate_grid3d(
    kagen_t* gen, unsigned long long grid_x, unsigned long long grid_y, unsigned long long grid_z, double p,
    bool periodic);
kagen_result_t* kagen_generate_grid3d_n(kagen_t* gen, unsigned long long n, double p, bool periodic);

kagen_result_t*
kagen_generate_kronecker(kagen_t* gen, unsigned long long n, unsigned long long m, bool directed, bool self_loops);

kagen_result_t* kagen_generate_rmat(
    kagen_t* gen, unsigned long long n, unsigned long long m, double a, double b, double c, bool directed,
    bool self_loops);

void kagen_build_vertex_distribution(kagen_result_t* result, kagen_index_t* dist, MPI_Comm comm);

#ifdef __cplusplus
}
#endif
