#pragma once

#include <mpi.h>
#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef unsigned long long kagen_index;
typedef long long          kagen_weight;

typedef struct kagen_obj    kagen_obj;
typedef struct kagen_result kagen_result;

struct kagen_edge {
    kagen_index source;
    kagen_index target;
};
typedef struct kagen_edge kagen_edge;

kagen_obj* kagen_create(MPI_Comm comm);
void       kagen_free(kagen_obj* gen);

void          kagen_result_vertex_range(kagen_result* result, kagen_index* begin, kagen_index* end);
kagen_edge*   kagen_result_edge_list(kagen_result* result, size_t* nedges);
kagen_index*  kagen_result_csr_xadj(kagen_result* result, size_t* nvertices);
kagen_index*  kagen_result_csr_adjncy(kagen_result* result, size_t* nedges);
kagen_weight* kagen_result_vertex_weights(kagen_result* result, size_t* size);
kagen_weight* kagen_result_edge_weights(kagen_result* result, size_t* size);
void          kagen_result_free(kagen_result* result);

void kagen_set_seed(kagen_obj* gen, int seed);
void kagen_enable_undirected_graph_verification(kagen_obj* gen);
void kagen_enable_basic_statistics(kagen_obj* gen);
void kagen_enable_advanced_statistics(kagen_obj* gen);
void kagen_enable_output(kagen_obj* gen, bool header);
void kagen_use_hp_floats(kagen_obj* gen, bool state);
void kagen_set_numer_of_chunks(kagen_obj* gen, unsigned long long k);
void kagen_use_edge_list_representation(kagen_obj* gen);
void kagen_use_csr_representation(kagen_obj* gen);

kagen_result* kagen_generate_from_option_string(kagen_obj* gen, const char* options);

kagen_result* kagen_generate_directed_gnm(kagen_obj* gen, unsigned long long n, unsigned long long m, bool self_loops);
kagen_result*
kagen_generate_undirected_gnm(kagen_obj* gen, unsigned long long n, unsigned long long m, bool self_loops);
kagen_result* kagen_generate_directed_gnp(kagen_obj* gen, unsigned long long n, double p, bool self_loops);
kagen_result* kagen_generate_undirected_gnp(kagen_obj* gen, unsigned long long n, double p, bool self_loops);

kagen_result* kagen_generate_rgg2d(kagen_obj* gen, unsigned long long n, double r);
kagen_result* kagen_generate_rgg2d_nm(kagen_obj* gen, unsigned long long n, unsigned long long m);
kagen_result* kagen_generate_rgg2d_mr(kagen_obj* gen, unsigned long long m, double r);

kagen_result* kagen_generate_rgg3d(kagen_obj* gen, unsigned long long n, double r);
kagen_result* kagen_generate_rgg3d_nm(kagen_obj* gen, unsigned long long n, unsigned long long m);
kagen_result* kagen_generate_rgg3d_mr(kagen_obj* gen, unsigned long long m, double r);

kagen_result* kagen_generate_rdg2d(kagen_obj* gen, unsigned long long n, bool periodic);
kagen_result* kagen_generate_rdg2d_m(kagen_obj* gen, unsigned long long m, bool periodic);
kagen_result* kagen_generate_rdg3d(kagen_obj* gen, unsigned long long n);
kagen_result* kagen_generate_rdg3d_m(kagen_obj* gen, unsigned long long m);

kagen_result*
kagen_generate_ba(kagen_obj* gen, unsigned long long n, unsigned long long d, bool directed, bool self_loops);
kagen_result*
kagen_generate_ba_nm(kagen_obj* gen, unsigned long long n, unsigned long long m, bool directed, bool self_loops);
kagen_result*
kagen_generate_ba_md(kagen_obj* gen, unsigned long long m, unsigned long long d, bool directed, bool self_loops);

kagen_result* kagen_generate_rhg(kagen_obj* gen, double gamma, unsigned long long n, double d);
kagen_result* kagen_generate_rhg_nm(kagen_obj* gen, double gamma, unsigned long long n, unsigned long long m);
kagen_result* kagen_generate_rhg_md(kagen_obj* gen, double gamma, unsigned long long m, double d);

kagen_result*
kagen_generate_grid2d(kagen_obj* gen, unsigned long long grid_x, unsigned long long grid_y, double p, bool periodic);
kagen_result* kagen_generate_grid2d_n(kagen_obj* gen, unsigned long long n, double p, bool periodic);
kagen_result* kagen_generate_grid3d(
    kagen_obj* gen, unsigned long long grid_x, unsigned long long grid_y, unsigned long long grid_z, double p,
    bool periodic);
kagen_result* kagen_generate_grid3d_n(kagen_obj* gen, unsigned long long n, double p, bool periodic);

kagen_result*
kagen_generate_kronecker(kagen_obj* gen, unsigned long long n, unsigned long long m, bool directed, bool self_loops);

kagen_result* kagen_generate_rmat(
    kagen_obj* gen, unsigned long long n, unsigned long long m, double a, double b, double c, bool directed,
    bool self_loops);

void kagen_build_vertex_distribution(kagen_result* result, kagen_index* dist, MPI_Comm comm);

#ifdef __cplusplus
}
#endif
