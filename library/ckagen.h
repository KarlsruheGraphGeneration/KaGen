#pragma once

#include <mpi.h>
#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct kagen_gen_t;
typedef struct kagen_gen_t kagen_gen;

struct kagen_result_t;
typedef struct kagen_result_t kagen_result;

struct kagen_edge_t {
    unsigned long long target, source;
};
typedef struct kagen_edge_t kagen_edge;

kagen_gen* kagen_create(MPI_Comm comm);
void       kagen_free(kagen_gen* gen);

void        kagen_result_vertex_range(kagen_result* result, unsigned long long* begin, unsigned long long* end);
kagen_edge* kagen_result_edge_list(kagen_result* result, size_t* nedges);
void        kagen_result_free(kagen_result* result);

void kagen_set_seed(kagen_gen* gen, int seed);
void kagen_enable_undirected_graph_verification(kagen_gen* gen);
void kagen_enable_basic_statistics(kagen_gen* gen);
void kagen_enable_advanced_statistics(kagen_gen* gen);
void kagen_enable_output(kagen_gen* gen, bool header);
void kagen_use_hp_floats(kagen_gen* gen, bool state);
void kagen_set_numer_of_chunks(kagen_gen* gen, unsigned long long k);

kagen_result* kagen_generate_directed_gnm(kagen_gen* gen, unsigned long long n, unsigned long long m, bool self_loops);
kagen_result*
kagen_generate_undirected_gnm(kagen_gen* gen, unsigned long long n, unsigned long long m, bool self_loops);
kagen_result* kagen_generate_directed_gnp(kagen_gen* gen, unsigned long long n, double p, bool self_loops);
kagen_result* kagen_generate_undirected_gnp(kagen_gen* gen, unsigned long long n, double p, bool self_loops);

kagen_result* kagen_generate_rgg2d(kagen_gen* gen, unsigned long long n, double r);
kagen_result* kagen_generate_rgg2d_nm(kagen_gen* gen, unsigned long long n, unsigned long long m);
kagen_result* kagen_generate_rgg2d_mr(kagen_gen* gen, unsigned long long m, double r);

kagen_result* kagen_generate_rgg3d(kagen_gen* gen, unsigned long long n, double r);
kagen_result* kagen_generate_rgg3d_nm(kagen_gen* gen, unsigned long long n, unsigned long long m);
kagen_result* kagen_generate_rgg3d_mr(kagen_gen* gen, unsigned long long m, double r);

kagen_result* kagen_generate_rdg2d(kagen_gen* gen, unsigned long long n, bool periodic);
kagen_result* kagen_generate_rdg2d_m(kagen_gen* gen, unsigned long long m, bool periodic);
kagen_result* kagen_generate_rdg3d(kagen_gen* gen, unsigned long long n);
kagen_result* kagen_generate_rdg3d_m(kagen_gen* gen, unsigned long long m);

kagen_result*
kagen_generate_ba(kagen_gen* gen, unsigned long long n, unsigned long long d, bool directed, bool self_loops);
kagen_result*
kagen_generate_ba_nm(kagen_gen* gen, unsigned long long n, unsigned long long m, bool directed, bool self_loops);
kagen_result*
kagen_generate_ba_md(kagen_gen* gen, unsigned long long m, unsigned long long d, bool directed, bool self_loops);

kagen_result* kagen_generate_rhg(kagen_gen* gen, double gamma, unsigned long long n, double d);
kagen_result* kagen_generate_rhg_nm(kagen_gen* gen, double gamma, unsigned long long n, unsigned long long m);
kagen_result* kagen_generate_rhg_md(kagen_gen* gen, double gamma, unsigned long long m, double d);

kagen_result*
kagen_generate_grid2d(kagen_gen* gen, unsigned long long grid_x, unsigned long long grid_y, double p, bool periodic);
kagen_result* kagen_generate_grid2d_n(kagen_gen* gen, unsigned long long n, double p, bool periodic);
kagen_result* kagen_generate_grid3d(
    kagen_gen* gen, unsigned long long grid_x, unsigned long long grid_y, unsigned long long grid_z, double p,
    bool periodic);
kagen_result* kagen_generate_grid3d_n(kagen_gen* gen, unsigned long long n, double p, bool periodic);

kagen_result*
kagen_generate_kronecker(kagen_gen* gen, unsigned long long n, unsigned long long m, bool directed, bool self_loops);

kagen_result* kagen_generate_rmat(
    kagen_gen* gen, unsigned long long n, unsigned long long m, double a, double b, double c, bool directed,
    bool self_loops);

void kagen_build_vertex_distribution(kagen_result* result, unsigned long long* dist, MPI_Comm comm);
void kagen_build_csr(kagen_result* result, unsigned long long* xadj, unsigned long long* adjncy);

#ifdef __cplusplus
}
#endif
