#include <ckagen.h>
#include <stdio.h>

#include <mpi.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    kagen_t*        gen    = kagen_create(MPI_COMM_WORLD);
    kagen_result_t* result = kagen_generate_rgg2d_nm(gen, 16, 32);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    unsigned long long from, to;
    kagen_result_vertex_range(result, &from, &to);
    printf("On PE %d [%lld, %lld): ", rank, from, to);
    size_t        nedges;
    kagen_edge_t* current_edge = kagen_result_edge_list(result, &nedges);
    for (unsigned i = 0; i < nedges; i++) {
        printf("%lld->%lld ", current_edge->source, current_edge->target);
        current_edge++;
    }
    printf("\n");

    kagen_result_free(result);
    kagen_free(gen);
    return MPI_Finalize();
}
