#include <kagen.h>
#include <stdio.h>

#include <mpi.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    kagen_obj*   gen   = kagen_create(MPI_COMM_WORLD);
    kagen_graph* graph = kagen_generate_from_option_string(gen, "rgg2d;n=16;m=32");

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    unsigned long long from, to;
    kagen_graph_vertex_range(graph, &from, &to);
    printf("On PE %d [%lld, %lld): ", rank, from, to);
    size_t      nedges;
    kagen_edge* current_edge = kagen_graph_edge_list(graph, &nedges);
    for (unsigned i = 0; i < nedges; i++) {
        printf("%lld->%lld ", current_edge->source, current_edge->target);
        current_edge++;
    }
    printf("\n");

    kagen_graph_free(graph);
    kagen_free(gen);
    return MPI_Finalize();
}
