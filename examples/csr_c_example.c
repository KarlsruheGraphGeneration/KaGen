#include <ckagen.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Generate 2D RGG graph with 16 nodes and edge radius 0.125
    kagen_gen*         gen   = kagen_create(MPI_COMM_WORLD);
    kagen_result*      graph = kagen_generate_rgg2d(gen, 16, 0.125);
    unsigned long long from, to;
    kagen_result_vertex_range(graph, &from, &to);
    printf("Vertices on PE %d: [%lld, %lld)\n", rank, from, to);

    // Get vertex distribution
    unsigned long long* vtxdist = malloc((size + 1) * sizeof(unsigned long long));
    kagen_build_vertex_distribution(graph, vtxdist, MPI_COMM_WORLD);
    printf("[PE%d] Vertex distribution: ", rank);
    for (int i = 0; i < size + 1; i++) {
        printf("%lld ", vtxdist[i]);
    }
    printf("\n");
    free(vtxdist);

    // Get CSR format
    size_t local_num_nodes = to - from;
    size_t local_num_edges;
    kagen_result_edge_list(graph, &local_num_edges);

    unsigned long long* xadj   = malloc((local_num_nodes + 1) * sizeof(unsigned long long));
    unsigned long long* adjncy = malloc((local_num_edges) * sizeof(unsigned long long));
    kagen_build_csr(graph, xadj, adjncy);

    printf("[PE%d] Xadj: ", rank);
    for (size_t i = 0; i < local_num_nodes + 1; i++) {
        printf("%lld ", xadj[i]);
    }
    printf("\n");

    printf("[PE%d] Adjncy: ", rank);
    for (size_t i = 0; i < local_num_edges; i++) {
        printf("%lld ", adjncy[i]);
    }
    printf("\n");

    free(xadj);
    free(adjncy);
    kagen_result_free(graph);
    kagen_free(gen);
    return MPI_Finalize();
}
