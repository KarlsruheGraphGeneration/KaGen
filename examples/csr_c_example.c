#include <kagen.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Generate 2D RGG graph with 16 nodes and edge radius 0.125
    kagen_obj* gen = kagen_create(MPI_COMM_WORLD);
    kagen_use_csr_representation(gen);

    kagen_graph*       graph = kagen_generate_rgg2d(gen, 16, 0.125);
    unsigned long long from, to;
    kagen_graph_vertex_range(graph, &from, &to);
    printf("Vertices on PE %d: [%lld, %lld)\n", rank, from, to);

    // Get vertex distribution
    kagen_index* vtxdist = (kagen_index*)malloc((size + 1) * sizeof(kagen_index));
    kagen_build_vertex_distribution(graph, vtxdist, MPI_COMM_WORLD);
    printf("[PE%d] Vertex distribution: ", rank);
    for (int i = 0; i < size + 1; i++) {
        printf("%lld ", vtxdist[i]);
    }
    printf("\n");
    free(vtxdist);

    // Get CSR format
    size_t local_num_nodes;
    size_t local_num_edges;

    kagen_index* xadj   = kagen_graph_csr_xadj(graph, &local_num_nodes);
    kagen_index* adjncy = kagen_graph_csr_adjncy(graph, &local_num_edges);

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

    kagen_graph_free(graph);
    kagen_free(gen);
    return MPI_Finalize();
}
