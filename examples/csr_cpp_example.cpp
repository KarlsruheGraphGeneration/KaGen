#include <iostream>
#include <sstream>

#include <kagen.h>

#include <mpi.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Generate 2D RGG graph with 16 nodes and edge radius 0.125
    kagen::KaGen gen(MPI_COMM_WORLD);
    gen.UseCSRRepresentation();
    auto graph = gen.GenerateRGG2D(16, 0.125);
    std::cout << "Vertices on PE " << rank << ": [" << graph.vertex_range.first << ", " << graph.vertex_range.second
              << ")" << std::endl;

    // Get vertex distribution
    const auto vtxdist = kagen::BuildVertexDistribution<int>(graph, MPI_INT, MPI_COMM_WORLD);
    {
        std::stringstream out;
        out << "[PE" << rank << "] Vertex distribution: ";
        for (const auto& v: vtxdist) {
            out << v << " ";
        }
        std::cout << out.str() << std::endl;
    }

    // Get CSR format
    {
        std::stringstream out;
        out << "[PE" << rank << "] Xadj: ";
        for (const auto& v: graph.xadj) {
            out << v << " ";
        }
        std::cout << out.str() << std::endl;
    }
    {
        std::stringstream out;
        out << "[PE" << rank << "] Adjncy: ";
        for (const auto& v: graph.adjncy) {
            out << v << " ";
        }
        std::cout << out.str() << std::endl;
    }

    MPI_Finalize();
}
