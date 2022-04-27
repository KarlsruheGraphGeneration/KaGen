#include <iostream>
#include <sstream>

#include "kagen_library.h"

#include <mpi.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Generate 2D RGG graph with 16 nodes and edge radius 0.125
    auto graph = kagen::KaGen(MPI_COMM_WORLD).GenerateRGG2D(16, 0.125);
    std::cout << "Vertices on PE " << rank << ": [" << graph.vertex_range.first << ", " << graph.vertex_range.second
              << ")" << std::endl;

    // Get vertex distribution
    int* vtxdist;
    int  vtxdist_size;
    kagen::BuildVertexDistribution(graph, &vtxdist, &vtxdist_size, MPI_INT, MPI_COMM_WORLD);
    {
        std::stringstream out;
        out << "[PE" << rank << "] Vertex distribution (size=" << vtxdist_size << "): ";
        for (int i = 0; i < vtxdist_size; ++i) {
            out << vtxdist[i] << " ";
        }
        std::cout << out.str() << std::endl;
    }
    delete[] vtxdist;

    // Get CSR format
    int* xadj;
    int  xadj_size;
    int* adjncy;
    int  adjncy_size;
    kagen::BuildCSR(graph, &xadj, &xadj_size, &adjncy, &adjncy_size);
    {
        std::stringstream out;
        out << "[PE" << rank << "] Xadj (size=" << xadj_size << "): ";
        for (int i = 0; i < xadj_size; ++i) {
            out << xadj[i] << " ";
        }
        std::cout << out.str() << std::endl;
    }
    {
        std::stringstream out;
        out << "[PE" << rank << "] Adjncy (size=" << adjncy_size << "): ";
        for (int i = 0; i < adjncy_size; ++i) {
            out << adjncy[i] << " ";
        }
        std::cout << out.str() << std::endl;
    }
    delete[] xadj;
    delete[] adjncy;

    MPI_Finalize();
}
