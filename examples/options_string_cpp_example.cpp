#include <iostream>
#include <sstream>

#include <kagen.h>

#include <mpi.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    kagen::KaGen generator(MPI_COMM_WORLD);

    // Generate a RGG2D graph with 16 nodes and 32 edges
    const kagen::Graph        graph = generator.GenerateFromOptionString("rgg2d;n=16;m=32");
    const kagen::Edgelist&    edges = graph.edges;
    const kagen::VertexRange& range = graph.vertex_range;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    {
        std::stringstream ss;
        ss << "On PE " << rank << " [" << range.first << ", " << range.second << "): ";
        for (const auto& [u, v]: edges) {
            ss << u << "->" << v << " ";
        }
        ss << "\n";
        std::cout << ss.str();
    }

    MPI_Finalize();
}
