#include <kagen.h>
#include <mpi.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 3) {
        if (rank == 0) {
            std::cout << "Usage: ./streaming_example <graph> <chunks = 32>" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    const std::string graph   = argv[1];
    const kagen::PEID chunks  = std::atoi(argv[2]);

    if (rank == 0) {
        std::cout << "Graph: " << graph << ", chunks: " << chunks << std::endl;
    }
    
    const bool sequentialGeneration = true; 

    kagen::sKaGen gen(graph, chunks, MPI_COMM_WORLD);
    
    gen.Initialize();

    // this already works so gives me the expected range
    kagen::VertexRange my_expected_vertex_range = gen.EstimateVertexRange();


    if (rank == 0) {
        for (int pe = 0; pe < size; ++pe) {
            std::cout << "Vertices on PE " << std::setw(3) << pe << ": [" << gen.EstimateVertexRange(pe).first << ", "
                      << gen.EstimateVertexRange(pe).second << ")" << std::endl;
    }

        std::cout << "Generating " << std::flush;
    }

    while (gen.Continue()) {
        const kagen::StreamedGraph graph = gen.Next();

        graph.ForEachNode([&](kagen::SInt u, const std::vector<kagen::SInt>& neighbors) {
            std::cout << u << ":"; 
            for (kagen::SInt v : neighbors) {
                std::cout << " " << v; 
            }
            std::cout << std::endl;
        }, "ordered");
    }

    
    if (rank == 0) {
        std::cout << std::endl;
        std::cout << "Waiting for other PEs ..." << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
}