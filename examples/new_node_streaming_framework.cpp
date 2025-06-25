#include <kagen.h>
#include <mpi.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sys/resource.h>

long getMaxRSS() {
    struct rusage usage;

    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        return usage.ru_maxrss; // in kilobytes
    } else {
      std::cout << "Error getting resource usage information." << std::endl;
    }
    return -1; // error
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 4) {
        if (rank == 0) {
            std::cout << "Usage: ./streaming_example <graph> <chunks = 32> <numNodes>" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    const std::string graph   = argv[1];
    const kagen::PEID chunks  = std::atoi(argv[2]);
    const kagen::SInt numNodes = std::atoi(argv[3]);

    if (rank == 0) {
        std::cout << "Graph: " << graph << ", chunks: " << chunks << std::endl;
    }

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

    long unsigned int nrOfEdges = 0;
    //std::ofstream outFile("streamOut.txt");
    //if (!outFile) {
    //  std::cout << "Error: Could not open file for writing." << std::endl; 
    //  return 1; 
    //} 
    //outFile << numNodes << std::endl; 

    gen.StreamNodes([&](kagen::SInt u, const std::vector<kagen::SInt>& neighbors) {
        std::cout << u << ":";
        for (kagen::SInt v : neighbors) {
            std::cout << v << " ";
            nrOfEdges++;
        }
        std::cout << "\n";
    }, kagen::StreamingMode::ordered);

    
    if (rank == 0) {
        std::cout << std::endl;
        std::cout << "Waiting for other PEs ..." << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "Number of Edges: " << nrOfEdges << std::endl; 
    long maxRSS = getMaxRSS();
    if (rank == 0) {
        std::cout << "Max RSS: " << maxRSS << " KB" << std::endl;
    }
    MPI_Finalize();
}