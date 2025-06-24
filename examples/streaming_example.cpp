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

    if (argc < 4) {
        if (rank == 0) {
            std::cout << "Usage: ./streaming_example <graph> <chunks = 32> <out-dir>" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    const std::string graph   = argv[1];
    const kagen::PEID chunks  = std::atoi(argv[2]);
    const std::string out_dir = argv[3];

    if (rank == 0) {
        std::cout << "Graph: " << graph << ", chunks: " << chunks << std::endl;
    }

    kagen::sKaGen gen(graph, chunks, MPI_COMM_WORLD);
    gen.Initialize();

    kagen::VertexRange my_expected_vertex_range = gen.EstimateVertexRange();

    if (rank == 0) {
        for (int pe = 0; pe < size; ++pe) {
            std::cout << "Vertices on PE " << std::setw(3) << pe << ": [" << gen.EstimateVertexRange(pe).first << ", "
                      << gen.EstimateVertexRange(pe).second << ")" << std::endl;
        }

        std::cout << "Generating " << std::flush;
    }

    std::ofstream out(out_dir + "/" + std::to_string(rank) + ".edges", std::ios_base::trunc);

    while (gen.Continue()) {
        const kagen::StreamedGraph graph = gen.Next();

        std::vector<kagen::SInt> local_edges;
        graph.ForEachEdge([&](const auto from, const auto to) {
            local_edges.push_back(from);
            local_edges.push_back(to);
        }, kagen::StreamingMode::all);
        out.write(reinterpret_cast<const char*>(local_edges.data()), local_edges.size() * sizeof(kagen::SInt));
        local_edges.clear();

        if (rank == 0) {
            std::cout << "." << std::flush;
        }
    }

    if (rank == 0) {
        std::cout << std::endl;
        std::cout << "Waiting for other PEs ..." << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
}
