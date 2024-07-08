#include <kagen.h>
#include <mpi.h>

#include <fstream>
#include <iostream>
#include <numeric>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        if (rank == 0) {
            std::cout << "Usage: ./streaming_example <graph> <chunks = 32>" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    const std::string graph  = argv[1];
    const kagen::PEID chunks = (argc < 3 ? 32 : std::atoi(argv[2]));

    if (rank == 0) {
        std::cout << "Graph: " << graph << ", chunks: " << chunks << std::endl;
    }

    kagen::sKaGen gen(graph, chunks, MPI_COMM_WORLD);
    gen.Initialize();

    if (rank == 0) {
        std::cout << "Generating " << std::flush;
    }

    std::vector<kagen::SInt> local_edges;
    while (gen.Continue()) {
        const kagen::StreamedGraph graph = gen.Next();

        graph.ForEachEdge([&local_edges](const auto from, const auto to) {
            local_edges.push_back(from);
            local_edges.push_back(to);
        });

        if (rank == 0) {
            std::cout << "." << std::flush;
        }
    }

    if (rank == 0) {
        std::cout << std::endl;
        std::cout << "Waiting for other PEs ..." << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "Dumping edges from root ..." << std::endl;
    }

    const int num_local_edges = static_cast<int>(local_edges.size());

    std::vector<int> recvcounts(size);
    std::vector<int> recvdispls(size);

    MPI_Gather(&num_local_edges, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::exclusive_scan(recvcounts.begin(), recvcounts.end(), recvdispls.begin(), 0);

    std::vector<kagen::SInt> global_edges(recvdispls.back() + recvcounts.back());
    MPI_Gatherv(
        local_edges.data(), num_local_edges, KAGEN_MPI_SINT, global_edges.data(), recvcounts.data(), recvdispls.data(),
        KAGEN_MPI_SINT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::ofstream out("dumped.plain-edgelist", std::ios_base::trunc);
        for (std::size_t i = 0; i < global_edges.size(); i += 2) {
            out << global_edges[i] << " " << global_edges[i + 1] << "\n";
        }
    }

    MPI_Finalize();
}
