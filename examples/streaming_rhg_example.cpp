#include <kagen.h>
#include <mpi.h>

#include <fstream>
#include <numeric>

int main(int argc, char* argv[]) {
    const std::string graph  = "rhg;N=20;M=22";
    const kagen::PEID chunks = 32;

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    kagen::sKaGen gen(graph, chunks, MPI_COMM_WORLD);
    gen.Initialize();

    std::vector<kagen::SInt> local_edges;
    while (gen.Continue()) {
        gen.Next().ForEachEdge([&local_edges](const auto from, const auto to) {
            local_edges.push_back(from);
            local_edges.push_back(to);
        });
    }

    const int        num_local_edges = static_cast<int>(local_edges.size());
    std::vector<int> recvcounts(size);
    MPI_Gather(&num_local_edges, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<int> recvdispls(size);
    std::exclusive_scan(recvcounts.begin(), recvcounts.end(), recvdispls.begin(), 0);
    std::vector<kagen::SInt> global_edges(recvdispls.back() + recvcounts.back());
    MPI_Gatherv(
        local_edges.data(), num_local_edges, KAGEN_MPI_SINT, global_edges.data(), recvcounts.data(), recvdispls.data(),
        KAGEN_MPI_SINT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::fstream out("rhg.edgelist");
        for (std::size_t i = 0; i < global_edges.size(); i += 2) {
            out << global_edges[i] << " " << global_edges[i + 1] << "\n";
        }
    }

    MPI_Finalize();
}
