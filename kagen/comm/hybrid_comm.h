#pragma once

#include "kagen/comm/comm.h"

#include <barrier>
#include <vector>

namespace kagen {

// Per-thread parameters posted before each collective call.
// Thread 0 reads all threads' params after the arrive barrier and dispatches work.
struct HybridThreadParams {
    const void*  sendbuf    = nullptr;
    void*        recvbuf    = nullptr;
    int          count      = 0;
    CommDatatype datatype   = CommDatatype::BYTE;
    CommOp       op         = CommOp::SUM;
    CommPEID     root       = 0;
    int          sendcount  = 0;
    int          recvcount  = 0;
    const int*   sendcounts = nullptr;
    const int*   sdispls    = nullptr;
    const int*   recvcounts = nullptr;
    const int*   rdispls    = nullptr;
    CommDatatype sendtype   = CommDatatype::BYTE;
    CommDatatype recvtype   = CommDatatype::BYTE;
};

// State shared by all T HybridComm instances on the same MPI rank.
// Create one HybridCommShared per MPI rank, then construct T HybridComm objects from it.
struct HybridCommShared {
    int                         num_threads;
    std::barrier<>              arrive;
    std::barrier<>              depart;
    std::vector<HybridThreadParams> params;
    std::vector<char>           scratch; // scratch buffer written exclusively by thread 0

    explicit HybridCommShared(int T) : num_threads(T), arrive(T), depart(T), params(T) {}

    // Non-copyable, non-movable (barrier is not movable after construction)
    HybridCommShared(const HybridCommShared&)            = delete;
    HybridCommShared& operator=(const HybridCommShared&) = delete;
};

// A Comm implementation for one thread in a hybrid MPI+threads execution.
//
// Rank() = real_comm.Rank() * num_threads + thread_id
// Size() = real_comm.Size() * num_threads
//
// Collectives coordinate T threads on the same MPI rank using two std::barrier phases:
//   1. arrive barrier: all threads have written their params
//   2. thread 0 performs the operation (local reduction + MPI collective via real_comm)
//   3. depart barrier: thread 0 has written results to all threads' recvbufs
class HybridComm : public Comm {
public:
    // thread_id: 0-based index within [0, num_threads)
    // real_comm: the underlying MPI-level Comm (one per MPI rank, shared)
    // shared:    shared state for all T threads on this rank
    HybridComm(int thread_id, Comm& real_comm, HybridCommShared& shared);

    CommPEID Rank() const override;
    CommPEID Size() const override;
    void     Barrier() override;
    double   Wtime() const override;
    void     Abort(int errorcode) override;

    void Allreduce(
        const void* sendbuf, void* recvbuf, int count, CommDatatype datatype, CommOp op) override;

    void Reduce(
        const void* sendbuf, void* recvbuf, int count, CommDatatype datatype, CommOp op,
        CommPEID root) override;

    void Bcast(void* buffer, int count, CommDatatype datatype, CommPEID root) override;

    void Allgather(
        const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, int recvcount,
        CommDatatype recvtype) override;

    void Allgatherv(
        const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf,
        const int* recvcounts, const int* displs, CommDatatype recvtype) override;

    void Alltoall(
        const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, int recvcount,
        CommDatatype recvtype) override;

    void Alltoallv(
        const void* sendbuf, const int* sendcounts, const int* sdispls, CommDatatype sendtype,
        void* recvbuf, const int* recvcounts, const int* rdispls, CommDatatype recvtype) override;

    void Exscan(const void* sendbuf, void* recvbuf, int count, CommDatatype datatype, CommOp op) override;

    void Gather(
        const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, int recvcount,
        CommDatatype recvtype, CommPEID root) override;

private:
    int               thread_id_;
    Comm&             real_comm_;
    HybridCommShared& shared_;
};

} // namespace kagen
