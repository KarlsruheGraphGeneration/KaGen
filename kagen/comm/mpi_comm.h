#pragma once

#ifndef KAGEN_NOMPI

#include "kagen/comm/comm.h"

#include <mpi.h>

namespace kagen {
class MPIComm : public Comm {
public:
    explicit MPIComm(MPI_Comm comm);

    CommPEID Rank() const override;
    CommPEID Size() const override;
    void     Barrier() override;
    double   Wtime() const override;
    void     Abort(int errorcode) override;

    void Allreduce(const void* sendbuf, void* recvbuf, int count, CommDatatype datatype, CommOp op) override;

    void Reduce(
        const void* sendbuf, void* recvbuf, int count, CommDatatype datatype, CommOp op, CommPEID root) override;

    void Bcast(void* buffer, int count, CommDatatype datatype, CommPEID root) override;

    void Allgather(
        const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, int recvcount,
        CommDatatype recvtype) override;

    void Allgatherv(
        const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, const int* recvcounts,
        const int* displs, CommDatatype recvtype) override;

    void Alltoall(
        const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, int recvcount,
        CommDatatype recvtype) override;

    void Alltoallv(
        const void* sendbuf, const int* sendcounts, const int* sdispls, CommDatatype sendtype, void* recvbuf,
        const int* recvcounts, const int* rdispls, CommDatatype recvtype) override;

    void Exscan(const void* sendbuf, void* recvbuf, int count, CommDatatype datatype, CommOp op) override;

    void Gather(
        const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, int recvcount, CommDatatype recvtype,
        CommPEID root) override;

    // Expose underlying MPI communicator for code that still needs it (IO, etc.)
    MPI_Comm GetMPIComm() const;

private:
    MPI_Comm mpi_comm_;
};
} // namespace kagen

#endif // KAGEN_NOMPI
