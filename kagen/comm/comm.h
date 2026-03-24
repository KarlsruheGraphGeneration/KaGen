#pragma once

#include "kagen/comm/comm_types.h"

namespace kagen {
// Match kagen.h's PEID = int
using CommPEID = int;

class Comm {
public:
    virtual ~Comm() = default;

    virtual CommPEID Rank() const         = 0;
    virtual CommPEID Size() const         = 0;
    virtual void     Barrier()            = 0;
    virtual double   Wtime() const        = 0;
    virtual void     Abort(int errorcode) = 0;

    virtual void Allreduce(
        const void* sendbuf, void* recvbuf, int count, CommDatatype datatype, CommOp op) = 0;

    virtual void
    Reduce(const void* sendbuf, void* recvbuf, int count, CommDatatype datatype, CommOp op, CommPEID root) = 0;

    virtual void Bcast(void* buffer, int count, CommDatatype datatype, CommPEID root) = 0;

    virtual void Allgather(
        const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, int recvcount,
        CommDatatype recvtype) = 0;

    virtual void Allgatherv(
        const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, const int* recvcounts,
        const int* displs, CommDatatype recvtype) = 0;

    virtual void Alltoall(
        const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, int recvcount,
        CommDatatype recvtype) = 0;

    virtual void Alltoallv(
        const void* sendbuf, const int* sendcounts, const int* sdispls, CommDatatype sendtype, void* recvbuf,
        const int* recvcounts, const int* rdispls, CommDatatype recvtype) = 0;

    virtual void Exscan(const void* sendbuf, void* recvbuf, int count, CommDatatype datatype, CommOp op) = 0;

    virtual void Gather(
        const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, int recvcount, CommDatatype recvtype,
        CommPEID root) = 0;
};
} // namespace kagen
