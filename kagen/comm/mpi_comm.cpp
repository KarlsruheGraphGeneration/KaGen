#ifndef KAGEN_NOMPI
#include "kagen/comm/mpi_comm.h"

#include <mpi.h>

#include <cstdlib>
#include <stdexcept>

namespace kagen {

namespace {
MPI_Datatype ToMPIDatatype(CommDatatype type) {
    switch (type) {
        case CommDatatype::INT:
            return MPI_INT;
        case CommDatatype::UNSIGNED:
            return MPI_UNSIGNED;
        case CommDatatype::LONG_LONG:
            return MPI_LONG_LONG;
        case CommDatatype::UNSIGNED_LONG_LONG:
            return MPI_UNSIGNED_LONG_LONG;
        case CommDatatype::DOUBLE:
            return MPI_DOUBLE;
        case CommDatatype::LONG_DOUBLE:
            return MPI_LONG_DOUBLE;
        case CommDatatype::BYTE:
            return MPI_BYTE;
        case CommDatatype::UINT64:
            return MPI_UINT64_T;
        case CommDatatype::INT64:
            return MPI_INT64_T;
        case CommDatatype::C_BOOL:
            return MPI_C_BOOL;
        case CommDatatype::DERIVED:
            throw std::logic_error("DERIVED CommDatatype cannot be converted to MPI_Datatype directly");
    }
    throw std::logic_error("unknown CommDatatype");
}

MPI_Op ToMPIOp(CommOp op) {
    switch (op) {
        case CommOp::SUM:
            return MPI_SUM;
        case CommOp::MAX:
            return MPI_MAX;
        case CommOp::MIN:
            return MPI_MIN;
        case CommOp::LOR:
            return MPI_LOR;
    }
    throw std::logic_error("unknown CommOp");
}

// Translate COMM_IN_PLACE to MPI_IN_PLACE; otherwise keep pointer as-is.
const void* ToMPISendbuf(const void* sendbuf) {
    return (sendbuf == COMM_IN_PLACE) ? MPI_IN_PLACE : sendbuf;
}
} // namespace

MPIComm::MPIComm(MPI_Comm comm) : mpi_comm_(comm) {}

CommPEID MPIComm::Rank() const {
    int rank;
    MPI_Comm_rank(mpi_comm_, &rank);
    return rank;
}

CommPEID MPIComm::Size() const {
    int size;
    MPI_Comm_size(mpi_comm_, &size);
    return size;
}

void MPIComm::Barrier() {
    MPI_Barrier(mpi_comm_);
}

double MPIComm::Wtime() const {
    return MPI_Wtime();
}

void MPIComm::Abort(int errorcode) {
    MPI_Abort(mpi_comm_, errorcode);
}

void MPIComm::Allreduce(const void* sendbuf, void* recvbuf, int count, CommDatatype datatype, CommOp op) {
    MPI_Allreduce(ToMPISendbuf(sendbuf), recvbuf, count, ToMPIDatatype(datatype), ToMPIOp(op), mpi_comm_);
}

void MPIComm::Reduce(
    const void* sendbuf, void* recvbuf, int count, CommDatatype datatype, CommOp op, CommPEID root) {
    MPI_Reduce(ToMPISendbuf(sendbuf), recvbuf, count, ToMPIDatatype(datatype), ToMPIOp(op), root, mpi_comm_);
}

void MPIComm::Bcast(void* buffer, int count, CommDatatype datatype, CommPEID root) {
    MPI_Bcast(buffer, count, ToMPIDatatype(datatype), root, mpi_comm_);
}

void MPIComm::Allgather(
    const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, int recvcount, CommDatatype recvtype) {
    if (sendbuf == COMM_IN_PLACE) {
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, recvbuf, recvcount, ToMPIDatatype(recvtype), mpi_comm_);
    } else {
        MPI_Allgather(sendbuf, sendcount, ToMPIDatatype(sendtype), recvbuf, recvcount, ToMPIDatatype(recvtype),
                      mpi_comm_);
    }
}

void MPIComm::Allgatherv(
    const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, const int* recvcounts, const int* displs,
    CommDatatype recvtype) {
    if (sendbuf == COMM_IN_PLACE) {
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, recvbuf, recvcounts, displs, ToMPIDatatype(recvtype),
                       mpi_comm_);
    } else {
        MPI_Allgatherv(sendbuf, sendcount, ToMPIDatatype(sendtype), recvbuf, recvcounts, displs,
                       ToMPIDatatype(recvtype), mpi_comm_);
    }
}

void MPIComm::Alltoall(
    const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, int recvcount, CommDatatype recvtype) {
    MPI_Alltoall(sendbuf, sendcount, ToMPIDatatype(sendtype), recvbuf, recvcount, ToMPIDatatype(recvtype), mpi_comm_);
}

void MPIComm::Alltoallv(
    const void* sendbuf, const int* sendcounts, const int* sdispls, CommDatatype sendtype, void* recvbuf,
    const int* recvcounts, const int* rdispls, CommDatatype recvtype) {
    MPI_Alltoallv(
        sendbuf, sendcounts, sdispls, ToMPIDatatype(sendtype), recvbuf, recvcounts, rdispls, ToMPIDatatype(recvtype),
        mpi_comm_);
}

void MPIComm::Exscan(const void* sendbuf, void* recvbuf, int count, CommDatatype datatype, CommOp op) {
    MPI_Exscan(sendbuf, recvbuf, count, ToMPIDatatype(datatype), ToMPIOp(op), mpi_comm_);
}

void MPIComm::Gather(
    const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, int recvcount, CommDatatype recvtype,
    CommPEID root) {
    MPI_Gather(
        sendbuf, sendcount, ToMPIDatatype(sendtype), recvbuf, recvcount, ToMPIDatatype(recvtype), root, mpi_comm_);
}

MPI_Comm MPIComm::GetMPIComm() const {
    return mpi_comm_;
}

} // namespace kagen

#endif // KAGEN_NOMPI
