#include "kagen/comm/seq_comm.h"

#include <chrono>
#include <cstdlib>
#include <cstring>

namespace kagen {

CommPEID SeqComm::Rank() const {
    return 0;
}

CommPEID SeqComm::Size() const {
    return 1;
}

void SeqComm::Barrier() {
    // No-op: single process
}

double SeqComm::Wtime() const {
    using Clock = std::chrono::steady_clock;
    static const auto start = Clock::now();
    return std::chrono::duration<double>(Clock::now() - start).count();
}

void SeqComm::Abort(int errorcode) {
    std::exit(errorcode);
}

// Helper: byte size of a contiguous block
static std::size_t ByteSize(int count, CommDatatype datatype) {
    return static_cast<std::size_t>(count) * CommDatatypeSize(datatype);
}

void SeqComm::Allreduce(const void* sendbuf, void* recvbuf, int count, CommDatatype datatype, CommOp /*op*/) {
    if (sendbuf == COMM_IN_PLACE) {
        // recvbuf already holds the single value — no-op
        return;
    }
    std::memcpy(recvbuf, sendbuf, ByteSize(count, datatype));
}

void SeqComm::Reduce(
    const void* sendbuf, void* recvbuf, int count, CommDatatype datatype, CommOp op, CommPEID /*root*/) {
    Allreduce(sendbuf, recvbuf, count, datatype, op);
}

void SeqComm::Bcast(void* /*buffer*/, int /*count*/, CommDatatype /*datatype*/, CommPEID /*root*/) {
    // No-op: single process, data already present
}

void SeqComm::Allgather(
    const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, int /*recvcount*/,
    CommDatatype /*recvtype*/) {
    if (sendbuf == COMM_IN_PLACE) {
        return;
    }
    std::memcpy(recvbuf, sendbuf, ByteSize(sendcount, sendtype));
}

void SeqComm::Allgatherv(
    const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, const int* /*recvcounts*/,
    const int* displs, CommDatatype recvtype) {
    if (sendbuf == COMM_IN_PLACE) {
        return;
    }
    const std::size_t recv_elem_size = CommDatatypeSize(recvtype);
    std::memcpy(
        static_cast<char*>(recvbuf) + static_cast<std::size_t>(displs[0]) * recv_elem_size, sendbuf,
        ByteSize(sendcount, sendtype));
}

void SeqComm::Alltoall(
    const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, int /*recvcount*/,
    CommDatatype /*recvtype*/) {
    if (sendbuf == COMM_IN_PLACE) {
        return;
    }
    std::memcpy(recvbuf, sendbuf, ByteSize(sendcount, sendtype));
}

void SeqComm::Alltoallv(
    const void* sendbuf, const int* sendcounts, const int* sdispls, CommDatatype sendtype, void* recvbuf,
    const int* recvcounts, const int* rdispls, CommDatatype recvtype) {
    if (sendbuf == COMM_IN_PLACE) {
        return;
    }
    const std::size_t send_elem_size = CommDatatypeSize(sendtype);
    const std::size_t recv_elem_size = CommDatatypeSize(recvtype);
    std::memcpy(
        static_cast<char*>(recvbuf) + static_cast<std::size_t>(rdispls[0]) * recv_elem_size,
        static_cast<const char*>(sendbuf) + static_cast<std::size_t>(sdispls[0]) * send_elem_size,
        static_cast<std::size_t>(sendcounts[0]) * send_elem_size);
    (void)recvcounts;
}

void SeqComm::Exscan(const void* /*sendbuf*/, void* recvbuf, int count, CommDatatype datatype, CommOp /*op*/) {
    // Rank 0: exclusive prefix sum has no preceding ranks — result is the identity element.
    // Zeroing is the identity for SUM and LOR; for MAX/MIN it is not strictly correct, but
    // the callers in KaGen only use Exscan to compute offsets (SUM), so zeroing is safe.
    std::memset(recvbuf, 0, ByteSize(count, datatype));
}

void SeqComm::Gather(
    const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, int /*recvcount*/,
    CommDatatype /*recvtype*/, CommPEID /*root*/) {
    if (sendbuf == COMM_IN_PLACE) {
        return;
    }
    std::memcpy(recvbuf, sendbuf, ByteSize(sendcount, sendtype));
}

} // namespace kagen
