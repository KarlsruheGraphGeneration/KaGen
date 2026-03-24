#include "kagen/comm/hybrid_comm.h"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <stdexcept>

namespace kagen {

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static std::size_t ByteSize(int count, CommDatatype dtype) {
    return static_cast<std::size_t>(count) * CommDatatypeSize(dtype);
}

// Element-wise apply op: dst[i] op= src[i]  for i in [0, count)
template <typename T>
static void ApplyOpTyped(void* dst, const void* src, int count, CommOp op) {
    T*       d = static_cast<T*>(dst);
    const T* s = static_cast<const T*>(src);
    switch (op) {
        case CommOp::SUM:
            for (int i = 0; i < count; ++i) d[i] += s[i];
            break;
        case CommOp::MAX:
            for (int i = 0; i < count; ++i) d[i] = std::max(d[i], s[i]);
            break;
        case CommOp::MIN:
            for (int i = 0; i < count; ++i) d[i] = std::min(d[i], s[i]);
            break;
        case CommOp::LOR:
            for (int i = 0; i < count; ++i) d[i] = static_cast<T>(d[i] || s[i]);
            break;
    }
}

static void ApplyOp(void* dst, const void* src, int count, CommDatatype dtype, CommOp op) {
    switch (dtype) {
        case CommDatatype::INT:
            ApplyOpTyped<int>(dst, src, count, op);
            break;
        case CommDatatype::UNSIGNED:
            ApplyOpTyped<unsigned>(dst, src, count, op);
            break;
        case CommDatatype::LONG_LONG:
            ApplyOpTyped<long long>(dst, src, count, op);
            break;
        case CommDatatype::UNSIGNED_LONG_LONG:
            ApplyOpTyped<unsigned long long>(dst, src, count, op);
            break;
        case CommDatatype::DOUBLE:
            ApplyOpTyped<double>(dst, src, count, op);
            break;
        case CommDatatype::LONG_DOUBLE:
            ApplyOpTyped<long double>(dst, src, count, op);
            break;
        case CommDatatype::BYTE:
            ApplyOpTyped<unsigned char>(dst, src, count, op);
            break;
        case CommDatatype::UINT64:
            ApplyOpTyped<std::uint64_t>(dst, src, count, op);
            break;
        case CommDatatype::INT64:
            ApplyOpTyped<std::int64_t>(dst, src, count, op);
            break;
        case CommDatatype::C_BOOL:
            ApplyOpTyped<bool>(dst, src, count, op);
            break;
        case CommDatatype::DERIVED:
            break;
    }
}

// ---------------------------------------------------------------------------
// HybridComm
// ---------------------------------------------------------------------------

HybridComm::HybridComm(int thread_id, Comm& real_comm, HybridCommShared& shared)
    : thread_id_(thread_id), real_comm_(real_comm), shared_(shared) {}

CommPEID HybridComm::Rank() const {
    return real_comm_.Rank() * shared_.num_threads + thread_id_;
}

CommPEID HybridComm::Size() const {
    return real_comm_.Size() * shared_.num_threads;
}

double HybridComm::Wtime() const {
    return real_comm_.Wtime();
}

void HybridComm::Abort(int errorcode) {
    real_comm_.Abort(errorcode);
}

// ---------------------------------------------------------------------------
// Barrier
// ---------------------------------------------------------------------------

void HybridComm::Barrier() {
    shared_.arrive.arrive_and_wait();
    if (thread_id_ == 0) {
        real_comm_.Barrier();
    }
    shared_.depart.arrive_and_wait();
}

// ---------------------------------------------------------------------------
// Allreduce
// All T threads reduce locally, then thread 0 does MPI Allreduce, then
// thread 0 broadcasts the result to all threads' recvbufs.
// ---------------------------------------------------------------------------

void HybridComm::Allreduce(
    const void* sendbuf, void* recvbuf, int count, CommDatatype dtype, CommOp op) {
    auto& p   = shared_.params[thread_id_];
    p.sendbuf = (sendbuf == COMM_IN_PLACE) ? recvbuf : sendbuf;
    p.recvbuf = recvbuf;
    p.count   = count;
    p.datatype = dtype;
    p.op      = op;

    shared_.arrive.arrive_and_wait();

    if (thread_id_ == 0) {
        const int    T      = shared_.num_threads;
        const std::size_t nb = ByteSize(count, dtype);
        shared_.scratch.resize(nb);

        std::memcpy(shared_.scratch.data(), shared_.params[0].sendbuf, nb);
        for (int t = 1; t < T; ++t)
            ApplyOp(shared_.scratch.data(), shared_.params[t].sendbuf, count, dtype, op);

        real_comm_.Allreduce(COMM_IN_PLACE, shared_.scratch.data(), count, dtype, op);

        for (int t = 0; t < T; ++t)
            std::memcpy(shared_.params[t].recvbuf, shared_.scratch.data(), nb);
    }

    shared_.depart.arrive_and_wait();
}

// ---------------------------------------------------------------------------
// Reduce
// ---------------------------------------------------------------------------

void HybridComm::Reduce(
    const void* sendbuf, void* recvbuf, int count, CommDatatype dtype, CommOp op, CommPEID root) {
    const int    T           = shared_.num_threads;
    const int    root_mpi    = root / T;
    const int    root_thread = root % T;

    auto& p    = shared_.params[thread_id_];
    p.sendbuf  = (sendbuf == COMM_IN_PLACE) ? recvbuf : sendbuf;
    p.recvbuf  = recvbuf;
    p.count    = count;
    p.datatype = dtype;
    p.op       = op;
    p.root     = root;

    shared_.arrive.arrive_and_wait();

    if (thread_id_ == 0) {
        const std::size_t nb = ByteSize(count, dtype);
        shared_.scratch.resize(nb * 2);
        char* local_sum  = shared_.scratch.data();
        char* mpi_recvbuf = shared_.scratch.data() + nb;

        std::memcpy(local_sum, shared_.params[0].sendbuf, nb);
        for (int t = 1; t < T; ++t)
            ApplyOp(local_sum, shared_.params[t].sendbuf, count, dtype, op);

        real_comm_.Reduce(local_sum, mpi_recvbuf, count, dtype, op, root_mpi);

        if (real_comm_.Rank() == root_mpi)
            std::memcpy(shared_.params[root_thread].recvbuf, mpi_recvbuf, nb);
    }

    shared_.depart.arrive_and_wait();
}

// ---------------------------------------------------------------------------
// Bcast
// Thread 0 copies root_thread's buffer, broadcasts via MPI, then copies
// the result to every thread's buffer.
// ---------------------------------------------------------------------------

void HybridComm::Bcast(void* buffer, int count, CommDatatype dtype, CommPEID root) {
    const int T           = shared_.num_threads;
    const int root_mpi    = root / T;
    const int root_thread = root % T;

    auto& p    = shared_.params[thread_id_];
    p.sendbuf  = buffer;
    p.recvbuf  = buffer;
    p.count    = count;
    p.datatype = dtype;
    p.root     = root;

    shared_.arrive.arrive_and_wait();

    if (thread_id_ == 0) {
        const std::size_t nb = ByteSize(count, dtype);
        shared_.scratch.resize(nb);

        if (real_comm_.Rank() == root_mpi)
            std::memcpy(shared_.scratch.data(), shared_.params[root_thread].sendbuf, nb);

        real_comm_.Bcast(shared_.scratch.data(), count, dtype, root_mpi);

        for (int t = 0; t < T; ++t)
            std::memcpy(shared_.params[t].recvbuf, shared_.scratch.data(), nb);
    }

    shared_.depart.arrive_and_wait();
}

// ---------------------------------------------------------------------------
// Allgather
// Each virtual PE contributes sendcount elements.
// Thread 0 packs all T threads' sends, calls MPI Allgather (each MPI rank
// sends T*sendcount), then copies the full result to every thread's recvbuf.
// ---------------------------------------------------------------------------

void HybridComm::Allgather(
    const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, int recvcount,
    CommDatatype recvtype) {
    const std::size_t recv_elem_size = CommDatatypeSize(recvtype);

    auto& p    = shared_.params[thread_id_];
    // COMM_IN_PLACE: contribution is already at recvbuf[Rank()*recvcount]
    p.sendbuf  = (sendbuf == COMM_IN_PLACE)
                     ? (static_cast<char*>(recvbuf) + static_cast<std::size_t>(Rank()) * recvcount * recv_elem_size)
                     : sendbuf;
    p.recvbuf   = recvbuf;
    p.sendcount = sendcount;
    p.recvcount = recvcount;
    p.sendtype  = sendtype;
    p.recvtype  = recvtype;

    shared_.arrive.arrive_and_wait();

    if (thread_id_ == 0) {
        const int    T              = shared_.num_threads;
        const int    P              = real_comm_.Size();
        const std::size_t send_elem_size = CommDatatypeSize(sendtype);
        const std::size_t local_send_bytes = static_cast<std::size_t>(T) * sendcount * send_elem_size;
        const std::size_t total_recv_bytes = static_cast<std::size_t>(P) * T * recvcount * recv_elem_size;

        shared_.scratch.resize(local_send_bytes + total_recv_bytes);
        char* mpi_send = shared_.scratch.data();
        char* mpi_recv = shared_.scratch.data() + local_send_bytes;

        for (int t = 0; t < T; ++t) {
            std::memcpy(
                mpi_send + static_cast<std::size_t>(t) * sendcount * send_elem_size,
                shared_.params[t].sendbuf,
                static_cast<std::size_t>(sendcount) * send_elem_size);
        }

        // Each MPI rank contributes T*sendcount; each receives P*T*recvcount total
        real_comm_.Allgather(mpi_send, T * sendcount, sendtype, mpi_recv, T * recvcount, recvtype);

        // Result layout: [rank0_t0, rank0_t1, ..., rank0_t(T-1), rank1_t0, ...]
        // which is exactly virtual rank order r*T+t — copy full result to every thread
        for (int t = 0; t < T; ++t)
            std::memcpy(shared_.params[t].recvbuf, mpi_recv, total_recv_bytes);
    }

    shared_.depart.arrive_and_wait();
}

// ---------------------------------------------------------------------------
// Allgatherv
// ---------------------------------------------------------------------------

void HybridComm::Allgatherv(
    const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf,
    const int* recvcounts, const int* displs, CommDatatype recvtype) {
    const std::size_t recv_elem_size = CommDatatypeSize(recvtype);

    auto& p    = shared_.params[thread_id_];
    p.sendbuf  = (sendbuf == COMM_IN_PLACE)
                     ? (static_cast<char*>(recvbuf) + static_cast<std::size_t>(displs[Rank()]) * recv_elem_size)
                     : sendbuf;
    p.recvbuf   = recvbuf;
    p.sendcount = sendcount;
    p.sendtype  = sendtype;
    p.recvtype  = recvtype;
    p.recvcounts = recvcounts;
    p.rdispls   = displs; // reuse rdispls for allgatherv displs

    shared_.arrive.arrive_and_wait();

    if (thread_id_ == 0) {
        const int T = shared_.num_threads;
        const int P = real_comm_.Size();
        const int V = P * T;
        const std::size_t send_elem_size = CommDatatypeSize(sendtype);

        // Local send total
        int local_send_count = 0;
        for (int t = 0; t < T; ++t)
            local_send_count += shared_.params[t].sendcount;

        // Total receive count
        int total_recv_count = 0;
        for (int v = 0; v < V; ++v)
            total_recv_count += recvcounts[v];

        // MPI-level recvcounts/displs: one entry per MPI rank covering its T threads
        std::vector<int> mpi_recvcounts(P), mpi_displs(P);
        for (int r = 0; r < P; ++r) {
            mpi_recvcounts[r] = 0;
            for (int t = 0; t < T; ++t)
                mpi_recvcounts[r] += recvcounts[r * T + t];
        }
        mpi_displs[0] = 0;
        for (int r = 1; r < P; ++r)
            mpi_displs[r] = mpi_displs[r - 1] + mpi_recvcounts[r - 1];

        shared_.scratch.resize(
            static_cast<std::size_t>(local_send_count) * send_elem_size +
            static_cast<std::size_t>(total_recv_count) * recv_elem_size);
        char* mpi_send = shared_.scratch.data();
        char* mpi_recv = shared_.scratch.data() + static_cast<std::size_t>(local_send_count) * send_elem_size;

        // Pack local threads in order
        {
            char* dst = mpi_send;
            for (int t = 0; t < T; ++t) {
                const std::size_t nb = static_cast<std::size_t>(shared_.params[t].sendcount) * send_elem_size;
                std::memcpy(dst, shared_.params[t].sendbuf, nb);
                dst += nb;
            }
        }

        real_comm_.Allgatherv(mpi_send, local_send_count, sendtype, mpi_recv, mpi_recvcounts.data(),
                              mpi_displs.data(), recvtype);

        // Scatter results: for each virtual PE v, copy mpi_recv data to all threads' recvbufs
        // mpi_recv layout: [rank0_t0_data, rank0_t1_data, ..., rank(P-1)_t(T-1)_data]
        // (sequential in virtual rank order within each MPI rank block)
        for (int r = 0; r < P; ++r) {
            const char* src = mpi_recv + static_cast<std::size_t>(mpi_displs[r]) * recv_elem_size;
            for (int t = 0; t < T; ++t) {
                int vp = r * T + t;
                const std::size_t nb = static_cast<std::size_t>(recvcounts[vp]) * recv_elem_size;
                for (int dest_t = 0; dest_t < T; ++dest_t) {
                    std::memcpy(
                        static_cast<char*>(shared_.params[dest_t].recvbuf) +
                            static_cast<std::size_t>(displs[vp]) * recv_elem_size,
                        src, nb);
                }
                src += nb;
            }
        }
    }

    shared_.depart.arrive_and_wait();
}

// ---------------------------------------------------------------------------
// Alltoall — implemented via Alltoallv with uniform counts
// ---------------------------------------------------------------------------

void HybridComm::Alltoall(
    const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, int recvcount,
    CommDatatype recvtype) {
    const int V = Size();
    std::vector<int> sendcounts(V, sendcount);
    std::vector<int> sdispls(V);
    std::vector<int> recvcounts_v(V, recvcount);
    std::vector<int> rdispls(V);
    for (int i = 0; i < V; ++i) {
        sdispls[i] = i * sendcount;
        rdispls[i] = i * recvcount;
    }
    Alltoallv(sendbuf, sendcounts.data(), sdispls.data(), sendtype, recvbuf, recvcounts_v.data(),
              rdispls.data(), recvtype);
}

// ---------------------------------------------------------------------------
// Alltoallv
//
// Two phases:
//   A. Intra-rank: memcpy between threads on the same MPI rank.
//   B. Inter-rank: pack per-rank buffers, call MPI Alltoallv, unpack.
//
// Packing order for data from rank R to rank R':
//   for t_s in [0, T):    // source thread on rank R
//     for t_d in [0, T):  // dest thread on rank R'
//       params[t_s].sendbuf[sdispls[R'*T+t_d]] * sendcounts[R'*T+t_d] elements
//
// Unpack at rank R receiving from rank R' uses the same iteration order.
// ---------------------------------------------------------------------------

void HybridComm::Alltoallv(
    const void* sendbuf, const int* sendcounts, const int* sdispls, CommDatatype sendtype,
    void* recvbuf, const int* recvcounts, const int* rdispls, CommDatatype recvtype) {
    auto& p     = shared_.params[thread_id_];
    p.sendbuf   = sendbuf;
    p.recvbuf   = recvbuf;
    p.sendcounts = sendcounts;
    p.sdispls   = sdispls;
    p.recvcounts = recvcounts;
    p.rdispls   = rdispls;
    p.sendtype  = sendtype;
    p.recvtype  = recvtype;

    shared_.arrive.arrive_and_wait();

    if (thread_id_ == 0) {
        const int T = shared_.num_threads;
        const int P = real_comm_.Size();
        const int R = real_comm_.Rank();

        const std::size_t send_elem_size = CommDatatypeSize(sendtype);
        const std::size_t recv_elem_size = CommDatatypeSize(recvtype);

        // Phase A: intra-rank copies (source and dest on same MPI rank R)
        for (int t_s = 0; t_s < T; ++t_s) {
            const auto& ps = shared_.params[t_s];
            for (int t_d = 0; t_d < T; ++t_d) {
                const int dst_vp = R * T + t_d; // destination virtual PE
                const int src_vp = R * T + t_s; // source virtual PE
                const auto& pd   = shared_.params[t_d];
                const int cnt    = ps.sendcounts[dst_vp];
                if (cnt > 0) {
                    std::memcpy(
                        static_cast<char*>(pd.recvbuf) +
                            static_cast<std::size_t>(pd.rdispls[src_vp]) * recv_elem_size,
                        static_cast<const char*>(ps.sendbuf) +
                            static_cast<std::size_t>(ps.sdispls[dst_vp]) * send_elem_size,
                        static_cast<std::size_t>(cnt) * send_elem_size);
                }
            }
        }

        // Phase B: inter-rank via MPI Alltoallv
        // Compute per-rank send/recv counts
        std::vector<int> mpi_sendcounts(P, 0), mpi_recvcounts(P, 0);
        for (int r = 0; r < P; ++r) {
            if (r == R) continue;
            for (int t_s = 0; t_s < T; ++t_s)
                for (int t_d = 0; t_d < T; ++t_d)
                    mpi_sendcounts[r] += shared_.params[t_s].sendcounts[r * T + t_d];
            for (int t_d = 0; t_d < T; ++t_d)
                for (int t_s = 0; t_s < T; ++t_s)
                    mpi_recvcounts[r] += shared_.params[t_d].recvcounts[r * T + t_s];
        }

        std::vector<int> mpi_sdispls(P, 0), mpi_rdispls(P, 0);
        for (int r = 1; r < P; ++r) {
            mpi_sdispls[r] = mpi_sdispls[r - 1] + mpi_sendcounts[r - 1];
            mpi_rdispls[r] = mpi_rdispls[r - 1] + mpi_recvcounts[r - 1];
        }
        const int total_send = (P > 0) ? mpi_sdispls[P - 1] + mpi_sendcounts[P - 1] : 0;
        const int total_recv = (P > 0) ? mpi_rdispls[P - 1] + mpi_recvcounts[P - 1] : 0;

        shared_.scratch.resize(
            static_cast<std::size_t>(total_send) * send_elem_size +
            static_cast<std::size_t>(total_recv) * recv_elem_size);
        char* mpi_sendbuf = shared_.scratch.data();
        char* mpi_recvbuf = shared_.scratch.data() + static_cast<std::size_t>(total_send) * send_elem_size;

        // Pack: for each remote rank r, pack in (t_s, t_d) order
        for (int r = 0; r < P; ++r) {
            if (r == R) continue;
            char* dst = mpi_sendbuf + static_cast<std::size_t>(mpi_sdispls[r]) * send_elem_size;
            for (int t_s = 0; t_s < T; ++t_s) {
                const auto& ps = shared_.params[t_s];
                for (int t_d = 0; t_d < T; ++t_d) {
                    const int dst_vp = r * T + t_d;
                    const int cnt    = ps.sendcounts[dst_vp];
                    if (cnt > 0) {
                        std::memcpy(
                            dst,
                            static_cast<const char*>(ps.sendbuf) +
                                static_cast<std::size_t>(ps.sdispls[dst_vp]) * send_elem_size,
                            static_cast<std::size_t>(cnt) * send_elem_size);
                        dst += static_cast<std::size_t>(cnt) * send_elem_size;
                    }
                }
            }
        }

        real_comm_.Alltoallv(
            mpi_sendbuf, mpi_sendcounts.data(), mpi_sdispls.data(), sendtype, mpi_recvbuf,
            mpi_recvcounts.data(), mpi_rdispls.data(), recvtype);

        // Unpack: for each remote rank r, unpack in (t_s, t_d) order
        // (same iteration order as the packing done by rank r)
        for (int r = 0; r < P; ++r) {
            if (r == R) continue;
            const char* src = mpi_recvbuf + static_cast<std::size_t>(mpi_rdispls[r]) * recv_elem_size;
            for (int t_s = 0; t_s < T; ++t_s) {   // source thread on remote rank r
                for (int t_d = 0; t_d < T; ++t_d) { // dest thread on our rank R
                    const int src_vp = r * T + t_s;
                    const auto& pd   = shared_.params[t_d];
                    const int cnt    = pd.recvcounts[src_vp];
                    if (cnt > 0) {
                        std::memcpy(
                            static_cast<char*>(pd.recvbuf) +
                                static_cast<std::size_t>(pd.rdispls[src_vp]) * recv_elem_size,
                            src,
                            static_cast<std::size_t>(cnt) * recv_elem_size);
                        src += static_cast<std::size_t>(cnt) * recv_elem_size;
                    }
                }
            }
        }
    }

    shared_.depart.arrive_and_wait();
}

// ---------------------------------------------------------------------------
// Exscan
//
// For virtual PE r*T+t:
//   result = MPI_exscan(local_sum)[r] + local_prefix_sum_of_threads_0_to_(t-1)
//
// Thread 0 computes local_sum = sum of all T threads, calls real_comm_.Exscan,
// then adds the per-thread prefix sums to produce per-thread results.
// ---------------------------------------------------------------------------

void HybridComm::Exscan(const void* sendbuf, void* recvbuf, int count, CommDatatype dtype, CommOp op) {
    auto& p    = shared_.params[thread_id_];
    p.sendbuf  = sendbuf;
    p.recvbuf  = recvbuf;
    p.count    = count;
    p.datatype = dtype;
    p.op       = op;

    shared_.arrive.arrive_and_wait();

    if (thread_id_ == 0) {
        const int T = shared_.num_threads;
        const std::size_t nb = ByteSize(count, dtype);
        shared_.scratch.resize(nb * 3);
        char* local_sum  = shared_.scratch.data();
        char* mpi_offset = shared_.scratch.data() + nb;
        char* running    = shared_.scratch.data() + nb * 2;

        // Compute local_sum = sum of all T threads
        std::memset(local_sum, 0, nb);
        for (int t = 0; t < T; ++t)
            ApplyOp(local_sum, shared_.params[t].sendbuf, count, dtype, op);

        // MPI exscan gives the sum of all preceding MPI ranks
        std::memset(mpi_offset, 0, nb);
        real_comm_.Exscan(local_sum, mpi_offset, count, dtype, op);
        // MPI rank 0's exscan result is implementation-defined; force it to zero (identity)
        if (real_comm_.Rank() == 0)
            std::memset(mpi_offset, 0, nb);

        // Compute per-thread results:
        //   thread 0: result = mpi_offset
        //   thread t: result = mpi_offset + sum(thread_0..thread_(t-1))
        std::memcpy(running, mpi_offset, nb);
        std::memcpy(shared_.params[0].recvbuf, running, nb);
        for (int t = 1; t < T; ++t) {
            ApplyOp(running, shared_.params[t - 1].sendbuf, count, dtype, op);
            std::memcpy(shared_.params[t].recvbuf, running, nb);
        }
    }

    shared_.depart.arrive_and_wait();
}

// ---------------------------------------------------------------------------
// Gather
// Thread 0 packs all T threads' sendbufs, calls MPI Gather (each MPI rank
// sends T*sendcount), then copies the full result to the root thread's recvbuf.
// ---------------------------------------------------------------------------

void HybridComm::Gather(
    const void* sendbuf, int sendcount, CommDatatype sendtype, void* recvbuf, int recvcount,
    CommDatatype recvtype, CommPEID root) {
    const int T           = shared_.num_threads;
    const int root_mpi    = root / T;
    const int root_thread = root % T;

    auto& p     = shared_.params[thread_id_];
    p.sendbuf   = (sendbuf == COMM_IN_PLACE) ? recvbuf : sendbuf;
    p.recvbuf   = recvbuf;
    p.sendcount = sendcount;
    p.recvcount = recvcount;
    p.sendtype  = sendtype;
    p.recvtype  = recvtype;
    p.root      = root;

    shared_.arrive.arrive_and_wait();

    if (thread_id_ == 0) {
        const int    P              = real_comm_.Size();
        const std::size_t send_elem_size = CommDatatypeSize(sendtype);
        const std::size_t recv_elem_size = CommDatatypeSize(recvtype);
        const std::size_t local_send_bytes = static_cast<std::size_t>(T) * sendcount * send_elem_size;
        const std::size_t total_recv_bytes = static_cast<std::size_t>(P) * T * recvcount * recv_elem_size;

        shared_.scratch.resize(local_send_bytes + total_recv_bytes);
        char* local_send  = shared_.scratch.data();
        char* mpi_recvbuf = shared_.scratch.data() + local_send_bytes;

        for (int t = 0; t < T; ++t) {
            std::memcpy(
                local_send + static_cast<std::size_t>(t) * sendcount * send_elem_size,
                shared_.params[t].sendbuf,
                static_cast<std::size_t>(sendcount) * send_elem_size);
        }

        real_comm_.Gather(local_send, T * sendcount, sendtype, mpi_recvbuf, T * recvcount, recvtype,
                          root_mpi);

        if (real_comm_.Rank() == root_mpi)
            std::memcpy(shared_.params[root_thread].recvbuf, mpi_recvbuf, total_recv_bytes);
    }

    shared_.depart.arrive_and_wait();
}

} // namespace kagen
