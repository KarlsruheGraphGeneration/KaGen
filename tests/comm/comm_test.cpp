// Unit tests for the Comm hierarchy.
//
// Coverage:
//   SeqComm     – single-process, all operations
//   VirtualComm – rank/size override inheriting SeqComm semantics
//   HybridComm  – multi-thread over SeqComm (single MPI rank, T in {1,2,4,8})
//   HybridComm  – multi-thread over MPIComm  (P MPI ranks × T threads, !KAGEN_NOMPI)
//   MPIComm     – MPI wrapper basics          (!KAGEN_NOMPI)

#include "kagen/comm/hybrid_comm.h"
#include "kagen/comm/seq_comm.h"
#include "kagen/comm/virtual_comm.h"

#ifndef KAGEN_NOMPI
    #include "kagen/comm/mpi_comm.h"
    #include <mpi.h>
#endif

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <atomic>
#include <numeric>
#include <thread>
#include <vector>

using namespace kagen;
using ::testing::ElementsAreArray;

// ---------------------------------------------------------------------------
// Helper: run body(thread_id, num_threads, HybridComm&) in T std::threads.
// ---------------------------------------------------------------------------
template <typename F>
static void RunHybrid(int T, Comm& underlying, F&& body) {
    HybridCommShared shared(T);
    std::vector<std::thread> threads;
    threads.reserve(T);
    for (int t = 0; t < T; ++t) {
        threads.emplace_back([&, t]() {
            HybridComm comm(t, underlying, shared);
            body(t, T, comm);
        });
    }
    for (auto& th : threads)
        th.join();
}

template <typename F>
static void RunHybrid(int T, F&& body) {
    SeqComm seq;
    RunHybrid(T, seq, std::forward<F>(body));
}

// ===========================================================================
// SeqComm
// ===========================================================================

TEST(SeqComm, RankAndSize) {
    SeqComm comm;
    EXPECT_EQ(comm.Rank(), 0);
    EXPECT_EQ(comm.Size(), 1);
}

TEST(SeqComm, AllreduceSumCopy) {
    SeqComm comm;
    int send = 42, recv = 0;
    comm.Allreduce(&send, &recv, 1, CommDatatype::INT, CommOp::SUM);
    EXPECT_EQ(recv, 42);
}

TEST(SeqComm, AllreduceSumInPlace) {
    SeqComm comm;
    int val = 7;
    comm.Allreduce(COMM_IN_PLACE, &val, 1, CommDatatype::INT, CommOp::SUM);
    EXPECT_EQ(val, 7); // single process – no-op
}

TEST(SeqComm, AllreduceMaxArray) {
    SeqComm comm;
    std::vector<int> send = {3, 1, 4};
    std::vector<int> recv(3, 0);
    comm.Allreduce(send.data(), recv.data(), 3, CommDatatype::INT, CommOp::MAX);
    EXPECT_THAT(recv, ElementsAreArray({3, 1, 4}));
}

TEST(SeqComm, ReduceCopiesLongLong) {
    SeqComm comm;
    long long send = 100LL, recv = 0;
    comm.Reduce(&send, &recv, 1, CommDatatype::LONG_LONG, CommOp::SUM, 0);
    EXPECT_EQ(recv, 100LL);
}

TEST(SeqComm, BcastIsNoOp) {
    SeqComm comm;
    int val = 99;
    comm.Bcast(&val, 1, CommDatatype::INT, 0);
    EXPECT_EQ(val, 99);
}

TEST(SeqComm, AllgatherCopy) {
    SeqComm comm;
    int send = 5, recv = 0;
    comm.Allgather(&send, 1, CommDatatype::INT, &recv, 1, CommDatatype::INT);
    EXPECT_EQ(recv, 5);
}

TEST(SeqComm, AllgatherInPlaceIsNoOp) {
    SeqComm comm;
    int recv = 5;
    comm.Allgather(COMM_IN_PLACE, 1, CommDatatype::INT, &recv, 1, CommDatatype::INT);
    EXPECT_EQ(recv, 5);
}

TEST(SeqComm, AllgathervCopiesAtDisplacement) {
    SeqComm comm;
    std::vector<int> send = {7, 8};
    std::vector<int> recv(4, 0);
    int recvcount = 2, displ = 2;
    comm.Allgatherv(send.data(), 2, CommDatatype::INT, recv.data(), &recvcount, &displ, CommDatatype::INT);
    EXPECT_THAT(recv, ElementsAreArray({0, 0, 7, 8}));
}

TEST(SeqComm, AlltoallCopy) {
    SeqComm comm;
    std::vector<int> send = {1, 2, 3};
    std::vector<int> recv(3, 0);
    comm.Alltoall(send.data(), 3, CommDatatype::INT, recv.data(), 3, CommDatatype::INT);
    EXPECT_THAT(recv, ElementsAreArray({1, 2, 3}));
}

TEST(SeqComm, AlltoallvCopiesAtDisplacement) {
    SeqComm comm;
    std::vector<int> send = {10, 20, 30};
    std::vector<int> recv(5, 0);
    int sc = 3, sd = 0, rc = 3, rd = 2; // write to recv[2..4]
    comm.Alltoallv(send.data(), &sc, &sd, CommDatatype::INT, recv.data(), &rc, &rd, CommDatatype::INT);
    EXPECT_THAT(recv, ElementsAreArray({0, 0, 10, 20, 30}));
}

TEST(SeqComm, ExscanIsZero) {
    SeqComm comm;
    int send = 42, recv = -1;
    comm.Exscan(&send, &recv, 1, CommDatatype::INT, CommOp::SUM);
    EXPECT_EQ(recv, 0);
}

TEST(SeqComm, ExscanArrayIsAllZero) {
    SeqComm comm;
    std::vector<int> send = {1, 2, 3};
    std::vector<int> recv(3, -1);
    comm.Exscan(send.data(), recv.data(), 3, CommDatatype::INT, CommOp::SUM);
    EXPECT_THAT(recv, ElementsAreArray({0, 0, 0}));
}

TEST(SeqComm, GatherCopy) {
    SeqComm comm;
    std::vector<int> send = {1, 2};
    std::vector<int> recv(2, 0);
    comm.Gather(send.data(), 2, CommDatatype::INT, recv.data(), 2, CommDatatype::INT, 0);
    EXPECT_THAT(recv, ElementsAreArray({1, 2}));
}

// ===========================================================================
// VirtualComm
// ===========================================================================

TEST(VirtualComm, RankAndSize) {
    VirtualComm comm(3, 8);
    EXPECT_EQ(comm.Rank(), 3);
    EXPECT_EQ(comm.Size(), 8);
}

TEST(VirtualComm, InheritsSeqAllreduce) {
    VirtualComm comm(2, 4);
    int send = 10, recv = 0;
    comm.Allreduce(&send, &recv, 1, CommDatatype::INT, CommOp::SUM);
    EXPECT_EQ(recv, 10);
}

TEST(VirtualComm, InheritsSeqExscan) {
    VirtualComm comm(2, 4);
    int send = 5, recv = -1;
    comm.Exscan(&send, &recv, 1, CommDatatype::INT, CommOp::SUM);
    EXPECT_EQ(recv, 0); // single-process; always zero
}

// ===========================================================================
// HybridComm over SeqComm  (single MPI rank, T threads)
// ===========================================================================

class HybridSeqTest : public ::testing::TestWithParam<int> {};

INSTANTIATE_TEST_SUITE_P(
    ThreadCounts, HybridSeqTest,
    ::testing::Values(1, 2, 4, 8),
    [](const ::testing::TestParamInfo<int>& info) {
        return std::to_string(info.param) + "threads";
    });

TEST_P(HybridSeqTest, RankAndSize) {
    const int T = GetParam();
    RunHybrid(T, [&](int t, int nthreads, HybridComm& comm) {
        EXPECT_EQ(comm.Rank(), t);
        EXPECT_EQ(comm.Size(), nthreads);
    });
}

TEST_P(HybridSeqTest, AllreduceSumInt) {
    const int T = GetParam();
    std::vector<int> results(T, 0);
    RunHybrid(T, [&](int t, int, HybridComm& comm) {
        int send = t + 1, recv = 0;
        comm.Allreduce(&send, &recv, 1, CommDatatype::INT, CommOp::SUM);
        results[t] = recv;
    });
    const int expected = T * (T + 1) / 2;
    for (int r : results)
        EXPECT_EQ(r, expected);
}

TEST_P(HybridSeqTest, AllreduceSumInPlace) {
    const int T = GetParam();
    std::vector<int> vals(T);
    std::iota(vals.begin(), vals.end(), 1); // 1, 2, ..., T
    RunHybrid(T, [&](int t, int, HybridComm& comm) {
        comm.Allreduce(COMM_IN_PLACE, &vals[t], 1, CommDatatype::INT, CommOp::SUM);
    });
    const int expected = T * (T + 1) / 2;
    for (int v : vals)
        EXPECT_EQ(v, expected);
}

TEST_P(HybridSeqTest, AllreduceMaxInt) {
    const int T = GetParam();
    std::vector<int> results(T, 0);
    RunHybrid(T, [&](int t, int, HybridComm& comm) {
        int send = t * 10, recv = 0;
        comm.Allreduce(&send, &recv, 1, CommDatatype::INT, CommOp::MAX);
        results[t] = recv;
    });
    const int expected = (T - 1) * 10;
    for (int r : results)
        EXPECT_EQ(r, expected);
}

TEST_P(HybridSeqTest, AllreduceMinInt) {
    const int T = GetParam();
    std::vector<int> results(T, -1);
    RunHybrid(T, [&](int t, int, HybridComm& comm) {
        int send = t * 10 + 100, recv = -1;
        comm.Allreduce(&send, &recv, 1, CommDatatype::INT, CommOp::MIN);
        results[t] = recv;
    });
    for (int r : results)
        EXPECT_EQ(r, 100);
}

TEST_P(HybridSeqTest, AllreduceSumArray) {
    const int T = GetParam();
    const int N = 3;
    std::vector<std::vector<int>> results(T, std::vector<int>(N, 0));
    RunHybrid(T, [&](int t, int, HybridComm& comm) {
        std::vector<int> send(N);
        for (int i = 0; i < N; ++i)
            send[i] = (t + 1) * (i + 1);
        comm.Allreduce(send.data(), results[t].data(), N, CommDatatype::INT, CommOp::SUM);
    });
    // expected[i] = (i+1) * sum(1..T)
    for (int t = 0; t < T; ++t)
        for (int i = 0; i < N; ++i)
            EXPECT_EQ(results[t][i], (i + 1) * T * (T + 1) / 2) << "t=" << t << " i=" << i;
}

TEST_P(HybridSeqTest, ReduceToRootThread) {
    const int T = GetParam();
    // Use a separate recvbuf per thread; only thread 0's gets written.
    std::vector<int> recvbufs(T, -1);
    RunHybrid(T, [&](int t, int, HybridComm& comm) {
        int send = t + 1;
        comm.Reduce(&send, &recvbufs[t], 1, CommDatatype::INT, CommOp::SUM, /*root=*/0);
    });
    EXPECT_EQ(recvbufs[0], T * (T + 1) / 2);
}

TEST_P(HybridSeqTest, BcastFromRootThread) {
    const int T = GetParam();
    std::vector<int> bufs(T, -1);
    bufs[0] = 42;
    RunHybrid(T, [&](int t, int, HybridComm& comm) {
        comm.Bcast(&bufs[t], 1, CommDatatype::INT, /*root=*/0);
    });
    for (int v : bufs)
        EXPECT_EQ(v, 42);
}

TEST_P(HybridSeqTest, BcastFromNonZeroRoot) {
    const int T = GetParam();
    if (T < 2) GTEST_SKIP() << "need T>=2 for a non-zero root";
    std::vector<int> bufs(T, -1);
    bufs[1] = 77; // root is virtual rank 1
    RunHybrid(T, [&](int t, int, HybridComm& comm) {
        comm.Bcast(&bufs[t], 1, CommDatatype::INT, /*root=*/1);
    });
    for (int v : bufs)
        EXPECT_EQ(v, 77);
}

TEST_P(HybridSeqTest, ExscanSumInt) {
    const int T = GetParam();
    std::vector<int> results(T, -1);
    RunHybrid(T, [&](int t, int, HybridComm& comm) {
        int send = t + 1, recv = -1;
        comm.Exscan(&send, &recv, 1, CommDatatype::INT, CommOp::SUM);
        results[t] = recv;
    });
    // exclusive scan of [1,2,...,T] → [0,1,3,6,...]
    int prefix = 0;
    for (int t = 0; t < T; ++t) {
        EXPECT_EQ(results[t], prefix) << "t=" << t;
        prefix += t + 1;
    }
}

TEST_P(HybridSeqTest, AllgatherUniform) {
    const int T = GetParam();
    std::vector<std::vector<int>> results(T, std::vector<int>(T, -1));
    RunHybrid(T, [&](int t, int, HybridComm& comm) {
        int send = t;
        comm.Allgather(&send, 1, CommDatatype::INT, results[t].data(), 1, CommDatatype::INT);
    });
    std::vector<int> expected(T);
    std::iota(expected.begin(), expected.end(), 0);
    for (int t = 0; t < T; ++t)
        EXPECT_THAT(results[t], ElementsAreArray(expected)) << "t=" << t;
}

TEST_P(HybridSeqTest, AllgathervVariableCounts) {
    const int T = GetParam();
    // Thread t sends (t+1) ints all equal to t.
    // recvcounts[v]=v+1; displs[v]=v*(v+1)/2
    const int total = T * (T + 1) / 2;
    std::vector<int> recvcounts(T), displs(T);
    for (int v = 0; v < T; ++v) {
        recvcounts[v] = v + 1;
        displs[v]     = v * (v + 1) / 2;
    }
    std::vector<std::vector<int>> results(T, std::vector<int>(total, -1));
    RunHybrid(T, [&](int t, int, HybridComm& comm) {
        std::vector<int> send(t + 1, t);
        comm.Allgatherv(
            send.data(), t + 1, CommDatatype::INT,
            results[t].data(), recvcounts.data(), displs.data(), CommDatatype::INT);
    });
    // expected: for virtual PE v, (v+1) copies of v
    std::vector<int> expected;
    for (int v = 0; v < T; ++v)
        for (int i = 0; i <= v; ++i)
            expected.push_back(v);
    for (int t = 0; t < T; ++t)
        EXPECT_THAT(results[t], ElementsAreArray(expected)) << "t=" << t;
}

TEST_P(HybridSeqTest, AlltoallUniform) {
    const int T = GetParam();
    // Thread t sends sendbuf[dest] = t * T + dest
    std::vector<std::vector<int>> sendbufs(T, std::vector<int>(T));
    std::vector<std::vector<int>> recvbufs(T, std::vector<int>(T, -1));
    for (int t = 0; t < T; ++t)
        for (int d = 0; d < T; ++d)
            sendbufs[t][d] = t * T + d;
    RunHybrid(T, [&](int t, int, HybridComm& comm) {
        comm.Alltoall(
            sendbufs[t].data(), 1, CommDatatype::INT,
            recvbufs[t].data(), 1, CommDatatype::INT);
    });
    // Thread t recvbuf[src] = src * T + t
    for (int t = 0; t < T; ++t)
        for (int src = 0; src < T; ++src)
            EXPECT_EQ(recvbufs[t][src], src * T + t) << "t=" << t << " src=" << src;
}

TEST_P(HybridSeqTest, AlltoallvNonUniform) {
    const int T = GetParam();
    // Thread t sends value (t * 100 + dest) to thread dest (1 element each)
    std::vector<std::vector<int>> sendbufs(T, std::vector<int>(T));
    std::vector<std::vector<int>> recvbufs(T, std::vector<int>(T, -1));
    std::vector<int> sendcounts(T, 1), sdispls(T), recvcounts(T, 1), rdispls(T);
    for (int i = 0; i < T; ++i) {
        sdispls[i] = i;
        rdispls[i] = i;
    }
    for (int t = 0; t < T; ++t)
        for (int dest = 0; dest < T; ++dest)
            sendbufs[t][dest] = t * 100 + dest;
    RunHybrid(T, [&](int t, int, HybridComm& comm) {
        comm.Alltoallv(
            sendbufs[t].data(), sendcounts.data(), sdispls.data(), CommDatatype::INT,
            recvbufs[t].data(), recvcounts.data(), rdispls.data(), CommDatatype::INT);
    });
    // Thread t recvbuf[src] = src * 100 + t
    for (int t = 0; t < T; ++t)
        for (int src = 0; src < T; ++src)
            EXPECT_EQ(recvbufs[t][src], src * 100 + t) << "t=" << t << " src=" << src;
}

TEST_P(HybridSeqTest, GatherToRootThread) {
    const int T = GetParam();
    std::vector<int> rootbuf(T, -1);
    RunHybrid(T, [&](int t, int, HybridComm& comm) {
        int send = t * 10;
        // All threads pass rootbuf; only root thread 0's copy gets written.
        comm.Gather(&send, 1, CommDatatype::INT, rootbuf.data(), 1, CommDatatype::INT, /*root=*/0);
    });
    for (int t = 0; t < T; ++t)
        EXPECT_EQ(rootbuf[t], t * 10) << "position " << t;
}

TEST_P(HybridSeqTest, Barrier) {
    const int T = GetParam();
    std::atomic<int> counter{0};
    RunHybrid(T, [&](int, int, HybridComm& comm) {
        ++counter;
        comm.Barrier();
        EXPECT_EQ(counter.load(), T);
    });
}

// ===========================================================================
// HybridComm over MPIComm  (P MPI ranks × T threads per rank)
// Tests inter-rank communication with different P×T combinations.
//   CORES 1:  1 MPI rank × {1,2,4} threads  (1×1, 1×2, 1×4)
//   CORES 2:  2 MPI ranks × {1,2,4} threads  (2×1, 2×2, 2×4)
//   CORES 4:  4 MPI ranks × {1,2,4} threads  (4×1, 4×2, 4×4)
// Notably, CORES 4 T=1 gives 4×1, CORES 2 T=2 gives 2×2, CORES 1 T=4 gives 1×4.
// ===========================================================================

#ifndef KAGEN_NOMPI

class HybridMPITest : public ::testing::TestWithParam<int> {};

INSTANTIATE_TEST_SUITE_P(
    ThreadCounts, HybridMPITest,
    ::testing::Values(1, 2, 4),
    [](const ::testing::TestParamInfo<int>& info) {
        return std::to_string(info.param) + "threads";
    });

TEST_P(HybridMPITest, RankAndSize) {
    const int T = GetParam();
    MPIComm mpi(MPI_COMM_WORLD);
    const int P = mpi.Size();
    const int R = mpi.Rank();
    RunHybrid(T, mpi, [&](int t, int, HybridComm& comm) {
        EXPECT_EQ(comm.Rank(), R * T + t);
        EXPECT_EQ(comm.Size(), P * T);
    });
}

TEST_P(HybridMPITest, AllreduceSumAcrossAllVirtualPEs) {
    const int T = GetParam();
    MPIComm mpi(MPI_COMM_WORLD);
    const int V = mpi.Size() * T;
    std::vector<long long> results(T, 0);
    RunHybrid(T, mpi, [&](int t, int, HybridComm& comm) {
        long long send = comm.Rank() + 1; // virtual rank +1
        comm.Allreduce(&send, &results[t], 1, CommDatatype::LONG_LONG, CommOp::SUM);
    });
    const long long expected = (long long)V * (V + 1) / 2;
    for (int t = 0; t < T; ++t)
        EXPECT_EQ(results[t], expected) << "t=" << t;
}

TEST_P(HybridMPITest, ExscanSumAcrossAllVirtualPEs) {
    // Each virtual PE contributes 1; exscan result for PE r equals r.
    const int T = GetParam();
    MPIComm mpi(MPI_COMM_WORLD);
    const int R = mpi.Rank();
    std::vector<long long> results(T, -1);
    RunHybrid(T, mpi, [&](int t, int, HybridComm& comm) {
        long long send = 1, recv = -1;
        comm.Exscan(&send, &recv, 1, CommDatatype::LONG_LONG, CommOp::SUM);
        results[t] = recv;
    });
    for (int t = 0; t < T; ++t)
        EXPECT_EQ(results[t], (long long)(R * T + t)) << "R=" << R << " t=" << t;
}

TEST_P(HybridMPITest, AllgatherAcrossAllVirtualPEs) {
    const int T = GetParam();
    MPIComm mpi(MPI_COMM_WORLD);
    const int V = mpi.Size() * T;
    std::vector<std::vector<int>> recvbufs(T, std::vector<int>(V, -1));
    RunHybrid(T, mpi, [&](int t, int, HybridComm& comm) {
        int send = (int)comm.Rank();
        comm.Allgather(&send, 1, CommDatatype::INT, recvbufs[t].data(), 1, CommDatatype::INT);
    });
    std::vector<int> expected(V);
    std::iota(expected.begin(), expected.end(), 0);
    for (int t = 0; t < T; ++t)
        EXPECT_THAT(recvbufs[t], ElementsAreArray(expected)) << "t=" << t;
}

TEST_P(HybridMPITest, AlltoallUniformAcrossAllVirtualPEs) {
    // Virtual PE r sends to virtual PE d: value = r * V + d
    // After Alltoall (sendcount=1): virtual PE r recvbuf[src] = src * V + r
    const int T = GetParam();
    MPIComm mpi(MPI_COMM_WORLD);
    const int V = mpi.Size() * T;
    const int R = mpi.Rank();

    std::vector<std::vector<int>> sendbufs(T, std::vector<int>(V));
    std::vector<std::vector<int>> recvbufs(T, std::vector<int>(V, -1));
    for (int t = 0; t < T; ++t) {
        const int r = R * T + t;
        for (int d = 0; d < V; ++d)
            sendbufs[t][d] = r * V + d;
    }
    RunHybrid(T, mpi, [&](int t, int, HybridComm& comm) {
        comm.Alltoall(
            sendbufs[t].data(), 1, CommDatatype::INT,
            recvbufs[t].data(), 1, CommDatatype::INT);
    });
    for (int t = 0; t < T; ++t) {
        const int r = R * T + t;
        for (int src = 0; src < V; ++src)
            EXPECT_EQ(recvbufs[t][src], src * V + r) << "R=" << R << " t=" << t << " src=" << src;
    }
}

// ===========================================================================
// MPIComm basics
// ===========================================================================

TEST(MPIComm, RankAndSize) {
    MPIComm comm(MPI_COMM_WORLD);
    EXPECT_GE(comm.Rank(), 0);
    EXPECT_GE(comm.Size(), 1);
    EXPECT_LT(comm.Rank(), comm.Size());
}

TEST(MPIComm, AllreduceSumCountsRanks) {
    MPIComm comm(MPI_COMM_WORLD);
    int send = 1, recv = 0;
    comm.Allreduce(&send, &recv, 1, CommDatatype::INT, CommOp::SUM);
    EXPECT_EQ(recv, comm.Size());
}

TEST(MPIComm, ExscanSumGivesRank) {
    MPIComm comm(MPI_COMM_WORLD);
    int send = 1, recv = -1;
    comm.Exscan(&send, &recv, 1, CommDatatype::INT, CommOp::SUM);
    // MPI leaves rank 0's receive buffer undefined; only assert for rank > 0.
    if (comm.Rank() > 0) {
        EXPECT_EQ(recv, comm.Rank());
    }
}

TEST(MPIComm, AllgatherCollectsRanks) {
    MPIComm comm(MPI_COMM_WORLD);
    const int P = comm.Size();
    int send = comm.Rank();
    std::vector<int> recv(P, -1);
    comm.Allgather(&send, 1, CommDatatype::INT, recv.data(), 1, CommDatatype::INT);
    std::vector<int> expected(P);
    std::iota(expected.begin(), expected.end(), 0);
    EXPECT_THAT(recv, ElementsAreArray(expected));
}

#endif // !KAGEN_NOMPI
