#pragma once

#include "kagen/comm/seq_comm.h"

namespace kagen {

// A SeqComm that reports a custom (virtual) rank and size.
// Used for intra-process threaded generation: each thread gets a VirtualComm
// whose Rank() = real_rank * num_threads + thread_id and
// whose Size() = real_size * num_threads.
// All collective operations degenerate to local copies / no-ops (inherited from SeqComm)
// because generation is communication-free — each thread generates its chunk independently.
class VirtualComm : public SeqComm {
public:
    VirtualComm(CommPEID rank, CommPEID size) : rank_(rank), size_(size) {}

    CommPEID Rank() const override {
        return rank_;
    }

    CommPEID Size() const override {
        return size_;
    }

private:
    CommPEID rank_;
    CommPEID size_;
};

} // namespace kagen
