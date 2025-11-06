#include "kagen/edge_range.h"

#include <cassert>
#include <memory>

namespace kagen {

EdgeRange::EdgeRange(const Edgelist& edgelist) noexcept
    : representation_(GraphRepresentation::EDGE_LIST),
      edgelist_ptr_(std::addressof(edgelist)),
      xadj_ptr_(nullptr),
      adjncy_ptr_(nullptr),
      vertex_base_(Vertex(0)) {}

EdgeRange::EdgeRange(const XadjArray& xadj, const AdjncyArray& adjncy, VertexRange vertex_range) noexcept
    : representation_(GraphRepresentation::CSR),
      edgelist_ptr_(nullptr),
      xadj_ptr_(std::addressof(xadj)),
      adjncy_ptr_(std::addressof(adjncy)),
      vertex_base_(vertex_range.first) {
    assert(xadj_ptr_->size() >= 1);
}

EdgeRange::EdgeRange(const Graph& graph) noexcept : EdgeRange(FromGraph(graph)) {}

EdgeRange EdgeRange::FromGraph(const Graph& g) noexcept {
    if (g.representation == GraphRepresentation::EDGE_LIST) {
        return EdgeRange(g.edges);
    } else {
        return EdgeRange(g.xadj, g.adjncy, g.vertex_range);
    }
}

EdgeRange::iterator EdgeRange::iterator::edgelist_begin(const EdgeRange* parent) noexcept {
    iterator it;
    it.parent_ = parent;
    it.idx_    = 0;
    it.load_current();
    return it;
}

EdgeRange::iterator EdgeRange::iterator::edgelist_end(const EdgeRange* parent) noexcept {
    iterator it;
    it.parent_ = parent;
    it.idx_    = parent->edgelist_ptr_->size();
    return it;
}

EdgeRange::iterator EdgeRange::iterator::csr_begin(const EdgeRange* parent) noexcept {
    iterator it;
    it.parent_ = parent;
    it.u_      = 0;
    it.off_    = 0;
    it.init_csr_begin(); // finds first valid (u,off) or sets to end
    it.load_current();
    return it;
}

EdgeRange::iterator EdgeRange::iterator::csr_end(const EdgeRange* parent) noexcept {
    iterator it;
    it.parent_ = parent;
    it.u_      = parent->xadj_ptr_->size() == 0 ? 0u : parent->xadj_ptr_->size() - 1;
    it.off_    = parent->adjncy_ptr_->size();
    return it;
}

std::size_t EdgeRange::iterator::edge_index() const noexcept {
    if (parent_->representation_ == GraphRepresentation::EDGE_LIST) {
        return idx_;
    } else {
        return off_;
    }
}

EdgeRange::iterator::value_type EdgeRange::iterator::operator*() const noexcept {
    return cur_;
}

EdgeRange::iterator& EdgeRange::iterator::operator++() {
    assert(parent_ != nullptr);

    if (parent_->representation_ == GraphRepresentation::EDGE_LIST) {
        ++idx_;
        load_current();
        return *this;
    }

    advance_to_next_valid_csr();
    load_current();
    return *this;
}

EdgeRange::iterator EdgeRange::iterator::operator++(int) {
    iterator tmp = *this;
    ++(*this);
    return tmp;
}

EdgeRange::iterator& EdgeRange::iterator::operator--() {
    assert(parent_ != nullptr);

    if (parent_->representation_ == GraphRepresentation::EDGE_LIST) {
        assert(idx_ > 0 && "Cannot decrement begin iterator");
        --idx_;
        load_current();
        return *this;
    }

    retreat_to_prev_valid_csr();
    load_current();
    return *this;
}

EdgeRange::iterator EdgeRange::iterator::operator--(int) {
    iterator tmp = *this;
    --(*this);
    return tmp;
}

bool EdgeRange::iterator::operator==(const iterator& other) const noexcept {
    // must belong to same parent to compare reliably
    if (parent_ != other.parent_)
        return false;
    if (parent_ == nullptr)
        return true; // both default constructed?
    if (parent_->representation_ != other.parent_->representation_)
        return false;
    if (parent_->representation_ == GraphRepresentation::EDGE_LIST) {
        return idx_ == other.idx_;
    } else {
        return (u_ == other.u_) && (off_ == other.off_);
    }
}

bool EdgeRange::iterator::operator!=(const iterator& other) const noexcept {
    return !(*this == other);
}

void EdgeRange::iterator::load_current() noexcept {
    if (parent_->representation_ == GraphRepresentation::EDGE_LIST) {
        const auto& elist = *parent_->edgelist_ptr_;
        if (idx_ < elist.size()) {
            cur_ = elist[idx_];
        }
        // else leave cur_ unspecified for end()
    } else {
        const auto&       xadj    = *parent_->xadj_ptr_;
        const auto&       adjncy  = *parent_->adjncy_ptr_;
        const std::size_t n_local = xadj.empty() ? 0 : (xadj.size() - 1);
        if (u_ >= n_local || off_ >= adjncy.size()) {
            // end iterator; cur_ remains unspecified
            return;
        }
        Vertex u_global = static_cast<Vertex>(u_) + parent_->vertex_base_;
        Vertex v_global = adjncy[off_];
        cur_.first      = u_global;
        cur_.second     = v_global;
    }
}

void EdgeRange::iterator::init_csr_begin() noexcept {
    const auto&       xadj    = *parent_->xadj_ptr_;
    const auto&       adjncy  = *parent_->adjncy_ptr_;
    const std::size_t n_local = xadj.empty() ? 0 : (xadj.size() - 1);
    if (n_local == 0) {
        // empty CSR: set end
        u_   = 0;
        off_ = 0;
        return;
    }
    off_ = xadj[0];
    u_   = 0;
    // skip vertices with empty adjacency ranges
    while (u_ < n_local && off_ >= xadj[u_ + 1]) {
        ++u_;
        if (u_ < n_local)
            off_ = xadj[u_];
    }
    if (u_ >= n_local) {
        // set to end state
        u_   = n_local;
        off_ = adjncy.size();
    }
}

void EdgeRange::iterator::advance_to_next_valid_csr() noexcept {
    // CSR: advance off_, move to next vertex if necessary
    ++off_;
    const auto&       xadj    = *parent_->xadj_ptr_;
    const std::size_t n_local = xadj.empty() ? 0 : (xadj.size() - 1);
    while (u_ < n_local && off_ >= xadj[u_ + 1]) {
        ++u_;
        if (u_ < n_local)
            off_ = xadj[u_];
    }
    if (u_ >= n_local) {
        // set canonical end state
        off_ = parent_->adjncy_ptr_->size();
    }
}

void EdgeRange::iterator::retreat_to_prev_valid_csr() noexcept {
    const auto&       xadj    = *parent_->xadj_ptr_;
    const auto&       adjncy  = *parent_->adjncy_ptr_;
    const std::size_t n_local = xadj.empty() ? 0 : (xadj.size() - 1);

    // If we're at end(), move to the last valid edge
    if (off_ >= adjncy.size()) {
        assert(adjncy.size() > 0 && "Cannot decrement begin iterator");
        off_ = adjncy.size() - 1;
        // Find which vertex this edge belongs to
        for (std::size_t v = 0; v < n_local; ++v) {
            if (xadj[v] <= off_ && off_ < xadj[v + 1]) {
                u_ = v;
                return;
            }
        }
        // Should not reach here if data is consistent
        assert(false && "Invalid state in retreat_to_prev_valid_csr");
        return;
    }

    // Check if we can move back within the current vertex's adjacency list
    if (off_ > 0 && xadj[u_] < off_) {
        --off_;
        return;
    }

    // Move to the previous vertex
    assert(u_ > 0 && "Cannot decrement begin iterator");
    --u_;

    // Find the last edge of the previous vertex (skip empty vertices going backward)
    while (u_ < n_local && xadj[u_] >= xadj[u_ + 1]) {
        if (u_ == 0) {
            // No valid edges before this point
            assert(false && "Cannot decrement begin iterator");
            return;
        }
        --u_;
    }

    // Set off_ to the last edge of vertex u_
    off_ = xadj[u_ + 1] - 1;
}

EdgeRange::iterator EdgeRange::begin() const noexcept {
    if (representation_ == GraphRepresentation::EDGE_LIST)
        return iterator::edgelist_begin(this);
    return iterator::csr_begin(this);
}

EdgeRange::iterator EdgeRange::end() const noexcept {
    if (representation_ == GraphRepresentation::EDGE_LIST)
        return iterator::edgelist_end(this);
    return iterator::csr_end(this);
}

std::size_t EdgeRange::size() const noexcept {
    if (representation_ == GraphRepresentation::EDGE_LIST) {
        return edgelist_ptr_->size();
    } else {
        return adjncy_ptr_->size();
    }
}

} // namespace kagen
