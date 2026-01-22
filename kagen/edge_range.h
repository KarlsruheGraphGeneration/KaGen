#pragma once

#include "kagen/kagen.h"

#include <cstddef>
#include <iterator>
#include <utility>

namespace kagen {

class EdgeRange {
public:
    using Vertex = SInt;
    using Edge   = std::pair<Vertex, Vertex>;

    explicit EdgeRange(const Edgelist& edgelist) noexcept;
    EdgeRange(const XadjArray& xadj, const AdjncyArray& adjncy, VertexRange vertex_range) noexcept;
    EdgeRange(const Graph& graph)noexcept;

    static EdgeRange FromGraph(const Graph& graph) noexcept;

    class iterator {
    public:
        using iterator_category = std::bidirectional_iterator_tag;
        using value_type        = Edge;
        using difference_type   = std::ptrdiff_t;

        iterator() noexcept = default;

        static iterator edgelist_begin(const EdgeRange* parent) noexcept;
        static iterator edgelist_end(const EdgeRange* parent) noexcept;

        // begin / end for CSR
        static iterator csr_begin(const EdgeRange* parent) noexcept;
        static iterator csr_end(const EdgeRange* parent) noexcept;

        std::size_t edge_index() const noexcept;

        value_type operator*() const noexcept;

        iterator& operator++();

        iterator operator++(int);

        iterator& operator--();

        iterator operator--(int);

        bool operator==(const iterator& other) const noexcept;
        bool operator!=(const iterator& other) const noexcept;

    private:
        const EdgeRange* parent_{nullptr};

        // EDGELIST state
        std::size_t idx_{0};

        // CSR state
        std::size_t u_{0};   // current vertex index (0..n_local-1)
        std::size_t off_{0}; // current offset into adjncy

        Edge cur_{};

        void load_current() noexcept;
        void init_csr_begin() noexcept;
        void advance_to_next_valid_csr() noexcept;
        void retreat_to_prev_valid_csr() noexcept;
    };

    iterator begin() const noexcept;
    iterator end() const noexcept;

    std::size_t size() const noexcept;

private:
    GraphRepresentation representation_;

    const Edgelist*    edgelist_ptr_;
    const XadjArray*   xadj_ptr_;
    const AdjncyArray* adjncy_ptr_;
    Vertex             vertex_base_;

    friend class iterator;
};

} // namespace kagen
