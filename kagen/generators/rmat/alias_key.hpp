/*******************************************************************************
 * rmat/alias_key.hpp
 *
 * Alias method implementations for key-weight-pair inputs
 *
 * Copyright (C) 2018 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once

#include "kagen/generators/rmat/memory.hpp"
#include "kagen/generators/rmat/timer.hpp"
#include "kagen/tlx/attribute_always_inline.hpp"

#include <cassert>
#include <limits>
#include <memory>
#include <numeric>
#include <ostream>
#include <utility>
#include <vector>

#ifndef NDEBUG
    #include <unordered_map>
#endif

namespace rmat {

template <typename key_type = uint32_t>
struct alias_key {
    using result_type = key_type;
    using pair        = std::pair<key_type, double>;

    static_assert((sizeof(double) / sizeof(key_type)) * sizeof(key_type) == sizeof(double), "key_type has weird size");
    static constexpr size_t ti_ptroffset = sizeof(double) / sizeof(key_type);

    // Alias table member item, used to be a std::tuple<double, key_type, key_type>
    struct tableitem {
        double   p; // probability
        key_type i; // item
        key_type a; // alia

        tableitem() : p(0), i(-1), a(-1) {}
        tableitem(double p_, key_type i_, key_type a_) : p(p_), i(i_), a(a_) {}
    };

    friend std::ostream& operator<<(std::ostream& os, const tableitem& item) {
        return os << '(' << item.p << ',' << item.i << ',' << item.a << ')';
    }

    static constexpr bool debug = false;
    static constexpr bool time  = false;

    alias_key() : table_(nullptr), work_(nullptr), size_(-1) {}

    template <typename Iterator>
    alias_key(Iterator w_begin, Iterator w_end) {
        init(w_end - w_begin);
        construct(w_begin, w_end);
    }

    void init(size_t size) {
        table_ = make_alloc_arr<tableitem>(size);
        work_  = make_alloc_arr<pair>(size);
        size_  = size;
    }

    template <typename Iterator>
    void construct(Iterator begin, Iterator end, bool is_dist = false) {
        timer t;
        if (end - begin != static_cast<ssize_t>(size_)) {
            /*
            sLOG1 << "Error: tried to construct alias table of incorrect size!"
                  << "Expected" << size_ << "got" << end - begin;
            */
            return;
        }
        timers_.clear();

#ifdef NDEBUG
        if (is_dist) {
            W_ = 1.0;
        } else {
#endif
            W_ = std::accumulate(begin, end, 0.0, [](double curr, auto b) { return curr + b.second; });
            // If the input is supposed to be a probability distribution and
            // we're in debug mode, check that the weights really sum to 1
            if (is_dist) {
                assert(std::abs(W_ - 1.0) < 1e-6);
                W_ = 1.0;
            }
#ifdef NDEBUG
        }
#endif

        W_n_ = W_ / size_;
        // LOG << "W = " << W_ << ", W/n = " << W_n_;
        timers_.push_back(t.get_and_reset());
        // LOGC(time) << "Step 0: preprocessing took " << timers_.back() << "ms";

        // Classify into small and large items
        ssize_t i_small = 0, i_large = size_ - 1;
        for (Iterator it = begin; it != end; ++it) {
            if (is_small(it->second)) {
                // LOG << *it << " is small";
                work_[i_small++] = *it;
            } else {
                // LOG << *it << " is large";
                work_[i_large--] = *it;
            }
        }
        assert(i_small > i_large); // must meet
        timers_.push_back(t.get_and_reset());
        // LOGC(time) << "Step 1: classification took " << timers_.back() << "ms";

        // Assign items
        auto s_begin = work_.get(), s_end = work_.get() + i_small, l_begin = s_end, l_end = work_.get() + size_,
             small = s_begin, large = l_begin;
        size_t table_idx = 0;
        // LOG_ARR(work_.get(), "work");
        while (small != s_end && large != l_end) {
            /*
            LOG << "table_[" << table_idx << "] = (" << small->second << ", " << small->first << ", " << large->first
                << ")";
            */
            table_[table_idx++] = tableitem(small->second, small->first, large->first);
            large->second       = (large->second + small->second) - W_n_;

            if (is_small(large)) {
                // sLOG << "large item became small, remaining piece" << *large << "moving to end of small items";
                if (large > s_end) {
                    std::swap(large, s_end);
                } else {
                    // LOG << "no need to swap, is first";
                }
                s_end++;
                // l_begin++; // not really needed, just for bookkeping
                large++;
            } else {
                // sLOG << "large item" << *large << "remained large, advancing small";
            }
            small++;
        }
        while (large != l_end) {
            size_t i_large = large->first;
            // LOG << "large item left at the end: " << *large;
            // LOG << "table_[" << table_idx << "] = (" << W_n_ << ", " << i_large << ", -1)";

            table_[table_idx++] = tableitem(W_n_, i_large, i_large);
            large++;
        }
        while (small != s_end) {
            size_t i_small = small->first;
            /*
            sLOG << "encountered some numerical instability!"
                 << "Item" << i_small << "weight" << *small;
            LOG << "table_[" << table_idx << "] = (" << W_n_ << ", " << i_small << ", -1)";
            */

            table_[table_idx++] = tableitem(W_n_, i_small, i_small);
            small++;
        }
        timers_.push_back(t.get_and_reset());
        // LOGC(time) << "Step 2: assignment took " << timers_.back() << "ms";
    }

    // Query given a uniformly distributed [0,1) random value
    key_type sample(double uniform) const {
        double     rand  = uniform * size_;
        size_t     index = rand;
        tableitem& item  = table_[index];

        // sanity check: if alias == -1, then prob == 1
        assert(item.a != static_cast<key_type>(-1) || std::abs(item.p - 1) < 1e-10);

        rand          = (rand - index) * W_n_;
        size_t offset = rand >= item.p;
        assert(item.a != static_cast<key_type>(-1) || offset == 0);
        // ugly hack
        return *(reinterpret_cast<key_type*>(&item) + ti_ptroffset + offset);
    }

    template <typename Iterator>
    void verify(Iterator begin, Iterator end) const {
        (void)begin;
        (void)end;
#ifndef NDEBUG
        assert(end - begin == static_cast<ssize_t>(size_));
        std::unordered_map<key_type, double> weights(size_);

        for (size_t i = 0; i < size_; i++) {
            auto [p, i1, i2] = table_[i];
            assert(p > 0);
            assert(i1 != static_cast<key_type>(-1));
            weights[i1] += p;
            if (p < W_n_) {
                assert(i2 != static_cast<key_type>(-1));
                weights[i2] += W_n_ - p;
            }
        }

        for (Iterator it = begin; it != end; ++it) {
            key_type key    = it->first;
            double   should = it->second;
            double   have   = weights[key];
            // sLOG << "key" << key << "should be" << should << "have" << have;
            assert(std::abs(should - have) < W_ * 1e-10);
        }
        // LOG1 << "Verification succeeded!";
#endif
    }

    std::vector<double> get_timers() const {
        return timers_;
    }

    double total_weight() const {
        return W_;
    }

    size_t size() const {
        return size_;
    }

protected:
    TLX_ATTRIBUTE_ALWAYS_INLINE
    bool is_small(double& d) const {
        return d <= W_n_;
    }

    template <typename Iterator>
    TLX_ATTRIBUTE_ALWAYS_INLINE bool is_small(Iterator it) const {
        return it->second <= W_n_;
    }

    alloc_arr_ptr<tableitem> table_;
    alloc_arr_ptr<pair>      work_;
    std::vector<double>      timers_;
    size_t                   size_;
    double                   W_, W_n_;
};

} // namespace rmat

