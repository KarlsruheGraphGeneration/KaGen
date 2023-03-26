/*******************************************************************************
 * sampling/methodB.hpp
 *
 * Copyright (C) 2016-2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef SAMPLING_METHOD_B_HEADER
    #define SAMPLING_METHOD_B_HEADER

    #include "kagen/sampling/benchmark.hpp"
    #include "kagen/sampling/config.hpp"
    #include "kagen/sampling/methodH.hpp"
    #include "kagen/sampling/methodR.hpp"
    #include "kagen/sampling/methodSH.hpp"
    #include "kagen/sampling/rng/select.hpp"
    #include "kagen/sampling/timer.hpp"

    #include <algorithm>
    #include <cassert>
    #include <cmath>
    #include <cstring>
    #include <memory>
    #include <random>
    #include <immintrin.h>

namespace sampling {

struct global_stats {
    void push_sum(double t) {
        s_sum.push(t);
    }
    void push_gen(double t) {
        s_gen.push(t);
    }
    void push_prefsum(double t) {
        s_prefsum.push(t);
    }
    void push_fix(double t) {
        s_fix.push(t);
    }
    statistics s_sum, s_gen, s_prefsum, s_fix;
};

template <typename Generator = rng::select_t>
class method_B {
public:
    using generator_t = Generator;

    method_B(size_t seed) : rng(seed), seed(seed) {}

    /**
     * @param dest output iterator
     * @param size size of array (last valid pos is dest + size - 1)
     * @param k number of samples to draw
     * @param universe value range of samples [0..universe)
     */
    template <typename int_t>
    void sample(
        std::vector<int_t>& output, size_t k, size_t universe, global_stats* stats = nullptr,
        const bool verbose = false) {
        ((void)verbose);

        double p;
        size_t ssize;
        std::tie(p, ssize) = calc_params(universe, k);
        output.resize(ssize);

        auto   pos = output.begin(), dest = output.begin();
        double t_gen(0.0), t_pref(0.0), t_fix(0.0);
        // size_t its = 0; // how many iterations it took
        timer  t, t_total;
        size_t usable_samples = 0;

        do {
            t.reset();
            rng.generate_geometric_block(p, output, ssize);
            t_gen += t.get_and_reset();

            inplace_prefix_sum_disp<true>(dest, dest + ssize);
            t_pref += t.get_and_reset();

            // find out how many samples are within the universe.
            // pos is the first one that's too large.
            pos            = std::lower_bound(dest, dest + ssize, universe);
            usable_samples = pos - dest;

            /*
            ++its;
            LOGC(verbose) << "\tIt " << its << ": got " << usable_samples
                          << " samples in range (" << k << " of " << ssize
                          << " required) => "
                          << (long long)usable_samples - (long long)k
                          << " to delete, "  << ssize - usable_samples
                          << " outside universe ignored (largest: "
                          << *(dest + ssize - 1) << ")";
            */
        } while (pos == (dest + ssize) || usable_samples < k);
        // first condition is that the whole universe is covered (i.e. ssize was
        // big enough) and second condition is that >= k samples lie within it

        t.reset();
        if (usable_samples > k) {
            // pick k out of the pos-dest-1 elements
    #ifdef FIX_STABLE
            pos = fix_stable(dest, pos, k);
    #else
            pos = fix(dest, pos, k);
    #endif
            // pos is now the past-the-end iterator of the sample indices
            assert(pos - dest == (long)k);
        }
        t_fix += t.get_and_reset();

        double t_sum = t_total.get();
        if (stats != nullptr) {
            stats->push_sum(t_sum);
            stats->push_gen(t_gen);
            stats->push_prefsum(t_pref);
            stats->push_fix(t_fix);
        }

        /*
        // Disabled: return timings as log string
        std::stringstream stream;
        stream << "INFO"
               << " time=" << t_sum
               << " k=" << k
               << " b=" << ssize
               << " restarts=" << its-1
               << " t_gen=" << t_gen
               << " t_prefsum=" << t_pref
               << " t_fix=" << t_fix
#ifdef FIX_STABLE
               << " fixer=stable";
#else
               << " fixer=fast";
#endif
        return stream.str();
        */
    }

private:
    // Formulas from "Sequential Random Sampling" by Ahrens and Dieter, 1985
    static std::pair<double, size_t> calc_params(size_t universe, size_t k /* samples */) {
        double r = sqrt(k);
        double a = sqrt(log(1 + k / (2 * M_PI)));
        a        = a + a * a / (3.0 * r);
        size_t b = k + size_t(4 * a * r);
        double p = (k + a * r) / universe;
        return std::make_pair(p, b);
    }

    // deviates generated are 0-based, i.e. number of failures *between* samples
    // we thus have to add one to every element to get the sample positions
    template <bool addone, typename It>
    static void inplace_prefix_sum(It begin, It end) {
        using value_type = typename std::iterator_traits<It>::value_type;

        if (begin == end)
            return;
        value_type sum = *begin;

        while (++begin != end) {
            sum += *begin;
            // addone is a template parameter, thus no branch is emitted
            if (addone)
                sum++;
            *begin = sum;
        }
    }

    #ifdef SAMPLING_HAVE_SSE2
    // Based on http://stackoverflow.com/a/32501562/3793885 by Peter Cordes
    // In-place prefix sum, optionally incrementing the sum between elements
    // (required for MKL geometric distribution...)
    template <bool addone, typename T>
    static void inplace_prefix_sum_sse(T* data, size_t n) {
        static_assert(std::is_same<T, int>::value, "can only do int");

        __m128i*  datavec   = (__m128i*)data;
        const int vec_elems = sizeof(*datavec) / sizeof(T);

        if (addone)
            data[0]--; // fix first element

        // don't start an iteration beyond this
        const __m128i* endp  = (__m128i*)(data + n - 2 * vec_elems);
        __m128i        carry = _mm_setzero_si128();
        const __m128i  ones  = _mm_set1_epi32(1);

        for (; datavec <= endp; datavec += 2) {
            __m128i x0 = _mm_load_si128(datavec + 0);
            __m128i x1 = _mm_load_si128(datavec + 1); // unroll / pipeline by 1

            if (addone) {
                x0 = _mm_add_epi32(x0, ones);
                x1 = _mm_add_epi32(x1, ones);
            }

            x0 = _mm_add_epi32(x0, _mm_slli_si128(x0, 4)); // sizeof(T)));
            x1 = _mm_add_epi32(x1, _mm_slli_si128(x1, 4)); // sizeof(T)));

            x0 = _mm_add_epi32(x0, _mm_slli_si128(x0, 8)); // 2*sizeof(T)));
            x1 = _mm_add_epi32(x1, _mm_slli_si128(x1, 8)); // 2*sizeof(T)));

            // add the previous sum (carry)
            x0 = _mm_add_epi32(x0, carry);
            // store first to allow destructive shuffle
            _mm_store_si128(datavec, x0);

            // add last element of x0 to all elements of x1 (carry)
            x1 = _mm_add_epi32(_mm_shuffle_epi32(x0, _MM_SHUFFLE(3, 3, 3, 3)), x1);
            _mm_store_si128(datavec + 1, x1);

            // broadcast the high element of x1 as carry for the next iteration
            carry = _mm_shuffle_epi32(x1, _MM_SHUFFLE(3, 3, 3, 3));
        }

        // handle the leftover elements
        T* ptr = (T*)datavec;
        if (ptr < data + n) {
            *ptr += *(ptr - 1);
            inplace_prefix_sum<addone>(ptr, data + n);
        }
    }

    // Dispatch prefix sum to vectorized implementation if possible
    template <bool addone, typename It, typename value_type = typename std::iterator_traits<It>::value_type>
    static typename std::enable_if<std::is_same<value_type, int>::value>::type
    inplace_prefix_sum_disp(It begin, It end) {
        // &(*begin) so that it works with iterators, eww
        inplace_prefix_sum_sse<addone>(&(*begin), end - begin);
    }

    #endif // SAMPLING_HAVE_SSE2

    // Fallback to non-vectorized implementation
    template <bool addone, typename It, typename value_type = typename std::iterator_traits<It>::value_type>
    static typename std::enable_if<!SAMPLING_SSE2 || !std::is_same<value_type, int>::value>::type
    inplace_prefix_sum_disp(It begin, It end) {
        inplace_prefix_sum<addone>(begin, end);
    }

    // unrolled version of inplace_prefix_sum - not really all that helpful
    template <bool addone, size_t unroll = 8, typename It>
    static void inplace_prefix_sum_unroll(It begin, It end) {
        using value_type = typename std::iterator_traits<It>::value_type;
        value_type sum   = 0;

        while (begin + unroll < end) {
            faux_unroll<unroll>::call([&](size_t i) {
                sum += *(begin + i);
                if (addone)
                    ++sum;
                *(begin + i) = sum;
            });
            begin += unroll;
        }

        while (begin < end) {
            sum += *begin;
            *begin++ = sum;
        }
    }

    // Remove 'to_remove' dummies from the range [begin, end) (end is implicit).
    // The holes array specifies the position of the dummies.  This allows for
    // faster compaction than 'std::remove' because we can memmove the region
    // between the holes into place.  Assumes that holes[to_remove] points to
    // the element-beyond-last ('end')
    template <typename It>
    static It compact(It begin, ssize_t* holes, size_t to_remove) {
        using value_type = typename std::iterator_traits<It>::value_type;

        It dest = begin + holes[0];
        // memmove aborts if src == dest, so we don't have to handle it
        for (size_t i = 0; i < to_remove; ++i) {
            assert(holes[i + 1] > holes[i]);
            size_t size = holes[i + 1] - holes[i] - 1;

            memmove(dest, begin + holes[i] + 1, size * sizeof(value_type));
            dest += size;
            assert(begin + holes[i + 1] - 1 == dest + i);
        }
        return dest;
    }

    // Sampling without replacement to pick position of holes.  Optimized for
    // this scenario, not general purpose.
    template <typename It>
    std::unique_ptr<size_t[]> pick_holes(It begin, It end, size_t k, bool sorted) {
        assert(end - begin > (long)k);
        const size_t to_remove = (end - begin) - k;

        // C++11 lacks std::make_unique :(
        std::unique_ptr<size_t[]> holes(new size_t[to_remove + 1]());
        holes[to_remove] = (end - begin);
        size_t hole_idx  = 0;

        const size_t basecase = 1024;
        // Only use for k up to 4MM, it gets slow after that
        if (sorted && k < (1 << 22)) {
            // SORTED hash sampling
            SortedHashSampling<>            hs((ULONG)++seed, to_remove);
            SeqDivideSampling<decltype(hs)> s(hs, basecase, (ULONG)seed);
            // end - begin - 1 because the range is inclusive
            s.sample(end - begin - 1, to_remove, [&](size_t pos) { holes[hole_idx++] = pos; });
        } else {
            HashSampling<>                  hs((ULONG)++seed, to_remove);
            SeqDivideSampling<decltype(hs)> s(hs, basecase, (ULONG)seed);
            // end - begin - 1 because the range is inclusive
            s.sample(end - begin - 1, to_remove, [&](size_t pos) { holes[hole_idx++] = pos; });
            if (sorted)
                std::sort(holes.get(), holes.get() + to_remove);
        }
        assert(hole_idx == to_remove);

        if (sorted) {
            assert(std::is_sorted(holes.get(), holes.get() + to_remove + 1));
        }

        return holes;
    }

    template <typename It>
    It fix_stable(It begin, It end, size_t k) {
        assert(end - begin > (long)k);

        // Get holes (sorted), then apply memmove compactor
        auto holes = pick_holes(begin, end, k, true);

        size_t to_remove = (end - begin) - k;
        return compact(begin, holes.get(), to_remove);
    }

    template <typename It>
    It fix(It begin, It end, size_t k) {
        assert(end - begin > (long)k);
        size_t to_remove = (end - begin) - k;

        // Holes should be sorted so we can process them in one sweep
        auto indices = pick_holes(begin, end, k, true);

        auto    last = end - 1;
        ssize_t pos  = to_remove;
        // handle case where the last element is to be removed
        // revert last postincrement even if loop doesn't match
        if (begin + indices[pos] == last) {
            --last;
            --pos;
        }
        while (pos > 0) {
            // std::iter_swap(begin + indices[pos--], last--);
            *(begin + indices[pos--]) = std::move(*last--);
        }

        return last + 1;
    }

    generator_t rng;
    size_t      seed;
};

} // namespace sampling

#endif // SAMPLING_METHOD_B_HEADER
