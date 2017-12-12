/*******************************************************************************
 * include/methodSH.hpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef METHOD_SH_HEADER
#define METHOD_SH_HEADER

#include "dSFMT.hpp"
#include "definitions.hpp"

#include <tlx/define.hpp>
#include <tlx/math.hpp>

#include <algorithm>
#include <iterator>
#include <limits>
#include <vector>

namespace sampling {

template <ULONG blocksize = (1ULL << 24),
          ULONG dummy = std::numeric_limits<ULONG>::max()>
class SortedHashSampling {
public:
    SortedHashSampling(ULONG seed, ULONG n) {
        // Modification: dSFMT
        _dSFMT::dsfmt_init_gen_rand(&dsfmt, seed);
        max_blocksize = std::max(std::min(n, blocksize),
                                 (ULONG)_dSFMT::dsfmt_get_min_array_size());
        max_blocksize += (max_blocksize & 0x1); // needs to be even
        randblock.reserve(max_blocksize);

        resizeTable(n);
    }

    void resizeTable(ULONG n) {
        // Table size
        table_lg = 3 + tlx::integer_log2_ceil(n);
        table_size = ipow(2, table_lg);
        hash_table.resize(table_size, dummy);

        // Offset for fast indexing
        offset = &(hash_table[0]);
    }

    // See SH subroutine in Ahrens and Dieter
    template <typename F>
    void sample(ULONG N, ULONG n, F &&callback) {
        ULONG variate, index, hash_elem;
        ULONG population_lg = tlx::integer_log2_ceil(N);
        ULONG address_mask = (table_lg >= population_lg) ? 0 : population_lg - table_lg;
        orig_n = n;

        // Modification: dSFMT
        ULONG curr_blocksize = std::max(std::min(n, blocksize),
                                        (ULONG)_dSFMT::dsfmt_get_min_array_size());
        curr_blocksize += (curr_blocksize & 0x1); // needs to be even
        curr_blocksize = std::min(curr_blocksize, max_blocksize);
        _dSFMT::dsfmt_fill_array_close_open(&dsfmt, &(randblock[0]), curr_blocksize);
        ULONG array_index = 0;
        // Modification: End

        while (n > 0) {
            while (true) {
                // Take sample

                // Modification: dSFMT
                if (array_index >= curr_blocksize) {
                    curr_blocksize = std::max(std::min(n, blocksize),
                                              (ULONG)_dSFMT::dsfmt_get_min_array_size());
                    curr_blocksize += (curr_blocksize & 0x1); // needs to be even
                    curr_blocksize = std::min(curr_blocksize, max_blocksize);
                    _dSFMT::dsfmt_fill_array_close_open(&dsfmt, &(randblock[0]), curr_blocksize);
                    array_index = 0;
                }
                variate = N * randblock[array_index++];
                // Modification: End

                index = variate >> address_mask;
                hash_elem = *(offset + index);

                // Table lookup
                if (TLX_LIKELY(hash_elem == dummy))
                    break; // done
                else if (hash_elem == variate)
                    continue; // already sampled
                else {
                increment:
                    ++index;
                    index &= (table_size - 1);
                    hash_elem = *(offset + index);
                    if (hash_elem == dummy)
                        break; // done
                    else if (hash_elem == variate)
                        continue;   // already sampled
                    goto increment; // keep incrementing
                }
            }
            // Add sample
            *(offset + index) = variate;
            n--;
        }

        // Condense
        ULONG i = 0;
        ULONG j = 0;
        while (i < orig_n) {
            while (*(offset + j) == dummy)
                j++;
            *(offset + i) = *(offset + j);
            if (i != j) *(offset + j) = dummy;
            i++;
            j++;
        }

        // Insertion sort
        ULONG tmp = 0;
        for (i = 1; i < orig_n; i++) {
            tmp = *(offset + i);
            for (j = i; j > 0 && *(offset + j - 1) > tmp; j--)
                *(offset + j) = *(offset + j - 1);
            *(offset + j) = tmp;
        }

        // Output in sorted order and clear
        for (i = 0; i < orig_n; i++) {
            callback(*(offset + i) + 1);
            *(offset + i) = dummy;
        }
    }

    void clear() {
        // Alternative
        // memset(offset, dummy, sizeof(ULONG)*table_size);
        std::fill(hash_table.begin(), hash_table.end(), dummy);
    }

private:
    _dSFMT::dsfmt_t dsfmt;

    std::vector<ULONG> hash_table;
    std::vector<double> randblock;

    ULONG table_lg, table_size, max_blocksize;
    ULONG *offset;
    ULONG orig_n;

    inline ULONG ipow(ULONG base, ULONG exp) {
        ULONG result = 1;
        while (exp) {
            if (exp & 1) result *= base;
            exp >>= 1;
            base *= base;
        }
        return result;
    }
};

} // namespace sampling

#endif // METHOD_SH_HEADER
