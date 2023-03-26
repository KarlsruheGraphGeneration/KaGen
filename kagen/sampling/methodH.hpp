/*******************************************************************************
 * sampling/methodH.hpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef SAMPLING_METHOD_H_HEADER
    #define SAMPLING_METHOD_H_HEADER

    #include "kagen/sampling/definitions.hpp"
    #include "kagen/sampling/rng/select.hpp"
    #include "kagen/tlx/integer_log2.hpp"
    #include "kagen/tlx/likely.hpp"

    #include <algorithm>
    #include <iterator>
    #include <limits>
    #include <vector>

namespace sampling {

template <typename Generator = rng::select_t>
class HashSampling {
public:
    using generator_t = Generator;

    HashSampling(
        ULONG seed, ULONG n, ULONG max_bs_param = (1ULL << 16), ULONG dummy = std::numeric_limits<ULONG>::max())
        : rng(seed),
          max_bs_param(max_bs_param),
          dummy(dummy) {
        // Calculate maximum blocksize and reserve space accordingly
        max_blocksize = std::max(std::min(n, max_bs_param), (ULONG)rng.minimum_reasonable_block_size());
        max_blocksize += (max_blocksize & 0x1); // needs to be even
        randblock.resize(max_blocksize);

        resizeTable(n);
    }

    void resizeTable(ULONG n) {
        // Table size
        table_lg   = 3 + tlx::integer_log2_ceil(n);
        table_size = 1UL << table_lg;
        hash_table.resize(table_size, dummy);
        indices.reserve(table_size);

        // Offset for fast indexing
        offset = &(hash_table[0]);
    }

    // See SH subroutine in Ahrens and Dieter
    template <typename F>
    void sample(ULONG N, ULONG n, F&& callback) {
        ULONG variate, index, hash_elem;
        ULONG population_lg = tlx::integer_log2_ceil(N);
        ULONG address_mask  = (table_lg >= population_lg) ? 0 : population_lg - table_lg;

        // Generate next random block
        ULONG curr_blocksize = std::max(std::min(n, max_bs_param), (ULONG)rng.minimum_reasonable_block_size());
        curr_blocksize += (curr_blocksize & 0x1); // needs to be even
        curr_blocksize = std::min(curr_blocksize, max_blocksize);
        rng.generate_block(randblock, curr_blocksize);
        ULONG array_index = 0;

        while (n > 0) {
            while (true) {
                // Take sample

                // Need a new random block?
                if (array_index >= curr_blocksize) {
                    curr_blocksize = std::max(std::min(n, max_bs_param), (ULONG)rng.minimum_reasonable_block_size());
                    curr_blocksize += (curr_blocksize & 0x1); // needs to be even
                    curr_blocksize = std::min(curr_blocksize, max_blocksize);
                    rng.generate_block(randblock, curr_blocksize);
                    array_index = 0;
                }
                variate = N * randblock[array_index++];

                index     = variate >> address_mask;
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
            callback(variate + 1);
            indices.push_back(index);
            n--;
        }

        clear();
    }

    void clear() {
        for (ULONG index: indices)
            hash_table[index] = dummy;
        indices.clear();
    }

private:
    generator_t rng;

    std::vector<ULONG>  hash_table, indices;
    std::vector<double> randblock;

    ULONG       table_lg, table_size, max_blocksize;
    const ULONG max_bs_param, dummy;
    ULONG*      offset;
};

} // namespace sampling

#endif // SAMPLING_METHOD_H_HEADER
