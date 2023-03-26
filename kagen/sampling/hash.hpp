/*******************************************************************************
 * sampling/hash.hpp
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 * Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef SAMPLING_HASH_HEADER
#define SAMPLING_HASH_HEADER

#include "kagen/sampling/config.hpp"
#include "kagen/sampling/definitions.hpp"
#include "kagen/sampling/spooky/spooky.h"

#include <cassert>
#include <string>

#ifdef SAMPLING_HAVE_SSE4_2
    #include <smmintrin.h> // crc32 instructions
#endif

#if defined(__GNUC__) && __GNUC__ >= 7
    #define TLX_ATTRIBUTE_FALLTHROUGH __attribute__((fallthrough))
#elif defined(__clang__)
    #define TLX_ATTRIBUTE_FALLTHROUGH [[clang::fallthrough]]
#else
    #define TLX_ATTRIBUTE_FALLTHROUGH
#endif

namespace sampling {

/*!
 * Hashing helper that decides what is hashed
 *
 * Defaults to pointer to the object and sizeof(its type). Override these values
 * for heap-allocated types. Some default overrides are provided.
 */
template <typename T>
struct HashDataSwitch {
    static const char* ptr(const T& x) {
        return reinterpret_cast<const char*>(&x);
    }
    static size_t size(const T&) {
        return sizeof(T);
    }
};

template <>
struct HashDataSwitch<std::string> {
    static const char* ptr(const std::string& s) {
        return s.c_str();
    }
    static size_t size(const std::string& s) {
        return s.length();
    }
};

#ifdef SAMPLING_HAVE_SSE4_2
/**
 * A CRC32C hasher using SSE4.2 intrinsics.
 *
 * Note that you need to provide specializations of HashDataSwitch if you want
 * to hash types with heap storage.
 */
template <typename ValueType>
struct HashCrc32Sse42 {
    // Hash data with Intel's CRC32C instructions
    // Copyright 2008,2009,2010 Massachusetts Institute of Technology.
    // For constant sizes, this is neatly optimized away at higher optimization
    // levels - only a mov (for initialization) and crc32 instructions remain
    static uint32_t hash_bytes(const void* data, size_t length, uint32_t crc = 0xffffffff) {
        const char* p_buf = (const char*)data;
        // The 64-bit crc32 instruction returns a 64-bit value (even though a
        // CRC32 hash has - well - 32 bits. Whatever.
        uint64_t crc_carry = crc;
        for (size_t i = 0; i < length / sizeof(uint64_t); i++) {
            crc_carry = _mm_crc32_u64(crc_carry, *(const uint64_t*)p_buf);
            p_buf += sizeof(uint64_t);
        }
        crc = (uint32_t)crc_carry; // discard the rest
        length &= 7;               // remaining length

        // ugly switch statement, faster than a loop-based version
        switch (length) {
            case 7:
                crc = _mm_crc32_u8(crc, *p_buf++);
                TLX_ATTRIBUTE_FALLTHROUGH;
            case 6:
                crc = _mm_crc32_u16(crc, *(const uint16_t*)p_buf);
                p_buf += 2;
                TLX_ATTRIBUTE_FALLTHROUGH;
            // case 5 is below: 4 + 1
            case 4:
                crc = _mm_crc32_u32(crc, *(const uint32_t*)p_buf);
                break;
            case 3:
                crc = _mm_crc32_u8(crc, *p_buf++);
                TLX_ATTRIBUTE_FALLTHROUGH;
            case 2:
                crc = _mm_crc32_u16(crc, *(const uint16_t*)p_buf);
                break;
            case 5:
                crc = _mm_crc32_u32(crc, *(const uint32_t*)p_buf);
                p_buf += 4;
                TLX_ATTRIBUTE_FALLTHROUGH;
            case 1:
                crc = _mm_crc32_u8(crc, *p_buf);
                break;
            case 0:
                break;
            default: // wat
                assert(false);
        }
        return crc;
    }

    static uint32_t hash(const ValueType& val, uint32_t crc = 0xffffffff) {
        const char* ptr  = HashDataSwitch<ValueType>::ptr(val);
        size_t      size = HashDataSwitch<ValueType>::size(val);
        return hash_bytes(ptr, size, crc);
    }

    uint32_t operator()(const ValueType& val, uint32_t crc = 0xffffffff) const {
        return hash(val, crc);
    }
};
#endif

// CRC32C, adapted from Evan Jones' BSD-licensed implementation at
// http://www.evanjones.ca/crc32c.html
uint32_t crc32_slicing_by_8(uint32_t crc, const void* data, size_t length);

/**
 * Fallback CRC32C implementation in software.
 *
 * Note that you need to provide specializations of HashDataSwitch if you want to
 * hash types with heap storage.
 */
template <typename ValueType>
struct HashCrc32Fallback {
    static uint32_t hash(const ValueType& val, uint32_t crc = 0xffffffff) {
        const char* ptr  = HashDataSwitch<ValueType>::ptr(val);
        size_t      size = HashDataSwitch<ValueType>::size(val);
        return crc32_slicing_by_8(crc, ptr, size);
    }

    uint32_t operator()(const ValueType& val, uint32_t crc = 0xffffffff) const {
        return hash(val, crc);
    }
};

// If SSE4.2 is available, use the hardware implementation, which is roughly
// four to five times faster than the software fallback (less for small sizes).
#ifdef SAMPLING_HAVE_SSE4_2
using CRCHash = HashCrc32Sse42<ULONG>;
#else
using CRCHash = HashCrc32Fallback<ULONG>;
#endif

/******************************************************************************/

/*
//! CityHash wrapper
class CityHash {
public:
    static ULONG hash(ULONG x) {
#ifdef SAMPLING_ENV64BIT
        ULONG hash = CityHash64WithSeed(&x, 8, SAMPLING_SEEDA);
#else
        ULONG hash = CityHash32(&x, 8);
#endif
        return hash;
    };

};
*/

/******************************************************************************/

//! Spooky hash wrapper
class Spooky {
public:
    static ULONG hash(ULONG x) {
#ifdef SAMPLING_ENV64BIT
        ULONG hash = SpookyHash::Hash64(&x, 8, SAMPLING_SEEDA);
#else
        ULONG hash = SpookyHash::Hash32(&x, 4, SAMPLING_SEEDA);
#endif
        return hash;
    }
};

} // namespace sampling

#endif // SAMPLING_HASH_HEADER
