#pragma once

#include <cstddef>
#include <cstdint>

namespace kagen {
enum class CommDatatype {
    INT,
    UNSIGNED,
    LONG_LONG,
    UNSIGNED_LONG_LONG,
    DOUBLE,
    LONG_DOUBLE,
    BYTE,
    UINT64,
    INT64,
    C_BOOL,
    DERIVED, // Placeholder for derived types; actual size stored in Comm
};

enum class CommOp {
    SUM,
    MAX,
    MIN,
    LOR,
};

// Sentinel for in-place operations (reinterpret_cast is not constexpr)
inline void* const COMM_IN_PLACE = reinterpret_cast<void*>(1);

// Opaque handle for derived datatypes created via TypeContiguous
using CommDerivedType = std::size_t;

inline std::size_t CommDatatypeSize(CommDatatype type) {
    switch (type) {
        case CommDatatype::INT:
            return sizeof(int);
        case CommDatatype::UNSIGNED:
            return sizeof(unsigned);
        case CommDatatype::LONG_LONG:
            return sizeof(long long);
        case CommDatatype::UNSIGNED_LONG_LONG:
            return sizeof(unsigned long long);
        case CommDatatype::DOUBLE:
            return sizeof(double);
        case CommDatatype::LONG_DOUBLE:
            return sizeof(long double);
        case CommDatatype::BYTE:
            return 1;
        case CommDatatype::UINT64:
            return sizeof(std::uint64_t);
        case CommDatatype::INT64:
            return sizeof(std::int64_t);
        case CommDatatype::C_BOOL:
            return sizeof(bool);
        case CommDatatype::DERIVED:
            return 0; // Must use CommDerivedType for size
    }
    return 0;
}
} // namespace kagen
