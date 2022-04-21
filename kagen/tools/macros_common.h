/*******************************************************************************
 * include/tools/macros_common.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _MACROS_COMMON_H_
#define _MACROS_COMMON_H_

namespace kagen {

// If DEBUG has not been defined yet then define it as false.
#ifndef DEBUG
    #define DEBUG 0
#endif

#if DEBUG
    #define VLOG(x)  x
    #define VVLOG(x) x
#else
    #define VLOG(x)
    #define VVLOG(x)
#endif

// A macro to define empty copy constructors and assignment operators.
// Use this after private: to hide them to the outside.
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
    TypeName(const TypeName&);             \
    void operator=(const TypeName&)

// A macro to disallow the assignment operators.
#define DISALLOW_ASSIGN(TypeName) void operator=(const TypeName&)

// Helper macro for STR.
#define ASSERT_H_XSTR(x) (#x)

// This macro allows to convert an expression to a string.
#define STR(x) ASSERT_H_XSTR(x)

} // namespace kagen
#endif
