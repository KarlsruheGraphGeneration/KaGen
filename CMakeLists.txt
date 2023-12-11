################################################################################
# CMakeLists.txt
#
# Root CMake build script for generator
#
# Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

cmake_minimum_required(VERSION 3.16)
project(kagen LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

################################################################################

option(KAGEN_NODEPS "Build KaGen without any dependencies." OFF)

option(KAGEN_USE_MARCH_NATIVE "Compile with -march=native." OFF)
option(KAGEN_USE_CGAL "If available, link against CGal to enable RDG generators." ON)
option(KAGEN_USE_SPARSEHASH "Build with Google Sparsehash. If turned off, fall back to std::unordered_map<>." ON)
option(KAGEN_USE_FAST_MATH "Use -ffast-math." OFF)
option(KAGEN_USE_MKL "Build with Intel MKL random generator." OFF)
option(KAGEN_USE_XXHASH "Build with xxHash. If turned off, path permutation will not be available." ON)

option(KAGEN_WARNINGS_ARE_ERRORS "Make compiler warnings compiler errors." OFF)

option(KAGEN_BUILD_TESTS "Build unit tests." OFF)
option(KAGEN_BUILD_APPS "Build binaries." ON)
option(KAGEN_BUILD_EXAMPLES "Build examples." ON)

################################################################################

include(FetchContent) 

FetchContent_Declare(googletest
    SYSTEM 
    GIT_REPOSITORY https://github.com/google/googletest.git 
    GIT_TAG release-1.12.1)

################################################################################

# If KAGEN_NODEPS is set, disable all dependency flags
if (KAGEN_NODEPS)
    message(STATUS "Building without any dependencies.")

    set(KAGEN_USE_CGAL OFF)
    set(KAGEN_USE_SPARSEHASH OFF)
    set(KAGEN_USE_MKL OFF)
    set(KAGEN_USE_XXHASH OFF)
    set(KAGEN_BUILD_TESTS OFF) # requires GoogleTest
endif ()

# Prohibit in-source builds
if ("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(SEND_ERROR "In-source builds are not allowed.")
endif ()

# Default to Release building for single-config generators
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Defaulting CMAKE_BUILD_TYPE to Release")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type")
endif ()

# Warning flags
list(APPEND KAGEN_WARNING_FLAGS
    "-W"
    "-Wall"
    "-Wextra"
    "-Wpedantic"
    "-Wno-unused-local-typedefs"
    )

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    list(APPEND KAGEN_WARNING_FLAGS
        "-Wextra-semi"
        "-fcolor-diagnostics"
        "-Wdeprecated"
        )
endif ()
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    list(APPEND KAGEN_WARNING_FLAGS
        "-Wnoexcept"
        "-Wsuggest-override"
        "-fdiagnostics-color=always"
        "-Wcast-qual"
        "-Winit-self"
        "-Woverloaded-virtual"
        "-Wredundant-decls"
        )
endif ()

if (KAGEN_WARNINGS_ARE_ERRORS)
    list(APPEND KAGEN_WARNING_FLAGS "-Werror")
endif ()

# Enable -march=native on Debug and Release builds
if (KAGEN_USE_MARCH_NATIVE)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-march=native" KAGEN_HAS_MARCH_NATIVE)
    if(KAGEN_HAS_MARCH_NATIVE)
        list(APPEND CMAKE_CXX_FLAGS "-march=native")
    endif ()
endif ()

# Remove -rdynamic from linker flags (smaller binaries which cannot be loaded
# with dlopen() -- something no one needs)
string(REGEX REPLACE "-rdynamic" "" CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_C_FLAGS}")
string(REGEX REPLACE "-rdynamic" "" CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS}")

# Enable UndefinedBehaviorSanitizer
if (OFF)
    list(APPEND CMAKE_CXX_FLAGS "-fsanitize=undefined")
endif ()

# Use -ffast-math
if (KAGEN_USE_FAST_MATH)
    list(APPEND CMAKE_CXX_FLAGS "-ffast-math")
endif ()

###############################################################################
# libmorton
###############################################################################
add_library(morton INTERFACE)
target_include_directories(morton SYSTEM INTERFACE extlib/libmorton/include)
list(APPEND KAGEN_LINK_LIBRARIES morton)

###############################################################################
# xxHash 
###############################################################################
if (KAGEN_USE_XXHASH)
    message(STATUS "Building with xxHash")

    set(${PROJECT_NAME}_XXHASH_DIR "${CMAKE_CURRENT_LIST_DIR}/extlib/xxHash")
    set(BUILD_SHARED_LIBS OFF)
    set(XXHASH_BUILD_ENABLE_INLINE_API ON)
    set(XXHASH_BUILD_XXHSUM OFF)
    add_subdirectory("${${PROJECT_NAME}_XXHASH_DIR}/cmake_unofficial" EXCLUDE_FROM_ALL)
    list(APPEND KAGEN_LINK_LIBRARIES xxHash::xxhash)
    add_definitions(-DKAGEN_XXHASH_FOUND)
endif ()

###############################################################################
# MPI
###############################################################################
set(MPI_DETERMINE_LIBRARY_VERSION TRUE) # needed for KaTestrophe
find_package(MPI REQUIRED)
list(APPEND KAGEN_LINK_LIBRARIES MPI::MPI_CXX)

###############################################################################
# pthread
###############################################################################
find_package(Threads REQUIRED)
list(APPEND KAGEN_LINK_LIBRARIES Threads::Threads)

###############################################################################
# Google Sparsehash
###############################################################################
if (KAGEN_USE_SPARSEHASH)
    message(STATUS "Building with Sparsehash")
    find_package(Sparsehash REQUIRED)
    list(APPEND KAGEN_LINK_LIBRARIES Sparsehash::Sparsehash)
    add_definitions(-DKAGEN_SPARSEHASH_FOUND)
endif ()

###############################################################################
# CGAL
###############################################################################
if (KAGEN_USE_CGAL)
    set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE CACHE BOOL "Do not warn about Debug mode")
    set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "Force CGAL to maintain CMAKE flags")
    find_package(CGAL QUIET)
    if (CGAL_FOUND)
        add_definitions(-DKAGEN_CGAL_FOUND) 
        include(${CGAL_USE_FILE})
    else ()
        message(STATUS "Could not find the CGAL library: Random Delaunay Graphs will not be available")
    endif ()
endif ()

###############################################################################
# Sampling library -> MKL
###############################################################################
if (KAGEN_USE_MKL)
    find_package(MKL)
    if (MKL_FOUND)
        message(STATUS "Building with MKL")
        list(APPEND KAGEN_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
        list(APPEND KAGEN_LINK_LIBRARIES ${MKL_LP_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_SEQUENTIAL_LIBRARY})
        add_definitions(-DSAMPLING_HAVE_MKL)
        add_definitions(-DRMAT_HAVE_MKL)
        add_definitions(-DKAGEN_MKL_FOUND)
    else ()
        message(STATUS "MKL requested but not found, building without MKL")
    endif ()
endif ()

################################################################################

add_subdirectory(kagen)

if (KAGEN_BUILD_APPS)
    add_subdirectory(app)
else ()
    message(STATUS "Apps disabled.")
endif ()

if (KAGEN_BUILD_EXAMPLES) 
    add_subdirectory(examples)
else ()
    message(STATUS "Examples disabled.")
endif ()

################################################################################

add_library(KaGen::KaGen ALIAS kagen)
add_library(KaGen::cKaGen ALIAS kagen) # @deprecated, use KaGen::KaGen

################################################################################

if (KAGEN_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()

