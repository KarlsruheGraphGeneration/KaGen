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
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/modules")

project(kagen
    VERSION 1.0.3
    DESCRIPTION "Communication-free Massively Distributed Graph Generators"
    HOMEPAGE_URL "https://github.com/KarlsruheGraphGeneration/KaGen")
set(CMAKE_CXX_STANDARD 20)

################################################################################

option(KAGEN_NODEPS "Build KaGen without any dependencies." OFF)

option(KAGEN_USE_MARCH_NATIVE "Compile with -march=native." OFF)
option(KAGEN_USE_CGAL "If available, link against CGal to enable RDG generators." ON)
option(KAGEN_USE_SPARSEHASH "Build with Google Sparsehash. If turned off, fall back to std::unordered_map<>." ON)
option(KAGEN_USE_FAST_MATH "Use -ffast-math." OFF)
option(KAGEN_USE_MKL "Build with Intel MKL random generator." OFF)
option(KAGEN_USE_XXHASH "Build with xxHash. If turned off, path permutation will not be available." ON)
option(KAGEN_USE_BUNDLED_GTEST "Use bundled GoogleTest (downloaded at configure time via FetchContent)." ON)

option(KAGEN_WARNINGS_ARE_ERRORS "Make compiler warnings compiler errors." OFF)

option(KAGEN_BUILD_TESTS "Build unit tests." OFF)
option(KAGEN_BUILD_APPS "Build binaries." ON)
option(KAGEN_BUILD_TOOLS "Build tools." ON)
option(KAGEN_BUILD_EXAMPLES "Build examples." ON)

option(INSTALL_KAGEN "Install KaGen." ON)

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
target_include_directories(morton SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/extlib/libmorton/include>)
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
        if (${CGAL_VERSION} VERSION_GREATER_EQUAL 5)
            list(APPEND KAGEN_LINK_LIBRARIES CGAL::CGAL)
        else()
            include(${CGAL_USE_FILE})
        endif()
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

################################################################################

if (NOT INSTALL_KAGEN)
    return()
endif ()

# If the installation prefix is not set by the user for Unix builds,
# assign a default value to ensure GNUInstallDirs selects the correct paths.
if (UNIX AND CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "/usr/" CACHE STRING "" FORCE)
endif ()

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

install(
    FILES ${PROJECT_SOURCE_DIR}/kagen/kagen.h
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    COMPONENT KaGen_Development
)

set(KAGEN_INSTALL_TARGETS kagen morton)
if (KAGEN_USE_XXHASH)
  list(APPEND KAGEN_INSTALL_TARGETS xxhash)
endif ()
if (KAGEN_BUILD_APPS)
  list(APPEND KAGEN_INSTALL_TARGETS KaGen)
endif ()

install(
    TARGETS ${KAGEN_INSTALL_TARGETS}
    EXPORT KaGenTargets
    RUNTIME COMPONENT KaGen_Runtime
    LIBRARY COMPONENT KaGen_Runtime NAMELINK_COMPONENT KaGen_Development
    ARCHIVE COMPONENT KaGen_Development
    INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

# Allow package maintainers to freely override the path for the configs
set(
    KAGEN_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/KaGen"
    CACHE STRING "CMake package config location relative to the install prefix"
)
set_property(CACHE KAGEN_INSTALL_CMAKEDIR PROPERTY TYPE PATH)
mark_as_advanced(KAGEN_INSTALL_CMAKEDIR)

install(
    FILES
        ${PROJECT_SOURCE_DIR}/cmake/modules/FindMKL.cmake
        ${PROJECT_SOURCE_DIR}/cmake/modules/FindSparsehash.cmake
    DESTINATION ${KAGEN_INSTALL_CMAKEDIR}
    COMPONENT KaGen_Development
)

configure_file(${PROJECT_SOURCE_DIR}/cmake/KaGenConfig.cmake.in KaGenConfig.cmake @ONLY)
install(
    FILES ${PROJECT_BINARY_DIR}/KaGenConfig.cmake
    DESTINATION ${KAGEN_INSTALL_CMAKEDIR}
    COMPONENT KaGen_Development
)

write_basic_package_version_file(
    "KaGenConfigVersion.cmake"
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY SameMajorVersion
)
install(
    FILES ${PROJECT_BINARY_DIR}/KaGenConfigVersion.cmake
    DESTINATION ${KAGEN_INSTALL_CMAKEDIR}
    COMPONENT KaGen_Development
)

install(
    EXPORT KaGenTargets
    NAMESPACE KaGen::
    DESTINATION ${KAGEN_INSTALL_CMAKEDIR}
    COMPONENT KaGen_Development
)
