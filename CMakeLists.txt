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

# prohibit in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif()

# default to Debug building for single-config generators
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message("Defaulting CMAKE_BUILD_TYPE to Release")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type")
endif()

option(KAGEN_USE_LTO "Compile with -flto (link-time optimization)." OFF)
option(KAGEN_USE_CGAL "If available, link against CGal to enable RDG generators." ON)
option(KAGEN_WARNINGS_ARE_ERRORS "Make compiler warnings compiler errors." OFF)
option(KAGEN_USE_FAST_MATH "Use -ffast-math." OFF)

################################################################################

# variables to collect compile-time definitions, include dirs, and libraries
set(KAGEN_DEFINITIONS "")
set(KAGEN_INCLUDE_DIRS 
  ${PROJECT_SOURCE_DIR}/extlib/sampling/sampling
  ${PROJECT_SOURCE_DIR}/extlib/libmorton/include/
  ${PROJECT_SOURCE_DIR}
  )

if(NOT MSVC)
  # enable warnings
  if(CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -W -Wall -Wextra -fPIC -fpermissive")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -W -Wall -Wextra -fPIC -fpermissive")
  endif()

  # enable more warnings
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpedantic")

  # remove unwanted warnings
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")

  # enable -march=native on Debug and Release builds
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-march=native" KAGEN_HAS_MARCH_NATIVE)
  if(KAGEN_HAS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  endif()

  # remove -rdynamic from linker flags (smaller binaries which cannot be loaded
  # with dlopen() -- something no one needs)
  string(REGEX REPLACE "-rdynamic" ""
    CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_C_FLAGS}")
  string(REGEX REPLACE "-rdynamic" ""
    CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS}")

  # warn on conversions
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wconversion -Werror")

  # enable AddressSanitizer
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

  # enable ThreadSanitizer
  if(OFF)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=thread -pie -fPIC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DKAGEN_HAVE_THREAD_SANITIZER=1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread -pie -fPIC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DKAGEN_HAVE_THREAD_SANITIZER=1")
  endif()

  # enable UndefinedBehaviorSanitizer
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined")

  # enable extra warnings on gcc
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wcast-qual")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Winit-self")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Woverloaded-virtual -Wredundant-decls")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wstrict-null-sentinel -Wstrict-overflow=5")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wtautological-compare")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fipa-pure-const -Wsuggest-attribute=const")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsign-conversion")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsign-promo")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wold-style-cast")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsuggest-override")
  endif()
  # enable extra warnings on clang
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wdeprecated")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wabstract-vbase-init")
  endif()

  # Clang < 3.6 0 (?) does not support debug output for auto return types yet.
  # try compiling a platform test for auto return types
  if(ON)
    include(CheckCXXSourceCompiles)
    set(OLD_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS} -g")

    check_cxx_source_compiles(
      "template <typename T> struct A { auto func(int i) { return 42 + i; } };
     int main() { A<int> a; return 0; }"
      KAGEN_CLANG_AUTO_RETURN_DEBUG_INFO)

    if (NOT KAGEN_CLANG_AUTO_RETURN_DEBUG_INFO)
      message(STATUS "compiler does not support -g debug info with auto returns")
      string(REPLACE "-g" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
      string(REPLACE "-g" "" CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
    endif()

    set(CMAKE_REQUIRED_FLAGS "${OLD_CMAKE_REQUIRED_FLAGS}")
  endif()
elseif(MSVC)
  # Force to always compile with W4
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
  # raise warnings as errors
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /WX")

  ### disable verbose warnings:
  # warning C4589: Constructor of abstract class '...' ignores initializer for
  # virtual base class '...' (false positive warnings)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4589")
  # warning C4127: conditional expression is constant
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4127")
  # warning C4458: declaration of '...' hides class member
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4458")
  # warning C4459: declaration of '...' hides global declaration
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4459")
  # warning C4702: unreachable code
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4702")
  # warning C4250: ABC inherits XYZ via dominance
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4250")
  # warning C4503: decorated name length exceeded, name was truncated
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4503")
  # disable lots of warnings about "unsecure" C runtime function
  list(APPEND KAGEN_DEFINITIONS "_CRT_SECURE_NO_WARNINGS")
  # disable "The POSIX name for this item is deprecated. Instead, use the ISO C
  # and C++ conformant name.", Nope. We will not.
  list(APPEND KAGEN_DEFINITIONS "_CRT_NONSTDC_NO_DEPRECATE")
  # disable lots of warnings about "unsecure" STL functions
  list(APPEND KAGEN_DEFINITIONS "_SCL_SECURE_NO_WARNINGS")
  # windef.h bizzarly defines min and max as macros, unless this is defined.
  list(APPEND KAGEN_DEFINITIONS "NOMINMAX")
endif()

if(KAGEN_USE_IWYU)
  # Generate clang compilation database
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

  find_package(PythonInterp)
  find_program(iwyu_tool_path NAMES iwyu_tool.py)
  if (iwyu_tool_path AND PYTHONINTERP_FOUND)
    add_custom_target(iwyu
      ALL      # Remove ALL if you don't iwyu to be run by default.
      COMMAND "${PYTHON_EXECUTABLE}" "${iwyu_tool_path}" -p "${CMAKE_BINARY_DIR}"
      COMMENT "Running include-what-you-use tool"
      VERBATIM)
  endif()
endif()

if(KAGEN_USE_LTO)
  # build with link-time optimization
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-flto CXX_HAS_LTO_FLAG)

  if(CMAKE_BUILD_TYPE MATCHES Release AND CXX_HAS_LTO_FLAG)
    find_program(CMAKE_GCC_AR
      NAMES ${_CMAKE_TOOLCHAIN_PREFIX}gcc-ar${_CMAKE_TOOLCHAIN_SUFFIX}
      HINTS ${_CMAKE_TOOLCHAIN_LOCATION})

    find_program(CMAKE_GCC_NM
      NAMES ${_CMAKE_TOOLCHAIN_PREFIX}gcc-nm
      HINTS ${_CMAKE_TOOLCHAIN_LOCATION})

    find_program(CMAKE_GCC_RANLIB
      NAMES ${_CMAKE_TOOLCHAIN_PREFIX}gcc-ranlib
      HINTS ${_CMAKE_TOOLCHAIN_LOCATION})

    if(CMAKE_GCC_AR AND CMAKE_GCC_NM AND CMAKE_GCC_RANLIB)
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto")
      set(CMAKE_AR "${CMAKE_GCC_AR}")
      set(CMAKE_NM "${CMAKE_GCC_NM}")
      set(CMAKE_RANLIB "${CMAKE_GCC_RANLIB}")
    else()
      message(WARNING "GCC indicates LTO support, but binutils wrappers could not be found. Disabling LTO.")
    endif()
  endif()
endif()

# @todo does not work with Clang 13+ on macOS
if(APPLE)
  # disable warnings about "ranlib: file: libsampling.a(...cpp.o) has no symbols"
  set(CMAKE_C_ARCHIVE_FINISH   "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
  set(CMAKE_CXX_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
endif()

if (KAGEN_WARNINGS_ARE_ERRORS)
    set(CMAKE_CXX_FLAGS "-Werror ${CMAKE_CXX_FLAGS}")
endif()

if (KAGEN_USE_FAST_MATH) 
    set(CMAKE_CXX_FLAGS "-ffast-math ${CMAKE_CXX_FLAGS}") 
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${KAGEN_DEFINITIONS}")
message(STATUS "KaGen CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

###############################################################################
# MPI
###############################################################################
find_package(MPI REQUIRED)
set(KAGEN_LINK_LIBRARIES MPI::MPI_CXX "${KAGEN_LINK_LIBRARIES}")
 
###############################################################################
# CGAL
###############################################################################
if (KAGEN_USE_CGAL)
  set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "Force CGAL to maintain CMAKE flags")
  find_package(CGAL QUIET)
  if(CGAL_FOUND)
    add_definitions(-DKAGEN_CGAL_FOUND) 
    include(${CGAL_USE_FILE})
  else()
    message(STATUS "Could not find CGAL library: Random Delaunay Graphs will not be available")
  endif()
endif()

###############################################################################
# Sampling library 
###############################################################################
# Do not use MKL by default
set(SAMPLING_USE_MKL OFF CACHE BOOL "Use MKL for sampling.")

set(CMAKE_POLICY_DEFAULT_CMP0074 NEW)
add_subdirectory(extlib/sampling)
set(KAGEN_LINK_LIBRARIES sampling spooky "${KAGEN_LINK_LIBRARIES}")

# Do not build other targets from sampling library
set_target_properties(run_methodD 
    run_methodSH 
    run_methodSR
    run_methodP
    run_methodR 
    run_methodH 
    run_methodB
    PROPERTIES EXCLUDE_FROM_ALL 1)

################################################################################

add_subdirectory(kagen)
add_subdirectory(app)
add_subdirectory(library)
add_subdirectory(examples)

################################################################################

add_library(KaGen::KaGen ALIAS kagen)
add_library(KaGen::ckagen ALIAS ckagen)

