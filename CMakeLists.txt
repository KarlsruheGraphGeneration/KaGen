################################################################################
# CMakeLists.txt
#
# Root CMake build script for sampler
#
# Copyright (C) 2017 Lorenz Hübschle-Schneider <lorenz@4z2.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

cmake_minimum_required(VERSION 2.8.12)

project(kagen)
#set_target_properties(kagen PROPERTIES LINKER_LANGUAGE CXX)

# prohibit in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif()

# default to Debug building for single-config generators
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message("Defaulting CMAKE_BUILD_TYPE to Debug")
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type")
endif()

option(KAGEN_USE_LTO "Compile with -flto (link-time optimization)." OFF)

option(KAGEN_OUTPUT_EDGES "Output edge list." OFF)

################################################################################

# variables to collect compile-time definitions, include dirs, and libraries
set(KAGEN_DEFINITIONS "")
set(KAGEN_INCLUDE_DIRS 
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/include/generators
  ${PROJECT_SOURCE_DIR}/include/generators/barabassi
  ${PROJECT_SOURCE_DIR}/include/generators/geometric
  ${PROJECT_SOURCE_DIR}/include/generators/geometric/delaunay
  ${PROJECT_SOURCE_DIR}/include/generators/geometric/rgg
  ${PROJECT_SOURCE_DIR}/include/generators/gnm
  ${PROJECT_SOURCE_DIR}/include/generators/gnp
  ${PROJECT_SOURCE_DIR}/include/generators/hyperbolic
  ${PROJECT_SOURCE_DIR}/include/io
  ${PROJECT_SOURCE_DIR}/include/tools
  ${PROJECT_SOURCE_DIR}/extlib/sampling/include
  )

if(NOT MSVC)
  # enable warnings
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -g -W -Wall -Wextra -fPIC")

  # enable more warnings
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpedantic")

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
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Winit-self -Wnoexcept")
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

if(KAGEN_OUTPUT_EDGES)
  set(KAGEN_DEFINITIONS "${KAGEN_DEFINITIONS} -DOUTPUT_EDGES")
endif()

if(APPLE)
  # disable warnings about "ranlib: file: libsampling.a(...cpp.o) has no symbols"
  set(CMAKE_C_ARCHIVE_FINISH   "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
  set(CMAKE_CXX_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
endif()

message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

###############################################################################

# set_join(VAR "foo" "bar" "abc") sets VAR="foo bar abc"
macro(set_join var)
  set(${var})
  foreach(v ${ARGN})
    set(${var} "${${var}} ${v}")
  endforeach()
  string(STRIP ${var} "${${var}}")
endmacro(set_join)

###############################################################################

# use tlx
add_subdirectory(extlib/sampling)
set(KAGEN_LINK_LIBRARIES sampling ${KAGEN_LINK_LIBRARIES})
set(KAGEN_LINK_LIBRARIES spooky ${KAGEN_LINK_LIBRARIES})

# find MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# find CGAL
set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "Force CGAL to maintain CMAKE flags")
find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})

################################################################################
### Build Macros

# macro for building main programs with correct libraries
macro(do_build TARGETNAME)
  message("\nINCLUDES: ${KAGEN_INCLUDE_DIRS}")
  add_executable(${TARGETNAME} ${ARGN})
  target_include_directories(${TARGETNAME} PUBLIC ${KAGEN_INCLUDE_DIRS})
  target_link_libraries(${TARGETNAME} ${KAGEN_LINK_LIBRARIES})
endmacro(do_build)

# macro for building main programs with correct libraries
macro(mpi_build TARGETNAME)
  add_executable(${TARGETNAME} ${ARGN})
  target_include_directories(${TARGETNAME} PUBLIC ${KAGEN_INCLUDE_DIRS})
  target_link_libraries(${TARGETNAME} ${KAGEN_LINK_LIBRARIES} ${MPI_LIBRARIES})
endmacro(mpi_build)

# macro for building a single-source thrill program
macro(build_prog PROGNAME)
  string(REPLACE "/" "_" TARGETNAME "${PROGNAME}") # replace slashes
  do_build(${TARGETNAME} ${PROGNAME}.cpp ${ARGN})
endmacro(build_prog)

macro(build_mpi_prog PROGNAME)
  string(REPLACE "/" "_" TARGETNAME "${PROGNAME}") # replace slashes
  mpi_build(${TARGETNAME} ${PROGNAME}.cpp ${ARGN})
endmacro(build_mpi_prog)

################################################################################

# descend into library source
# add_subdirectory(include)

# descend into apps
add_subdirectory(app)

################################################################################
