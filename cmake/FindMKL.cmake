# Adapted from https://github.com/Eyescale/CMake/blob/master/FindMKL.cmake
# licensed under a 3-clause BSD license

# - Try to find the Intel Math Kernel Library
#   Forked from: https://github.com/openmeeg/openmeeg/blob/master/macros/FindMKL.cmake
# Once done this will define
#
# MKL_FOUND - system has MKL
# MKL_ROOT_DIR - path to the MKL base directory
# MKL_INCLUDE_DIR - the MKL include directory
# MKL_LIBRARIES - MKL libraries
#
# There are few sets of libraries:
# Array indexes modes:
# LP - 32 bit indexes of arrays
# ILP - 64 bit indexes of arrays
# Threading:
# SEQUENTIAL - no threading
# INTEL - Intel threading library
# GNU - GNU threading library
# MPI support
# NOMPI - no MPI support

# linux
if(UNIX AND NOT APPLE)
  if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
    set(MKL_ARCH_DIR "intel64_lin")
  else()
    set(MKL_ARCH_DIR "ia32_lin")
  endif()
endif()

# macos
if(APPLE)
  set(MKL_ARCH_DIR "intel64")
endif()

IF(FORCE_BUILD_32BITS)
  set(MKL_ARCH_DIR "ia32")
ENDIF()

if (WIN32)
  if(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
    set(MKL_ARCH_DIR "intel64")
  else()
    set(MKL_ARCH_DIR "ia32")
  endif()
endif()

set (MKL_THREAD_VARIANTS SEQUENTIAL GNUTHREAD INTELTHREAD)
set (MKL_MODE_VARIANTS ILP LP)
set (MKL_MPI_VARIANTS NOMPI)

find_path(MKL_ROOT_DIR
  include/mkl_cblas.h
  PATHS
  $ENV{MKLROOT}
  $ENV{MKLDIR}
  /opt/intel/mkl/*/
  /opt/intel/cmkl/*/
  /Library/Frameworks/Intel_MKL.framework/Versions/Current/lib/universal
  "Program Files (x86)/Intel/ComposerXE-2011/mkl"
  )

# MESSAGE("MKL_ROOT_DIR : ${MKL_ROOT_DIR}") # for debug

find_path(MKL_INCLUDE_DIR
  mkl.h
  PATHS
  ${MKL_ROOT_DIR}/include
  ${INCLUDE_INSTALL_DIR}
  )

find_library(MKL_CORE_LIBRARY
  mkl_core
  PATHS
  ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}
  ${MKL_ROOT_DIR}/lib/
  )

# Threading libraries
find_library(MKL_SEQUENTIAL_LIBRARY
  mkl_sequential
  PATHS
  ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}
  ${MKL_ROOT_DIR}/lib/
  )

find_library(MKL_INTELTHREAD_LIBRARY
  mkl_intel_thread
  PATHS
  ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}
  ${MKL_ROOT_DIR}/lib/
  )

find_library(MKL_GNUTHREAD_LIBRARY
  mkl_gnu_thread
  PATHS
  ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}
  ${MKL_ROOT_DIR}/lib/
  )

# Intel Libraries
IF("${MKL_ARCH_DIR}" STREQUAL "ia32")
  find_library(MKL_LP_LIBRARY
    mkl_intel
    PATHS
    ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}
    ${MKL_ROOT_DIR}/lib/
    )

  find_library(MKL_ILP_LIBRARY
    mkl_intel
    PATHS
    ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}
    ${MKL_ROOT_DIR}/lib/
    )
else()
  find_library(MKL_LP_LIBRARY
    mkl_intel_lp64
    PATHS
    ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}
    ${MKL_ROOT_DIR}/lib/
    )

  find_library(MKL_ILP_LIBRARY
    mkl_intel_ilp64
    PATHS
    ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}
    ${MKL_ROOT_DIR}/lib/
    )
ENDIF()

foreach (MODEVAR ${MKL_MODE_VARIANTS})
  foreach (THREADVAR ${MKL_THREAD_VARIANTS})
    if (MKL_CORE_LIBRARY AND MKL_${MODEVAR}_LIBRARY AND MKL_${THREADVAR}_LIBRARY)
      set(MKL_${MODEVAR}_${THREADVAR}_LIBRARIES
        ${MKL_${MODEVAR}_LIBRARY} ${MKL_${THREADVAR}_LIBRARY} ${MKL_CORE_LIBRARY})
      # message("MKL_${MODEVAR}_${THREADVAR}_LIBRARIES: ${MKL_${MODEVAR}_${THREADVAR}_LIBRARIES}") # for debug
    endif()
  endforeach()
endforeach()

set(MKL_LIBRARIES ${MKL_LP_SEQUENTIAL_LIBRARIES})
#LINK_DIRECTORIES(${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}) # hack

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INCLUDE_DIR MKL_LIBRARIES)

mark_as_advanced(MKL_INCLUDE_DIR MKL_LIBRARIES
  MKL_CORE_LIBRARY MKL_LP_LIBRARY MKL_ILP_LIBRARY
  MKL_SEQUENTIAL_LIBRARY MKL_INTELTHREAD_LIBRARY MKL_GNUTHREAD_LIBRARY
  )
