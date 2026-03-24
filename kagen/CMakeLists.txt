################################################################################
# include/CMakeLists.txt
#
# CMake file for generator library
#
# Copyright (C) 2017 Sebastian Lamm <lamm@kit.edu>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

file(GLOB_RECURSE KAGEN_SOURCE_FILES CONFIGURE_DEPENDS 
	*.cpp *.h *.hpp)

# If we do not have CGal, exclude Delaunay generators
if (NOT CGAL_FOUND)
    list(FILTER KAGEN_SOURCE_FILES EXCLUDE REGEX "delaunay.*")
endif ()

# If we do not have MKL, exclude MKL headers
if (NOT MKL_FOUND)
    list(FILTER KAGEN_SOURCE_FILES EXCLUDE REGEX "mkl.*")
endif ()

# In NOMPI mode MPIComm is not compiled; SeqComm handles all single-process communication.
if (KAGEN_NOMPI)
    list(FILTER KAGEN_SOURCE_FILES EXCLUDE REGEX "mpi_comm\\.cpp")
endif ()

add_library(KaGen ${KAGEN_SOURCE_FILES})
target_compile_features(KaGen PRIVATE cxx_std_17)
target_link_libraries(KaGen PUBLIC ${KAGEN_LINK_LIBRARIES})
target_include_directories(KaGen PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)
target_include_directories(KaGen PRIVATE ${KAGEN_INCLUDE_DIRS})
target_compile_options(KaGen PRIVATE ${KAGEN_WARNING_FLAGS})

if (KAGEN_NOMPI)
    target_compile_definitions(KaGen PUBLIC KAGEN_NOMPI)
endif ()

# Allow library users to include kagen.h without kagen/ prefix
target_include_directories(KaGen INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
