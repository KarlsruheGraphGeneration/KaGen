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
if (NOT SAMPLING_HAVE_MKL)
    list(FILTER KAGEN_SOURCE_FILES EXCLUDE REGEX "mkl.*")
endif ()

add_library(kagen OBJECT ${KAGEN_SOURCE_FILES})
target_compile_features(kagen PRIVATE cxx_std_17)
target_link_libraries(kagen PUBLIC ${KAGEN_LINK_LIBRARIES})
target_include_directories(kagen PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(kagen PRIVATE ${KAGEN_INCLUDE_DIRS})
target_compile_options(kagen PRIVATE ${KAGEN_WARNING_FLAGS})

# Allow library users to include kagen.h without kagen/ prefix
target_include_directories(kagen INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS kagen DESTINATION lib)
install(FILES kagen.h DESTINATION include)
