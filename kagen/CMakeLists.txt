################################################################################
# include/CMakeLists.txt
#
# CMake file for generator library
#
# Copyright (C) 2017 Sebastian Lamm <lamm@kit.edu>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

# glob general sources
file(GLOB GENERATOR_SRCS
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.[ch]pp
  ${CMAKE_CURRENT_SOURCE_DIR}/generators/barabassi/*.[ch]pp
  ${CMAKE_CURRENT_SOURCE_DIR}/generators/geometric/*.[ch]pp
  ${CMAKE_CURRENT_SOURCE_DIR}/generators/geometric/delaunay/*.[ch]pp
  ${CMAKE_CURRENT_SOURCE_DIR}/generators/geometric/rgg/*.[ch]pp
  ${CMAKE_CURRENT_SOURCE_DIR}/generators/gnm/*.[ch]pp
  ${CMAKE_CURRENT_SOURCE_DIR}/generators/gnp/*.[ch]pp
  ${CMAKE_CURRENT_SOURCE_DIR}/generators/hyperbolic/*.[ch]pp
  # ${CMAKE_CURRENT_SOURCE_DIR}/generators/kronecker/*.[ch]pp
  ${CMAKE_CURRENT_SOURCE_DIR}/io/*.[ch]pp
  ${CMAKE_CURRENT_SOURCE_DIR}/tools/*.[ch]pp
  )

set(GENERATOR_INCS 
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/include/generators
  ${PROJECT_SOURCE_DIR}/include/generators/barabassi
  ${PROJECT_SOURCE_DIR}/include/generators/geometric
  ${PROJECT_SOURCE_DIR}/include/generators/geometric/delaunay
  ${PROJECT_SOURCE_DIR}/include/generators/geometric/rgg
  ${PROJECT_SOURCE_DIR}/include/generators/gnm
  ${PROJECT_SOURCE_DIR}/include/generators/gnp
  ${PROJECT_SOURCE_DIR}/include/generators/hyperbolic
  ${PROJECT_SOURCE_DIR}/include/generators/kronecker
  ${PROJECT_SOURCE_DIR}/include/io
  ${PROJECT_SOURCE_DIR}/include/tools
  )

message("\nLIBS: ${KAGEN_LINK_LIBRARIES}")
add_library(generators STATIC ${GENERATOR_SRCS})
set_target_properties(generators PROPERTIES LINKER_LANGUAGE CXX)
target_compile_definitions(generators PUBLIC ${KAGEN_DEFINITIONS})
target_include_directories(generators PUBLIC ${PROJECT_SOURCE_DIR})
# target_include_directories(generators PUBLIC ${GENERATOR_INCS})
target_include_directories(generators SYSTEM PUBLIC ${KAGEN_INCLUDE_DIRS})
target_link_libraries(generators ${KAGEN_LINK_LIBRARIES})

################################################################################
